#!/usr/bin/env python3
"""
Custom subsequence extraction script with advanced filtering capabilities.
Mimics seqtk subseq functionality with additional YAML-based filtering.
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Dict, List, Optional, Set


def load_yaml_config(yaml_file: Path) -> Dict:
    """Return parsed YAML config or an empty dict when the file is empty."""
    try:
        with open(yaml_file, "r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle) or {}
    except (yaml.YAMLError, OSError) as exc:
        logging.error("Failed to read YAML config '%s': %s", yaml_file, exc)
        print(f"Error loading YAML file {yaml_file}: {exc}", file=sys.stderr)
        sys.exit(1)
    return config


def load_filter_list(filter_file: Path) -> Set[str]:
    """Load sequence IDs from a newline-delimited list."""
    try:
        with open(filter_file, "r", encoding="utf-8") as handle:
            return {line.split(" ", 1)[0] for line in handle if line.strip()}
    except OSError as exc:
        logging.error("Failed to read filter list '%s': %s", filter_file, exc)
        print(f"Error loading filter list {filter_file}: {exc}", file=sys.stderr)
        sys.exit(1)


def find_pattern_positions(sequence: str, pattern: str) -> List[int]:
    """Locate pattern occurrences while ignoring gap characters in the sequence."""
    if not pattern:
        return []

    pattern = pattern.upper()
    seq_upper = sequence.upper()
    positions: List[int] = []
    pattern_len = len(pattern)

    for start in range(len(sequence)):
        if seq_upper[start] == "-":
            continue

        collected = []
        idx = start
        while len(collected) < pattern_len and idx < len(sequence):
            char = seq_upper[idx]
            if char != "-":
                collected.append(char)
            idx += 1

        if len(collected) == pattern_len and "".join(collected) == pattern:
            positions.append(start)

    return positions


def mask_leading_regions(sequence: str, pattern: str) -> str:
    """Mask characters from the start through the first pattern occurrence."""
    positions = find_pattern_positions(sequence, pattern)
    if not positions:
        return sequence

    mask_end = positions[0] + len(pattern)
    return "-" * mask_end + sequence[mask_end:]


def mask_trailing_regions(sequence: str, pattern: str) -> str:
    """Mask characters from the last pattern occurrence to the end."""
    positions = find_pattern_positions(sequence, pattern)
    if not positions:
        return sequence

    mask_start = positions[-1]
    return sequence[:mask_start] + "-" * (len(sequence) - mask_start)


def process_sequence(record: SeqRecord, gene_config: Dict, target_ids: Set[str]) -> Optional[SeqRecord]:
    """Apply simple filtering and masking rules to a sequence record."""
    seq_id = record.id
    sequence = str(record.seq)

    if target_ids and seq_id not in target_ids:
        logging.debug("Skipping %s; not present in filter list", seq_id)
        return None

    if seq_id in gene_config.get("to_remove", []):
        logging.debug("Masking entire sequence for %s", seq_id)
        return None

    chunk_cfg = gene_config.get("to_remove_chunk") or {}
    pattern = chunk_cfg.get("pattern")
    if pattern:
        if seq_id in chunk_cfg.get("leading", []):
            sequence = mask_leading_regions(sequence, pattern)
        elif seq_id in chunk_cfg.get("trailing", []):
            sequence = mask_trailing_regions(sequence, pattern)

    return SeqRecord(Seq(sequence), id=record.id, description=record.description)


def main():
    parser = argparse.ArgumentParser(
        description="Advanced subsequence extraction with YAML-based filtering",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example YAML configuration:
gene_name:
  to_remove:
    - id1
    - id2
  to_remove_chunk:
    pattern: "NNN"
    leading:
      - id3
      - id4
    trailing:
      - id5
      - id6
        """,
    )

    parser.add_argument("--sequence", "-s", type=Path, required=True, help="Input FASTA file containing sequences")

    parser.add_argument("--gene", "-g", type=str, help="Gene name to look up in YAML configuration")

    parser.add_argument("--yml", "-y", type=Path, help="YAML configuration file with filtering rules")

    parser.add_argument("--filter-list", "-f", type=Path, help="Optional filter list file (like seqtk subseq filter list)")

    parser.add_argument("--output", "-o", type=Path, help="Output FASTA file (default: stdout)")

    parser.add_argument("--debug", "-d", action="store_true", help="Enable debug logging output")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )

    logging.info("Starting subsequence processing")
    logging.debug("Arguments: %s", vars(args))

    # Prepare configuration
    gene_config: Dict = {}
    if args.yml:
        if not args.gene:
            parser.error("--gene is required when --yml is provided")
        logging.info("Loading configuration for gene: %s", args.gene)
        config = load_yaml_config(args.yml)
        gene_config = config.get(args.gene, {})
        if not gene_config:
            logging.warning("Gene '%s' not found in YAML configuration", args.gene)
            print(f"Warning: Gene '{args.gene}' not found in YAML configuration", file=sys.stderr)
    else:
        logging.info("No YAML configuration provided; using default subsequence behaviour")

    # Load filter list if provided
    target_ids: Set[str] = set()
    if args.filter_list:
        logging.info("Loading filter list from: %s", args.filter_list)
        target_ids = load_filter_list(args.filter_list)

    # Process sequences
    logging.info("Processing sequences from: %s", args.sequence)
    try:
        with open(args.sequence, "r", encoding="utf-8") as input_handle:
            sequences = SeqIO.parse(input_handle, "fasta")

            processed_sequences: List[SeqRecord] = []
            skipped = 0

            for record in sequences:
                processed_record = process_sequence(record, gene_config, target_ids)
                if processed_record is None:
                    skipped += 1
                    continue
                processed_sequences.append(processed_record)

            # Write output
            if args.output:
                logging.info("Writing %d sequences to: %s", len(processed_sequences), args.output)
                with open(args.output, "w", encoding="utf-8") as output_handle:
                    SeqIO.write(processed_sequences, output_handle, "fasta")
            else:
                logging.info("Writing %d sequences to stdout", len(processed_sequences))
                SeqIO.write(processed_sequences, sys.stdout, "fasta")

            logging.info("Skipped %d sequences due to filters", skipped)

    except (IOError, OSError) as e:
        logging.error("Error processing sequences: %s", e)
        print(f"Error processing sequences: {e}", file=sys.stderr)
        sys.exit(1)

    logging.info("Subsequence processing completed successfully")


if __name__ == "__main__":
    main()
