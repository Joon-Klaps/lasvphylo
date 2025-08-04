#!/usr/bin/env python3
"""
Custom subsequence extraction script with advanced filtering capabilities.
Mimics seqtk subseq functionality with additional YAML-based filtering.
"""

import argparse
import sys
import yaml
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Dict, List, Set


def load_yaml_config(yaml_file: Path) -> Dict:
    """Load and parse YAML configuration file."""
    try:
        with open(yaml_file, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
        return config
    except (yaml.YAMLError, IOError, OSError) as e:
        print(f"Error loading YAML file {yaml_file}: {e}", file=sys.stderr)
        sys.exit(1)


def load_filter_list(filter_file: Path) -> Set[str]:
    """Load sequence IDs from filter list file (one ID per line)."""
    try:
        with open(filter_file, "r", encoding="utf-8") as f:
            return {line.strip() for line in f if line.strip()}
    except (IOError, OSError) as e:
        print(f"Error loading filter list {filter_file}: {e}", file=sys.stderr)
        sys.exit(1)


def find_pattern_positions(sequence: str, pattern: str) -> List[int]:
    """Find all positions of a pattern in a sequence (non-overlapping)."""
    positions = []
    start = 0
    while True:
        pos = sequence.find(pattern, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + len(pattern)
    return positions


def mask_leading_regions(sequence: str, pattern: str) -> str:
    """Mask leading regions up to and including the pattern with '-'."""
    positions = find_pattern_positions(sequence, pattern)
    if not positions:
        return sequence

    # Mask from start to end of first pattern occurrence
    mask_end = positions[0] + len(pattern)
    masked_seq = "-" * mask_end + sequence[mask_end:]
    return masked_seq


def mask_trailing_regions(sequence: str, pattern: str) -> str:
    """Mask trailing regions from and including the pattern with '-'."""
    positions = find_pattern_positions(sequence, pattern)
    if not positions:
        return sequence

    # Mask from start of last pattern occurrence to end
    mask_start = positions[-1]
    masked_seq = sequence[:mask_start] + "-" * (len(sequence) - mask_start)
    return masked_seq


def process_sequence(record: SeqRecord, gene_config: Dict, target_ids: Set[str]) -> SeqRecord:
    """Process a single sequence record according to the configuration."""
    seq_id = record.id
    sequence = str(record.seq)

    # Check if this sequence should be filtered by target_ids (equivalent to seqtk filter_list)
    # Allow seq_id to be a substring of any element in target_ids
    # This enables matching short sequence IDs against longer filter list entries
    if target_ids:
        found_match = False
        for target_id in target_ids:
            if seq_id in target_id:
                found_match = True
                break
        if not found_match:
            return None

    # Check if entire sequence should be removed (replaced with gaps)
    to_remove = gene_config.get("to_remove", [])
    if seq_id in to_remove:
        masked_sequence = "-" * len(sequence)
        return SeqRecord(Seq(masked_sequence), id=record.id, description=record.description)

    # Check for chunk removal
    remove_chunk = gene_config.get("to_remove_chunk", {})
    if remove_chunk:
        pattern = remove_chunk.get("pattern", "")
        leading_ids = remove_chunk.get("leading", [])
        trailing_ids = remove_chunk.get("trailing", [])

        if pattern:
            if seq_id in leading_ids:
                sequence = mask_leading_regions(sequence, pattern)
            elif seq_id in trailing_ids:
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

    parser.add_argument("--gene", "-g", type=str, required=True, help="Gene name to look up in YAML configuration")

    parser.add_argument("--yml", "-y", type=Path, help="YAML configuration file with filtering rules")

    parser.add_argument("--filter-list", "-f", type=Path, help="Optional filter list file (like seqtk subseq filter list)")

    parser.add_argument("--output", "-o", type=Path, help="Output FASTA file (default: stdout)")

    args = parser.parse_args()

    # Load configuration
    config = load_yaml_config(args.yml)

    # Get gene-specific configuration
    if args.gene not in config:
        print(f"Warning: Gene '{args.gene}' not found in YAML configuration", file=sys.stderr)
        gene_config = {}
    else:
        gene_config = config[args.gene]

    # Load filter list if provided
    target_ids = set()
    if args.filter_list:
        target_ids = load_filter_list(args.filter_list)

    # Process sequences
    try:
        with open(args.sequence, "r", encoding="utf-8") as input_handle:
            sequences = SeqIO.parse(input_handle, "fasta")

            processed_sequences = []
            for record in sequences:
                processed_record = process_sequence(record, gene_config, target_ids)
                if processed_record is not None:
                    processed_sequences.append(processed_record)

            # Write output
            if args.output:
                with open(args.output, "w", encoding="utf-8") as output_handle:
                    SeqIO.write(processed_sequences, output_handle, "fasta")
            else:
                SeqIO.write(processed_sequences, sys.stdout, "fasta")

    except (IOError, OSError) as e:
        print(f"Error processing sequences: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
