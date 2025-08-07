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
from typing import Dict, List, Set, Optional


def setup_logging(debug: bool = False, log_file: Optional[str] = None) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if debug else logging.INFO

    # Create formatter
    formatter = logging.Formatter(fmt="%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

    # Get root logger
    logger = logging.getLogger()
    logger.setLevel(level)

    # Clear any existing handlers
    logger.handlers.clear()

    # Console handler
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler if log_file is specified
    if log_file:
        file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)


def load_yaml_config(yaml_file: Path) -> Dict:
    """Load and parse YAML configuration file."""
    logging.debug("Loading YAML configuration from: %s", yaml_file)
    try:
        with open(yaml_file, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
        logging.debug("Successfully loaded YAML config with keys: %s", list(config.keys()) if config else "None")
        return config
    except (yaml.YAMLError, IOError, OSError) as e:
        logging.error("Error loading YAML file %s: %s", yaml_file, e)
        print(f"Error loading YAML file {yaml_file}: {e}", file=sys.stderr)
        sys.exit(1)


def load_filter_list(filter_file: Path) -> Set[str]:
    """Load sequence IDs from filter list file (one ID per line)."""
    logging.debug("Loading filter list from: %s", filter_file)
    try:
        with open(filter_file, "r", encoding="utf-8") as f:
            filter_ids = {line.strip() for line in f if line.strip()}
        logging.debug("Successfully loaded %d filter IDs", len(filter_ids))
        logging.debug("Filter IDs: %s", list(filter_ids)[:10])  # Show first 10 for debugging
        return filter_ids
    except (IOError, OSError) as e:
        logging.error("Error loading filter list %s: %s", filter_file, e)
        print(f"Error loading filter list {filter_file}: {e}", file=sys.stderr)
        sys.exit(1)


def find_pattern_positions(sequence: str, pattern: str) -> List[int]:
    """Find all positions of a pattern in a sequence, accounting for gaps within the pattern."""
    logging.debug("Finding pattern '%s' in sequence of length %d", pattern, len(sequence))
    if not pattern:
        logging.debug("Empty pattern provided, returning empty list")
        return []

    positions = []
    pattern_len = len(pattern)
    sequence_len = len(sequence)

    # Try starting from each position in the sequence
    for start_pos in range(sequence_len):
        # Skip if current position is a gap
        if sequence[start_pos] == "-":
            continue

        # Collect non-gap characters starting from this position
        collected_chars = ""
        current_pos = start_pos

        while len(collected_chars) < pattern_len and current_pos < sequence_len:
            if sequence[current_pos] != "-":
                collected_chars += sequence[current_pos]
            current_pos += 1

        # Check if we collected enough characters and they match the pattern
        if len(collected_chars) == pattern_len and collected_chars.upper() == pattern.upper():
            positions.append(start_pos)

    logging.debug("Found pattern '%s' at positions: %s", pattern, positions)
    return positions


def mask_leading_regions(sequence: str, pattern: str) -> str:
    """Mask leading regions up to and including the pattern with '-'."""
    logging.debug("Masking leading regions for pattern '%s'", pattern)
    positions = find_pattern_positions(sequence, pattern)
    if not positions:
        logging.debug("No pattern positions found, returning original sequence")
        return sequence

    # Mask from start to end of first pattern occurrence
    mask_end = positions[0] + len(pattern)
    masked_seq = "-" * mask_end + sequence[mask_end:]
    logging.debug("Masked %d leading characters", mask_end)
    return masked_seq


def mask_trailing_regions(sequence: str, pattern: str) -> str:
    """Mask trailing regions from and including the pattern with '-'."""
    logging.debug("Masking trailing regions for pattern '%s'", pattern)
    positions = find_pattern_positions(sequence, pattern)
    if not positions:
        logging.debug("No pattern positions found, returning original sequence")
        return sequence

    # Mask from start of last pattern occurrence to end
    mask_start = positions[-1]
    masked_seq = sequence[:mask_start] + "-" * (len(sequence) - mask_start)
    logging.debug("Masked %d trailing characters starting from position %d", len(sequence) - mask_start, mask_start)
    return masked_seq


def process_sequence(record: SeqRecord, gene_config: Dict, target_ids: Set[str]) -> SeqRecord:
    """Process a single sequence record according to the configuration."""
    seq_id = record.id
    sequence = str(record.seq)
    logging.debug("Processing sequence: %s (length: %d)", seq_id, len(sequence))

    # Check if this sequence should be filtered by target_ids (equivalent to seqtk filter_list)
    # Allow seq_id to be a substring of any element in target_ids
    # This enables matching short sequence IDs against longer filter list entries
    if target_ids:
        logging.debug("Checking sequence ID '%s' against %d target IDs", seq_id, len(target_ids))
        found_match = False
        for target_id in target_ids:
            if seq_id == target_id:
                logging.debug("Found match: '%s' in target ID '%s'", seq_id, target_id)
                found_match = True
                break
        if not found_match:
            logging.debug("Sequence '%s' not found in target IDs, filtering out", seq_id)
            return None
        logging.debug("Sequence '%s' passed target ID filter", seq_id)

    # Check if entire sequence should be removed (replaced with gaps)
    to_remove = gene_config.get("to_remove", [])
    logging.debug("Checking if sequence '%s' is in to_remove list", seq_id)
    if seq_id in to_remove:
        logging.debug("Sequence '%s' found in to_remove list, masking entire sequence", seq_id)
        masked_sequence = "-" * len(sequence)
        return SeqRecord(Seq(masked_sequence), id=record.id, description=record.description)

    # Check for chunk removal
    remove_chunk = gene_config.get("to_remove_chunk", {})
    if remove_chunk:
        logging.debug("Processing chunk removal for sequence '%s'", seq_id)
        pattern = remove_chunk.get("pattern", "")
        leading_ids = remove_chunk.get("leading", [])
        trailing_ids = remove_chunk.get("trailing", [])

        if pattern:
            if seq_id in leading_ids:
                logging.debug("Applying leading region masking for sequence '%s'", seq_id)
                sequence = mask_leading_regions(sequence.upper(), pattern.upper())
            elif seq_id in trailing_ids:
                logging.debug("Applying trailing region masking for sequence '%s'", seq_id)
                sequence = mask_trailing_regions(sequence.upper(), pattern.upper())
            else:
                logging.debug("Sequence '%s' not in leading or trailing lists, no chunk masking applied", seq_id)

    logging.debug("Finished processing sequence '%s'", seq_id)
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

    parser.add_argument("--debug", "-d", action="store_true", help="Enable debug logging")

    parser.add_argument("--log-file", "-l", type=Path, help="Write debug logs to this file (in addition to console)")

    args = parser.parse_args()

    # Set up logging
    setup_logging(args.debug, str(args.log_file) if args.log_file else None)

    logging.info("Starting subsequence processing")
    if args.log_file:
        logging.info("Logging to file: %s", args.log_file)
    logging.debug("Arguments: %s", vars(args))

    # Load configuration
    logging.info("Loading configuration for gene: %s", args.gene)
    config = load_yaml_config(args.yml)

    # Get gene-specific configuration
    if args.gene not in config:
        logging.warning("Gene '%s' not found in YAML configuration", args.gene)
        print(f"Warning: Gene '{args.gene}' not found in YAML configuration", file=sys.stderr)
        gene_config = {}
    else:
        gene_config = config[args.gene]
        logging.debug("Gene config for '%s': %s", args.gene, gene_config)

    # Load filter list if provided
    target_ids = set()
    if args.filter_list:
        logging.info("Loading filter list from: %s", args.filter_list)
        target_ids = load_filter_list(args.filter_list)
    else:
        logging.debug("No filter list provided")

    # Process sequences
    logging.info("Processing sequences from: %s", args.sequence)
    try:
        with open(args.sequence, "r", encoding="utf-8") as input_handle:
            sequences = SeqIO.parse(input_handle, "fasta")

            processed_sequences = []
            total_count = 0
            filtered_count = 0

            for record in sequences:
                total_count += 1
                processed_record = process_sequence(record, gene_config, target_ids)
                if processed_record is not None:
                    processed_sequences.append(processed_record)
                else:
                    filtered_count += 1

            logging.info("Processed %d sequences total", total_count)
            logging.info("Filtered out %d sequences", filtered_count)
            logging.info("Kept %d sequences for output", len(processed_sequences))

            # Write output
            if args.output:
                logging.info("Writing output to: %s", args.output)
                with open(args.output, "w", encoding="utf-8") as output_handle:
                    SeqIO.write(processed_sequences, output_handle, "fasta")
                logging.info("Successfully wrote %d sequences to output file", len(processed_sequences))
            else:
                logging.debug("Writing output to stdout")
                SeqIO.write(processed_sequences, sys.stdout, "fasta")

    except (IOError, OSError) as e:
        logging.error("Error processing sequences: %s", e)
        print(f"Error processing sequences: {e}", file=sys.stderr)
        sys.exit(1)

    logging.info("Subsequence processing completed successfully")


if __name__ == "__main__":
    main()
