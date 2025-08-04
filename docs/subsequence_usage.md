# Custom Subsequence Extraction Script

This script extends the functionality of `seqtk subseq` with advanced YAML-based filtering capabilities for gene sequence processing.

## Features

1. **Standard subsequence extraction**: Like `seqtk subseq`, extracts specific sequences from a FASTA file using a filter list
   - **Flexible matching**: Sequence IDs can be substrings of filter list entries (e.g., `seq1` matches `sample_seq1_2024`)
2. **Complete sequence removal**: Replace entire sequences with gaps (`-`) for specified sequence IDs
3. **Pattern-based chunk removal**: Remove leading or trailing regions based on pattern matching

## Usage

```bash
subsequence.py --sequence input.fasta --gene pol --yml config.yml [--filter-list ids.txt] [--output output.fasta]
```

### Arguments

- `--sequence, -s`: Input FASTA file containing sequences (required)
- `--gene, -g`: Gene name to look up in YAML configuration (required)
- `--yml, -y`: YAML configuration file with filtering rules (required)
- `--filter-list, -f`: Optional filter list file (one sequence ID per line)
- `--output, -o`: Output FASTA file (default: stdout)

## YAML Configuration Format

```yaml
gene_name:
  to_remove:
    - sequence_id_1
    - sequence_id_2
  to_remove_chunk:
    pattern: "NNN"
    leading:
      - sequence_id_3
      - sequence_id_4
    trailing:
      - sequence_id_5
      - sequence_id_6
```

### Configuration Sections

#### `to_remove`

- Lists sequence IDs that should have their entire sequence replaced with gaps (`-`)
- Useful for removing contaminated or low-quality sequences while maintaining alignment structure

#### `to_remove_chunk`

- **`pattern`**: The nucleotide pattern to search for (exact string match, not regex)
- **`leading`**: Sequence IDs where nucleotides from the start up to and including the first pattern occurrence should be masked
- **`trailing`**: Sequence IDs where nucleotides from the last pattern occurrence to the end should be masked

## Examples

### Example 1: Basic usage with filter list

```bash
subsequence.py -s alignment.fasta -g pol -y config.yml -f sequence_ids.txt -o output.fasta
```

### Example 2: Process all sequences with YAML filtering only

```bash
subsequence.py -s alignment.fasta -g np -y config.yml -o output.fasta
```

### Example 3: Substring matching with filter list

If your filter list contains longer identifiers and your FASTA sequences have shorter IDs:

Filter list (`ids.txt`):

```
experiment_seq1_batch_A_2024
study_seq2_control_group
analysis_seq3_treatment_final
```

FASTA sequences:

```
>seq1
ATCGATCG...
>seq2
GCTAGCTA...
>seq3
TTAATTAA...
```

The script will match `seq1` → `experiment_seq1_batch_A_2024`, `seq2` → `study_seq2_control_group`, etc.

```bash
subsequence.py -s alignment.fasta -g pol -y config.yml -f ids.txt -o output.fasta
```

### Example Configuration

```yaml
pol:
  to_remove:
    - old_reference_1
    - contaminated_seq_2
  to_remove_chunk:
    pattern: "NNNNN"
    leading:
      - sample_A1 # Remove sequence from start to first "NNNNN"
      - sample_B2
    trailing:
      - sample_C3 # Remove sequence from last "NNNNN" to end
      - sample_D4

z:
  to_remove:
    - low_quality_seq
  to_remove_chunk:
    pattern: "NNN"
    leading:
      - partial_seq_1
    trailing:
      - partial_seq_2
```

## Integration with Nextflow

The script is designed to work with the custom Nextflow module `SUBSEQ_CUSTOM`. See `modules/local/subseq/main.nf` for integration details.

## Dependencies

- Python ≥ 3.8
- BioPython ≥ 1.79
- PyYAML ≥ 6.0

These are specified in `modules/local/subseq/environment.yml`.
