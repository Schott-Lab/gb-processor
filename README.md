# gb_processor

A command-line tool for batch processing [GenBank](https://www.ncbi.nlm.nih.gov/genbank/) (`.gb`) files. It sorts records by taxonomic lineage and extracts coding sequence (CDS) features into sanitized, BLAST-ready FASTA files.

Built on [Biopython](https://biopython.org/), `gb_processor` is designed for phylogenetics and comparative-genomics workflows where many GenBank files need to be reorganized and reduced to their protein-coding regions in a single pass.

---

## Features

- **Sort** records within each file alphabetically by their full taxonomy lineage, grouping un-annotated records under `Unknown`.
- **Extract** every CDS feature into a FASTA file, with headers carrying the parent record's description, accession, and lineage.
- **Combined mode** that sorts and then extracts in one step (the default).
- **Header sanitization** that collapses GenBank's duplicate accession blocks, strips the `PREDICTED:` prefix, removes parentheses, and converts headers into single whitespace-free tokens safe for downstream tools.
- **Flexible inputs** — accepts individual files, directories (scanned one level deep), or a mix of both. With no arguments, it scans the current directory after confirmation.
- **Safe-by-default cleanup** — originals are kept unless you explicitly opt in to archiving or deletion.
- **Resilient processing** — a failure on one feature or file is logged without aborting the rest of the batch.

---

## Requirements

- Python 3.9 or newer
- [Biopython](https://biopython.org/)

Install the dependency with:

```bash
pip install biopython
```

---

## Installation

No packaging step is required. Download `gb_processor.py` and, optionally, make it executable:

```bash
chmod +x gb_processor.py
```

You can then run it directly (`./gb_processor.py`) or via the interpreter (`python gb_processor.py`).

---

## Usage

```
gb_processor.py [MODE] [-d] [PATH ...]
```

If no `PATH` is given, the current working directory is scanned for `.gb` files. When invoked with no arguments at all, you are asked to confirm before processing begins.

### Modes

Exactly one mode may be selected. If none is given, `--all` is used.

| Flag | Aliases | Description |
| --- | --- | --- |
| `-a` | `--all` | Sort by taxonomy, then extract CDS features. **(default)** |
| `-s` | `--sort` | Sort records by taxonomy only. Produces a `.gb` file. |
| `-e` | `-x`, `--extract`, `--xtract` | Extract CDS features only. Produces a `.fas` file. |

### Options

| Flag | Aliases | Description |
| --- | --- | --- |
| `-d` | `-r`, `--delete`, `--remove` | After successful processing, archive originals (extract mode) or delete them (sort and all modes). |
| `--version` | | Print the version and exit. |
| `-h` | `--help` | Show the help message and exit. |

---

## Examples

Scan the current directory and run the default sort-then-extract workflow:

```bash
gb_processor.py
```

Sort and extract a single file:

```bash
gb_processor.py -a sample.gb
```

Sort every `.gb` file in a directory:

```bash
gb_processor.py -s data/
```

Extract CDS features from all matching files, then archive the originals:

```bash
gb_processor.py -e -d *.gb
```

---

## Output

`gb_processor` writes its results into subdirectories created alongside each processed input file:

| Directory | Contents |
| --- | --- |
| `processed_gb/` | Sorted GenBank files (`*_sorted.gb`) and, when `--delete` is used in extract mode, archived originals. |
| `cds_extracts_from_gb/` | Extracted FASTA files (`*_out.fas` in extract mode, `*_sorted.fas` in combined mode). |

### Cleanup behavior

The `--delete` flag behaves differently per mode so that source data is never lost without a recoverable copy:

- **Extract mode** — the original is *moved* into `processed_gb/` (archived).
- **Sort and all modes** — the original is *deleted*, since a sorted copy of the same records already exists in `processed_gb/`.

Without `--delete`, originals are always left untouched.

---

## Header sanitization

CDS FASTA headers are normalized to produce a single BLAST-safe token. The transformation reproduces the `sed` pipeline from the upstream Schott toolkit:

1. Drop everything between the first `)` and the next `@`, collapsing GenBank's duplicate accession block.
2. Strip the `PREDICTED:` prefix that NCBI applies to computationally predicted records.
3. Remove parentheses.
4. Replace spaces with underscores.

---

## Exit codes

| Code | Meaning |
| --- | --- |
| `0` | Success (or user cancelled at the confirmation prompt). |
| `1` | Setup failure — for example, no `.gb` files were found. |
| `2` | One or more files failed to process. |

---

## License

Released under the [MIT License](https://opensource.org/licenses/MIT).

## Author

Arshia Farajollahi
