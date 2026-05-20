#!/usr/bin/env python3
"""
GenBank Processor — sort records by taxonomy and extract CDS features to FASTA.

This module provides a command-line tool for batch processing of GenBank (.gb)
files. Three operations are supported:

    sort     Reorder records alphabetically by their taxonomy lineage.
    extract  Emit each CDS feature as a sanitized FASTA entry.
    all      Sort followed by extract (default).

Inputs may be individual files, directories (scanned non-recursively), or both.
With no arguments, the current working directory is scanned after the user
confirms.

Author:   Arshia Farajollahi
License:  MIT
Requires: Python 3.9+, Biopython
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path
from typing import Iterable

try:
    from Bio import SeqIO
except ImportError:
    sys.stderr.write(
        "Error: Biopython is required but not installed.\n"
        "Install it with: pip install biopython\n"
    )
    sys.exit(1)


__version__ = "1.0"

# Output directories created adjacent to each processed input file.
ARCHIVE_DIR_NAME = "processed_gb"
EXTRACT_DIR_NAME = "cds_extracts_from_gb"

# Header-cleaning pattern: drops everything between the first ')' and the
# subsequent '@', collapsing the duplicate accession that GenBank attaches
# to many description fields.
_HEADER_ACCESSION_PATTERN = re.compile(r"\).*@")


# ---------------------------------------------------------------------------
# Core processing
# ---------------------------------------------------------------------------


def sanitize_header(header: str) -> str:
    """Normalize a FASTA header line.

    Replicates the original sed pipeline from the upstream Schott toolkit:
    drop everything between the first ``)`` and the next ``@`` (collapsing
    the duplicate accession block), strip the ``PREDICTED:`` prefix used by
    NCBI for computationally predicted records, remove parentheses, and
    replace spaces with underscores so the header is a single BLAST-safe
    token with no whitespace.
    """
    header = _HEADER_ACCESSION_PATTERN.sub(") ", header)
    header = header.replace("PREDICTED: ", "")
    header = header.replace("(", "").replace(")", "")
    header = header.replace(" ", "_")
    return header


def sort_records_by_taxonomy(input_path: Path, output_path: Path) -> bool:
    """Sort GenBank records in *input_path* by taxonomy and write to *output_path*.

    Records lacking taxonomy annotations are grouped together under ``Unknown``.
    Returns True on success; on failure, the error is written to stderr and
    False is returned.
    """
    try:
        records = list(SeqIO.parse(str(input_path), "genbank"))
        records.sort(
            key=lambda r: "; ".join(r.annotations.get("taxonomy", ["Unknown"]))
        )
        SeqIO.write(records, str(output_path), "genbank")
        return True
    except Exception as exc:
        sys.stderr.write(f"Error sorting {input_path}: {exc}\n")
        return False


def extract_cds_to_fasta(input_path: Path, output_path: Path) -> bool:
    """Write a FASTA file containing every CDS feature found in *input_path*.

    Each entry's header carries the parent record's description, accession,
    and full taxonomy lineage. A failure on any single feature is logged but
    does not abort processing of the remaining features in the file.
    """
    try:
        with open(output_path, "w", newline="\n") as fasta_out:
            for record in SeqIO.parse(str(input_path), "genbank"):
                taxonomy = "; ".join(record.annotations.get("taxonomy", ["Unknown"]))
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    try:
                        sequence = str(feature.location.extract(record).seq)
                        header = sanitize_header(
                            f">{record.description} @{record.name}"
                        )
                        fasta_out.write(f"{header}\n{sequence}\n")
                    except Exception as feat_exc:
                        sys.stderr.write(
                            f"Warning: failed to extract a CDS feature "
                            f"in {record.id}: {feat_exc}\n"
                        )
        return True
    except Exception as exc:
        sys.stderr.write(f"Error extracting from {input_path}: {exc}\n")
        return False


# ---------------------------------------------------------------------------
# Cleanup helpers
# ---------------------------------------------------------------------------


def _archive_original(source: Path, archive_dir: Path) -> None:
    """Move *source* into *archive_dir*, creating the directory if needed."""
    archive_dir.mkdir(parents=True, exist_ok=True)
    try:
        shutil.move(str(source), str(archive_dir / source.name))
        print(f"  -> Original moved to {archive_dir.name}/")
    except Exception as exc:
        sys.stderr.write(f"Error moving {source}: {exc}\n")


def _delete_original(source: Path) -> None:
    """Remove *source* from disk, logging any error to stderr."""
    try:
        source.unlink()
        print("  -> Original deleted")
    except Exception as exc:
        sys.stderr.write(f"Error deleting {source}: {exc}\n")


# ---------------------------------------------------------------------------
# Per-file orchestration
# ---------------------------------------------------------------------------


def process_file(file_path: Path, mode: str, cleanup: bool) -> bool:
    """Run the requested operation on a single GenBank file.

    Args:
        file_path: Path to a ``.gb`` file.
        mode:      One of ``"extract"``, ``"sort"``, or ``"all"``.
        cleanup:   If True, the original input is archived (extract mode) or
                   deleted (sort and all modes) after successful processing.

    Returns:
        True if processing succeeded, False otherwise.
    """
    parent_dir = file_path.parent
    base_name = file_path.stem
    archive_dir = parent_dir / ARCHIVE_DIR_NAME
    extract_dir = parent_dir / EXTRACT_DIR_NAME

    if mode == "extract":
        extract_dir.mkdir(parents=True, exist_ok=True)
        print(f"Extracting CDS from {file_path.name}")
        fasta_out = extract_dir / f"{base_name}_out.fas"
        success = extract_cds_to_fasta(file_path, fasta_out)
        if success and cleanup:
            _archive_original(file_path, archive_dir)
        return success

    if mode == "sort":
        archive_dir.mkdir(parents=True, exist_ok=True)
        print(f"Sorting GenBank records in {file_path.name}")
        sorted_gb = archive_dir / f"{base_name}_sorted.gb"
        success = sort_records_by_taxonomy(file_path, sorted_gb)
        if success and cleanup:
            _delete_original(file_path)
        return success

    if mode == "all":
        extract_dir.mkdir(parents=True, exist_ok=True)
        archive_dir.mkdir(parents=True, exist_ok=True)
        print(f"Processing {file_path.name}")
        sorted_gb = archive_dir / f"{base_name}_sorted.gb"
        fasta_out = extract_dir / f"{base_name}_sorted.fas"
        success = sort_records_by_taxonomy(
            file_path, sorted_gb
        ) and extract_cds_to_fasta(sorted_gb, fasta_out)
        if success and cleanup:
            _delete_original(file_path)
        return success

    raise ValueError(f"Unknown mode: {mode}")


# ---------------------------------------------------------------------------
# Input resolution
# ---------------------------------------------------------------------------


def collect_gb_files(items: Iterable[str]) -> list[Path]:
    """Resolve a mixed list of file and directory paths into ``.gb`` files.

    Files are matched case-insensitively against the ``.gb`` extension.
    Directories are scanned one level deep. Missing paths and skipped files
    produce warnings on stdout.
    """
    collected: list[Path] = []
    for item in items:
        path = Path(item)
        if path.is_dir():
            matches = sorted(
                f for f in path.iterdir() if f.is_file() and f.suffix.lower() == ".gb"
            )
            if not matches:
                print(f"Warning: no .gb files found in directory: {item}")
            collected.extend(matches)
        elif path.is_file():
            if path.suffix.lower() == ".gb":
                collected.append(path)
            else:
                print(f"Warning: skipping non-GB file: {item}")
        else:
            print(f"Error: path not found: {item}")
    return collected


def scan_current_directory() -> list[Path]:
    """Return every ``.gb`` file directly inside the current working directory."""
    return sorted(
        f for f in Path(".").iterdir() if f.is_file() and f.suffix.lower() == ".gb"
    )


# ---------------------------------------------------------------------------
# CLI plumbing
# ---------------------------------------------------------------------------


def resolve_mode(args: argparse.Namespace) -> str:
    """Determine the processing mode from parsed arguments.

    Exits with an error message if more than one mutually exclusive mode
    flag is supplied.
    """
    selected = [
        name
        for name, enabled in (
            ("extract", args.extract),
            ("sort", args.sort),
            ("all", args.all),
        )
        if enabled
    ]
    if len(selected) > 1:
        sys.exit(
            f"Error: conflicting mode flags supplied "
            f"({', '.join(selected)}). "
            "Choose only one of --extract, --sort, or --all."
        )
    return selected[0] if selected else "all"


def confirm_scan(file_count: int, mode: str) -> bool:
    """Ask the user to confirm processing of auto-discovered files."""
    descriptions = {
        "extract": "CDS extraction",
        "sort": "taxonomy sorting",
        "all": "sorting followed by CDS extraction",
    }
    print(f"Found {file_count} .gb file(s) in the current directory.")
    response = input(f"Proceed with {descriptions[mode]}? (y/n): ").strip().lower()
    return response in {"y", "yes"}


def build_parser() -> argparse.ArgumentParser:
    """Construct and return the command-line argument parser."""
    parser = argparse.ArgumentParser(
        prog="gb_processor",
        description=(
            "Sort GenBank records by taxonomy and/or extract CDS features "
            "into sanitized FASTA files."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  gb_processor.py                   scan the current directory\n"
            "  gb_processor.py -a sample.gb      sort and extract one file\n"
            "  gb_processor.py -s data/          sort every .gb in data/\n"
            "  gb_processor.py -e -d *.gb        extract, then archive originals\n"
        ),
    )
    parser.add_argument(
        "-e",
        "-x",
        "--extract",
        "--xtract",
        action="store_true",
        help="extract CDS features only (produces .fas)",
    )
    parser.add_argument(
        "-s",
        "--sort",
        action="store_true",
        help="sort records by taxonomy only (produces .gb)",
    )
    parser.add_argument(
        "-a",
        "--all",
        action="store_true",
        help="sort followed by extract (default)",
    )
    parser.add_argument(
        "-d",
        "-r",
        "--delete",
        "--remove",
        action="store_true",
        help="archive or delete originals after successful processing",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        metavar="PATH",
        help="one or more .gb files and/or directories to process",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Program entry point.

    Args:
        argv: Optional argument list (used primarily for testing). Defaults
              to ``sys.argv[1:]``.

    Returns:
        A process exit code: 0 on success, 1 on setup failure, 2 if any file
        failed to process.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    invoked_without_arguments = len(argv) == 0

    parser = build_parser()
    args = parser.parse_args(argv)
    mode = resolve_mode(args)

    # Resolve the list of files to process: explicit inputs first, otherwise
    # fall back to scanning the current directory.
    if args.inputs:
        input_files = collect_gb_files(args.inputs)
    else:
        if invoked_without_arguments:
            print("No arguments provided. Scanning current directory for .gb files...")
        input_files = scan_current_directory()
        if not input_files:
            sys.stderr.write("Error: no .gb files found in current directory.\n")
            return 1
        if invoked_without_arguments and not confirm_scan(len(input_files), mode):
            print("Operation cancelled by user.")
            return 0

    if not input_files:
        sys.stderr.write("Error: no valid .gb files to process.\n")
        return 1

    failures = 0
    for file_path in input_files:
        if not process_file(file_path, mode, args.delete):
            failures += 1

    total = len(input_files)
    succeeded = total - failures
    print(f"Done. Processed {succeeded}/{total} file(s) successfully.")
    return 0 if failures == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
