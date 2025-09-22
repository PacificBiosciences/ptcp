#!/usr/bin/env python3
"""
Create input JSON files from BAM files for PTCP workflow.

This script processes BAM files in a directory, groups them by sample,
and generates JSON input files with BAM file paths organized into chunks.
"""

import argparse
import json
import sys
import logging
from pathlib import Path
from collections import defaultdict
import re
from typing import List, Dict, Tuple, Optional

SKIP_PATTERNS = [
    "unassigned",
    "sequencing_control",
]

HIFI_FAIL_PATTERN = re.compile(r"\.(hifi_reads|fail_reads)\.")
READS_PATTERN = re.compile(r"\.reads$")
FAIL_PATTERN = re.compile(r"fail", re.IGNORECASE)


def find_bam_files(directory: Path, max_depth: Optional[int] = None) -> List[Path]:
    if max_depth is None:
        return list(directory.rglob("*.bam"))

    bam_files = []

    def _search_recursive(path: Path, current_depth: int) -> None:
        if current_depth > max_depth:
            return

        for item in path.iterdir():
            if item.is_file() and item.suffix == ".bam":
                bam_files.append(item)
            elif item.is_dir() and current_depth < max_depth:
                _search_recursive(item, current_depth + 1)

    _search_recursive(directory, 0)
    return bam_files


def should_skip_file(file_path: Path) -> bool:
    """Check if a BAM file should be skipped based on skip patterns."""
    name_lower = file_path.name.lower()
    return any(pattern in name_lower for pattern in SKIP_PATTERNS)


def normalize_bam_name(bam_path: Path) -> str:
    """
    Normalize BAM file name by removing HiFi/fail suffixes and .reads extension.

    Args:
        bam_path: Path to the BAM file

    Returns:
        Normalized BAM name without extensions and suffixes
    """
    bam_name = HIFI_FAIL_PATTERN.sub(".", bam_path.name)
    bam_name = READS_PATTERN.sub("", bam_name)
    return bam_name.rstrip(".bam")


def group_bam_files(bam_files: List[Path]) -> Dict[str, List[Path]]:
    """
    Group BAM files by normalized sample name.

    Args:
        bam_files: List of all BAM file paths

    Returns:
        Dictionary mapping normalized sample names to lists of BAM files
    """
    bam_groups = defaultdict(list)

    logging.info(f"Processing {len(bam_files)} BAM files.")

    for bam_path in bam_files:
        if should_skip_file(bam_path):
            logging.debug(
                f"Skipping file due to matching SKIP_PATTERNS: {bam_path.name}"
            )
            continue

        normalized_name = normalize_bam_name(bam_path)
        bam_groups[normalized_name].append(bam_path)

    logging.info(f"Created BAM mapping with {len(bam_groups)} unique sample entries.")
    return dict(bam_groups)


def categorize_bam_files(bam_files: List[Path]) -> Tuple[List[str], List[str]]:
    """
    Categorize BAM files into HiFi and fail reads based on filename patterns.

    Args:
        bam_files: List of BAM file paths

    Returns:
        Tuple containing (hifi_reads, fail_reads) as lists of absolute path strings
    """
    hifi_reads = []
    fail_reads = []

    for bam_path in bam_files:
        abs_path = str(bam_path.resolve())
        if FAIL_PATTERN.search(bam_path.name):
            fail_reads.append(abs_path)
        else:
            hifi_reads.append(abs_path)

    return hifi_reads, fail_reads


def fill_template_chunk(
    template_path: Path,
    sample_sheet: Path,
    bam_groups: Dict[str, List[Path]],
    max_files: int,
    start_sample_idx: int = 0,
) -> Tuple[dict, int]:
    """
    Fill template JSON with BAM file information for a single chunk.

    Args:
        template_path: Path to the JSON template file
        sample_sheet: Path to the sample sheet CSV file
        bam_groups: Dictionary mapping sample names to BAM file lists
        max_files: Maximum number of BAM files per chunk
        start_sample_idx: Starting sample index for this chunk

    Returns:
        Tuple containing (filled_template_dict, next_sample_index)
    """
    logging.debug(f"Filling template chunk starting at sample index {start_sample_idx}")

    with open(template_path) as f:
        template = json.load(f)

    template["ptcp.sample_sheet"] = str(sample_sheet)
    template.setdefault("ptcp.hifi_reads", [])
    template.setdefault("ptcp.fail_reads", [])

    sample_names = sorted(bam_groups.keys())
    total_bam_count = 0
    samples_processed = 0

    for sample_name in sample_names[start_sample_idx:]:
        sample_bam_files = bam_groups[sample_name]
        sample_bam_count = len(sample_bam_files)

        # Check if adding this sample would exceed the file limit
        if total_bam_count + sample_bam_count > max_files:
            break

        # Categorize and add BAM files for this sample
        hifi_reads, fail_reads = categorize_bam_files(sample_bam_files)
        template["ptcp.hifi_reads"].extend(hifi_reads)
        template["ptcp.fail_reads"].extend(fail_reads)

        total_bam_count += sample_bam_count
        samples_processed += 1

    next_sample_idx = start_sample_idx + samples_processed
    logging.info(
        f"Chunk filled with {samples_processed} samples, {total_bam_count} BAM files total"
    )

    return template, next_sample_idx


def validate_chunk_size(template: dict, chunk_number: int) -> None:
    """
    Validate that the chunk size is reasonable and log warnings if needed.

    Args:
        template: Filled template dictionary
        chunk_number: Current chunk number for logging
    """
    total_files = len(template.get("ptcp.hifi_reads", [])) + len(
        template.get("ptcp.fail_reads", [])
    )

    logging.info(f"Chunk {chunk_number} contains {total_files} files")

    if total_files > 500:
        logging.warning(
            f"Chunk {chunk_number} contains {total_files} files, "
            f"which exceeds 500 and may cause processing issues"
        )


def process_chunks(
    template_path: Path,
    sample_sheet: Path,
    bam_groups: Dict[str, List[Path]],
    max_files: int,
) -> None:
    """
    Process all samples in chunks and output JSON to stdout.

    The reason to do this is that with certain WDL back-ends you may have to many files to bind, this is a way to constrain this

    Args:
        template_path: Path to the JSON template file
        sample_sheet: Path to the sample sheet CSV file
        bam_groups: Dictionary mapping sample names to BAM file lists
        max_files: Maximum number of BAM files per chunk
    """
    total_samples = len(bam_groups)
    current_sample_idx = 0
    chunk_count = 0

    logging.info(
        f"Processing {total_samples} samples in chunks of max {max_files} files"
    )

    while current_sample_idx < total_samples:
        chunk_count += 1
        logging.info(
            f"Processing chunk {chunk_count}, starting at sample index {current_sample_idx}"
        )

        filled_template, next_sample_idx = fill_template_chunk(
            template_path, sample_sheet, bam_groups, max_files, current_sample_idx
        )

        # Output the JSON for this chunk
        json.dump(filled_template, sys.stdout, indent=2)
        print()  # Add newline between chunks

        validate_chunk_size(filled_template, chunk_count)

        # Handle edge cases where no progress is made
        if current_sample_idx == next_sample_idx:
            total_files_in_chunk = len(
                filled_template.get("ptcp.hifi_reads", [])
            ) + len(filled_template.get("ptcp.fail_reads", []))

            if total_files_in_chunk > 0:
                logging.error(
                    f"No progress made at sample index {current_sample_idx} "
                    f"but {total_files_in_chunk} files were processed. "
                    f"This suggests a single sample exceeds max_files limit."
                )
                next_sample_idx += 1
            elif current_sample_idx < total_samples:
                logging.warning(
                    f"Sample at index {current_sample_idx} resulted in zero files. "
                    f"Advancing to next sample."
                )
                next_sample_idx += 1

        current_sample_idx = next_sample_idx

    logging.info(
        f"Successfully processed all {total_samples} samples in {chunk_count} chunk(s)"
    )


def validate_arguments(args: argparse.Namespace) -> None:
    """
    Validate command line arguments.

    Args:
        args: Parsed command line arguments

    Raises:
        SystemExit: If validation fails
    """
    if not args.data.exists():
        logging.error(f"Data directory does not exist: {args.data}")
        sys.exit(1)

    if not args.data.is_dir():
        logging.error(f"Data path is not a directory: {args.data}")
        sys.exit(1)

    if not args.sample_sheet.exists():
        logging.error(f"Sample sheet does not exist: {args.sample_sheet}")
        sys.exit(1)

    if not args.template.exists():
        logging.error(f"Template file does not exist: {args.template}")
        sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create input JSON file from BAM files for the PTCP workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --data /path/to/bams/ --sample_sheet samples.csv --template template.json > main.inputs.json
  %(prog)s --data /path/to/bams/ --sample_sheet samples.csv --template template.json --max_files 500 --verbose DEBUG
        """,
    )

    parser.add_argument(
        "--data",
        type=Path,
        help="Directory containing BAM files",
        required=True,
    )
    parser.add_argument(
        "--sample_sheet",
        type=Path,
        help="Path to the sample sheet CSV file",
        required=True,
    )
    parser.add_argument(
        "--template",
        type=Path,
        help="Path to template JSON file",
        required=True,
    )
    parser.add_argument(
        "--max-files",
        type=int,
        help="Maximum number of BAM files per chunk (default: %(default)s)",
        default=1000000,
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        help="Maximum recursion depth for finding BAM files (default: unlimited)",
        default=None,
    )
    parser.add_argument(
        "--verbose",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging level (default: %(default)s)",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.verbose.upper()),
        format="%(asctime)s - %(levelname)s - %(module)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logging.debug(f"Arguments: {args}")

    validate_arguments(args)

    bam_files = find_bam_files(args.data, args.max_depth)
    if not bam_files:
        logging.warning(f"No BAM files found in {args.data}")
        sys.exit(0)

    logging.info(f"Found {len(bam_files)} BAM files in {args.data}")

    bam_groups = group_bam_files(bam_files)
    if not bam_groups:
        logging.warning("No valid samples found after filtering")
        sys.exit(0)

    process_chunks(args.template, args.sample_sheet, bam_groups, args.max_files)


if __name__ == "__main__":
    main()
