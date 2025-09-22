#!/usr/bin/env python3
"""
Create a sample sheet from BAM files for PTCP.

This script processes BAM files in a directory, groups them by sample,
and generates a CSV sample sheet with sample information, note that sex defaults to female (F).
"""

import argparse
from pathlib import Path
from collections import defaultdict
import re
from typing import List, Tuple, Dict, Optional

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


def select_primary_bam(bam_files: List[Path]) -> str:
    """
    Select the primary BAM file name from a list of BAM files.
    Prefers non-fail files over fail files.

    Args:
        bam_files: List of BAM file paths for the same sample

    Returns:
        Name of the primary BAM file
    """
    # First try to find a non-fail file
    for bam_file in bam_files:
        if not FAIL_PATTERN.search(bam_file.stem):
            return bam_file.name
    # If all files contain "fail", return the first one
    return bam_files[0].name


def group_bam_files(bam_files: List[Path]) -> Dict[str, List[Path]]:
    """
    Group BAM files by normalized sample name.

    Args:
        bam_files: List of all BAM file paths

    Returns:
        Dictionary mapping normalized sample names to lists of BAM files
    """
    bam_groups = defaultdict(list)

    for bam_path in bam_files:
        if should_skip_file(bam_path):
            continue

        normalized_name = normalize_bam_name(bam_path)
        bam_groups[normalized_name].append(bam_path)

    return dict(bam_groups)


def validate_bam_groups(bam_groups: Dict[str, List[Path]]) -> None:
    """
    Validate that each sample has at most 2 BAM files.

    Args:
        bam_groups: Dictionary of grouped BAM files

    Raises:
        AssertionError: If any sample has more than 2 BAM files
    """
    for sample_id, bam_files in bam_groups.items():
        if len(bam_files) > 2:
            file_names = [f.name for f in bam_files]
            raise AssertionError(
                f"Found {len(bam_files)} BAM files for sample '{sample_id}': {file_names}. "
                f"Expected at most 2 BAM files per sample (one HiFi and one fail)"
            )


def create_sample_sheet(
    bam_files: List[Path],
) -> Tuple[List[str], Dict[str, List[Path]]]:
    """
    Create a sample sheet from BAM files.

    Args:
        bam_files: List of BAM file paths to process

    Returns:
        Tuple containing:
        - List of CSV lines for the sample sheet
        - Dictionary mapping sample IDs to their BAM files
    """
    bam_groups = group_bam_files(bam_files)
    validate_bam_groups(bam_groups)
    lines = ["bam_name,bam_id,sex"]
    for sample_id, sample_bam_files in bam_groups.items():
        primary_bam_name = select_primary_bam(sample_bam_files)
        lines.append(f"{primary_bam_name},{sample_id},F")
    return lines, bam_groups


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a sample sheet from BAM files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --data /path/to/bam/files
  %(prog)s --data ./data > sample_sheet.csv
        """,
    )
    parser.add_argument(
        "--data", type=Path, help="Directory containing BAM files", required=True
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        help="Maximum recursion depth for finding BAM files (default: unlimited)",
        default=None,
    )

    args = parser.parse_args()

    if not args.data.exists():
        parser.error(f"Directory does not exist: {args.data}")

    if not args.data.is_dir():
        parser.error(f"Path is not a directory: {args.data}")

    bam_files = find_bam_files(args.data, args.max_depth)

    if not bam_files:
        parser.error(f"No BAM files found in directory: {args.data}")

    lines, _ = create_sample_sheet(bam_files)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
