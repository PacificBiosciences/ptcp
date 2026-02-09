#!/usr/bin/env python3
"""
Create a sample sheet from BAM files for PTCP.

This script processes BAM files in a directory, groups them by sample,
and generates a CSV sample sheet with sample information, note that sex defaults to female (F).
By default it emits a two-column sheet (bam_id, sex); a legacy three-column sheet
including bam_name can be requested via CLI flag.
"""

__version__ = "1.0.0"

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

from utils import find_bam_files, group_bam_files, select_primary_bam


def validate_bam_groups(bam_groups: Dict[str, List[Path]]) -> None:
    """
    Validate that each sample has at most 2 BAM files.
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
    include_bam_name: bool = False,
) -> Tuple[List[str], Dict[str, List[Path]]]:
    """
    Create a sample sheet from BAM files.
    """
    bam_groups = group_bam_files(bam_files)
    validate_bam_groups(bam_groups)
    header = "bam_name,bam_id,sex" if include_bam_name else "bam_id,sex"
    lines = [header]
    for sample_id, sample_bam_files in bam_groups.items():
        if include_bam_name:
            primary_bam_name = select_primary_bam(sample_bam_files)
            lines.append(f"{primary_bam_name},{sample_id},F")
        else:
            lines.append(f"{sample_id},F")
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
    parser.add_argument(
        "--include-bam-name",
        action="store_true",
        help="Include bam_name column (legacy 3-column output)",
    )

    args = parser.parse_args()

    if not args.data.exists():
        parser.error(f"Directory does not exist: {args.data}")

    if not args.data.is_dir():
        parser.error(f"Path is not a directory: {args.data}")

    bam_files = find_bam_files(args.data, args.max_depth)

    if not bam_files:
        parser.error(f"No BAM files found in directory: {args.data}")

    lines, _ = create_sample_sheet(
        bam_files,
        include_bam_name=args.include_bam_name,
    )
    print("\n".join(lines))


if __name__ == "__main__":
    main()
