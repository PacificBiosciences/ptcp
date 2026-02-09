#!/usr/bin/env python3
"""
Shared helpers for PTCP input file handling.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re
from typing import Dict, List, Optional, Tuple

SKIP_PATTERNS = [
    "unassigned",
    "sequencing_control",
]

HIFI_FAIL_PATTERN = re.compile(r"\.(hifi_reads|fail_reads)\.")
READS_SUFFIX_PATTERN = re.compile(r"\.reads(\.bam)?$")
BAM_SUFFIX_PATTERN = re.compile(r"\.bam$")
FAIL_PATTERN = re.compile(r"fail", re.IGNORECASE)


def find_bam_files(directory: Path, max_depth: Optional[int] = None) -> List[Path]:
    if max_depth is None:
        return list(directory.rglob("*.bam"))

    bam_files: List[Path] = []

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
    """Strip hifi/fail suffixes and .reads/.bam extensions."""
    bam_name = HIFI_FAIL_PATTERN.sub(".", bam_path.name)
    bam_name = READS_SUFFIX_PATTERN.sub("", bam_name)
    bam_name = BAM_SUFFIX_PATTERN.sub("", bam_name)
    return bam_name


def group_bam_files(bam_files: List[Path]) -> Dict[str, List[Path]]:
    """
    Group BAM files by normalized sample name.
    """
    bam_groups: Dict[str, List[Path]] = defaultdict(list)

    for bam_path in bam_files:
        if should_skip_file(bam_path):
            continue

        normalized_name = normalize_bam_name(bam_path)
        bam_groups[normalized_name].append(bam_path)

    return dict(bam_groups)


def select_primary_bam(bam_files: List[Path]) -> str:
    """
    Select the primary BAM file name from a list of BAM files.
    Prefers non-fail files over fail files.
    """
    for bam_file in bam_files:
        if not FAIL_PATTERN.search(bam_file.stem):
            return bam_file.name
    return bam_files[0].name


def categorize_bam_files(bam_files: List[Path]) -> Tuple[List[str], List[str]]:
    """
    Categorize BAM files into HiFi and fail reads based on filename patterns.
    """
    hifi_reads: List[str] = []
    fail_reads: List[str] = []

    for bam_path in bam_files:
        abs_path = str(bam_path.resolve())
        if FAIL_PATTERN.search(bam_path.name):
            fail_reads.append(abs_path)
        else:
            hifi_reads.append(abs_path)

    return hifi_reads, fail_reads
