#!/usr/bin/env python3
"""
Collect and process input BAMs for PTCP workflows.

Validates BAMs, normalizes sample names, copies single BAMs, merges multiple BAMs
per sample, and writes ordered sample names and chunk BAM paths.
"""

__version__ = "1.0.0"

import argparse
import logging
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import List, Sequence, Union

try:
    from pbcore.io import BamReader  # type: ignore
except ImportError:  # used for test_collect_inputs.py
    BamReader = None

from utils import normalize_bam_name

logger = logging.getLogger(__name__)


def validate_bam_file(bam_file: Path) -> None:
    """Run a quickcheck on the BAM file using samtools."""
    command = ["samtools", "quickcheck", "-u", "-vvvvv", str(bam_file)]
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError as exc:
        raise RuntimeError("samtools is not installed or not in PATH.") from exc

    if result.returncode != 0:
        raise RuntimeError(
            f"samtools quickcheck failed for {bam_file} with error:\n"
            f"{result.stderr.strip()}"
        )


def collect_bam_files(
    hifi_bams: Sequence[str],
    fail_bams: Sequence[str],
    nproc: int,
    biosamples_out: Union[Path, str],
    chunk_bams_out: Union[Path, str],
) -> List[str]:
    """
    Validate, normalize, and copy/merge BAM files.

    Returns the ordered list of sample names written to biosamples_out.
    """
    if BamReader is None:
        raise ImportError("pbcore.io is required to process BAM files.")

    bam_map = defaultdict(list)
    missing_kinetics = []

    for bam_file_str in list(hifi_bams) + list(fail_bams):
        if not bam_file_str:
            continue
        bam_file = Path(bam_file_str)
        validate_bam_file(bam_file)

        bam_name = normalize_bam_name(bam_file)

        sample_names = set()
        bam = BamReader(bam_file)
        for rg in bam.readGroupTable:
            if len(rg.StringID.split("/")) < 2:
                logger.warning(
                    "The BAM file: %s possibly has not been through demultiplexing.",
                    bam_file,
                )
            sample_names.add(rg.SampleName)
        if len(sample_names) > 1:
            samples_str = ", ".join(sorted(sample_names))
            logger.warning(
                "Multiple samples found in a single BAM file: %s -> %s.",
                bam_file,
                samples_str,
            )

        programs = {prog["ID"] for prog in bam.peer.header.get("PG", [])}
        if "jasmine" not in programs:
            if not bam.hasBaseFeature("Ipd"):
                logger.warning(
                    "Warning: Analysis output will not contain methylation calls due to "
                    "missing kinetics in %s.",
                    bam_name,
                )
                missing_kinetics.append(bam_name)
            else:
                logger.warning("Warning: Jasmine has not been run on %s.", bam_name)
        bam_map[bam_name].append(bam_file)

    all_sample_names: List[str] = []
    chunk_bams: List[str] = []
    for i, (sample, bam_files) in enumerate(bam_map.items(), start=1):
        bam_name = f"sample{i:04d}_{sample}.unmapped.bam"
        if len(bam_files) == 1:
            logger.info("Found a single bam for sample %s", sample)
            shutil.copyfile(bam_files[0], bam_name)
        else:
            logger.info("Merging %d BAMs for sample %s", len(bam_files), sample)
            args = ["pbmerge", "-j", str(nproc), "-o", bam_name] + [
                str(b) for b in bam_files
            ]
            subprocess.check_call(args)
        all_sample_names.append(sample)
        chunk_bams.append(bam_name)

    with open(biosamples_out, "wt") as sample_out:
        sample_out.write("\n".join(all_sample_names))
    logger.info("Wrote list of sample names to %s", biosamples_out)
    with open(chunk_bams_out, "wt") as chunk_out:
        chunk_out.write("\n".join(chunk_bams))
    logger.info("Wrote list of chunk BAMs to %s", chunk_bams_out)
    return all_sample_names


def parse_bam_arg(value: str) -> List[str]:
    """Split a semicolon-separated string into a list, filtering empties."""
    return [part for part in value.split(";") if part]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Collect and process BAM inputs for PTCP workflows."
    )
    parser.add_argument(
        "--hifi-bams",
        required=True,
        help="Semicolon-separated list of HiFi BAMs",
    )
    parser.add_argument(
        "--fail-bams",
        default="",
        help="Semicolon-separated list of fail BAMs (optional)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to pass to pbmerge when merging",
    )
    parser.add_argument(
        "--output",
        default="biosamples.txt",
        help="Path to write sample names (default: biosamples.txt)",
    )
    parser.add_argument(
        "--chunk-bams",
        default="chunk_bams.txt",
        help="Path to write ordered chunk BAMs (default: chunk_bams.txt)",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO)",
    )

    args = parser.parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    hifi_bams = parse_bam_arg(args.hifi_bams)
    fail_bams = parse_bam_arg(args.fail_bams)

    collect_bam_files(
        hifi_bams=hifi_bams,
        fail_bams=fail_bams,
        nproc=args.threads,
        biosamples_out=Path(args.output),
        chunk_bams_out=Path(args.chunk_bams),
    )


if __name__ == "__main__":
    main()
