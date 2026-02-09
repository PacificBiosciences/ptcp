#!/usr/bin/env python3
"""Convert ptcp-qc JSON outputs (aggregate or per-sample) into CSV files."""

from __future__ import annotations

__version__ = "1.1.2"

import argparse
import csv
import json
import logging
import math
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


def accuracy_to_q(accuracy: float) -> Optional[float]:
    """
    Convert per-base accuracy (0..1) to Phred Q.
    Q = -10 * log10(P_error)
    """
    if accuracy == -1.0:
        return None
    if not (0.0 <= accuracy <= 1.0):
        raise ValueError(f"accuracy must be in [0, 1], got {accuracy!r}")
    if accuracy == 1.0:
        return float("inf")
    p_err = 1.0 - accuracy
    return -10.0 * math.log10(p_err)


AGGREGATE_OUTPUT_COLUMNS = [
    "Sample",
    "total_reads",
    "Median Mapped Read Length",
    "Median Mapped Read Quality",
    "Mean Target Coverage",
    "loci_median_coverage",
    "Percent of Targets >=10-fold Coverage",
    "Percent of Targets >=20-fold Coverage",
    "Percent of Targets >=30-fold Coverage",
    "Percent of Targets with Low Coverage (<5X)",
    "Percent of On-Target Reads",
    "Percent of Duplicate Reads",
    "on_target_total",
    "on_target_hifi",
    "on_target_fail",
    "off_target_total",
    "off_target_hifi",
    "off_target_fail",
    "unmapped_total",
    "unmapped_hifi",
    "unmapped_fail",
    "loci_count",
    "sex_ratios_ratio_x_auto",
    "sex_ratios_ratio_y_auto",
]

TRGT_OUTPUT_COLUMNS = [
    "sample",
    "locus",
    "allele",
    "read_count",
    "consensus_size",
    "min_size",
    "max_size",
    "repeat_unit",
    "motif_counts",
    "motif_spans",
]

PARAPHASE_OUTPUT_COLUMNS = [
    "sample",
    "region",
    "CN",
    "CN_adjusted",
    "inversion",
    "genotype_adjusted",
    "paraphase_sv_calls",
    "sawfish_sv_calls",
]

COVERAGE_OUTPUT_COLUMNS = [
    "sample",
    "region",
    "gene",
    "total reads",
    "hifi reads",
    "median read quality",
    "median read passes",
    "median read length",
]

REQUIRED_HEADER_FIELDS = (
    "report_type",
    "genome_version",
    "targets_bed",
    "ptcp-qc_version",
    "timestamp",
)


def setup_logging(verbosity: int) -> None:
    if verbosity == 0:
        level = logging.WARNING
    elif verbosity == 1:
        level = logging.INFO
    else:
        level = logging.DEBUG

    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )


@dataclass
class CsvTable:
    columns: List[str]
    rows: List[List[str]]


@dataclass(frozen=True)
class MissingPolicy:
    numeric: str = ""
    string: str = ""


@dataclass
class ParsedReport:
    path: Path
    data: Dict[str, Any]
    header: Dict[str, Any]
    report_type: str


def load_report(path: Path) -> ParsedReport:
    with path.open(encoding="utf-8") as handle:
        data = json.load(handle)

    header = extract_header(data, path)
    report_type = determine_report_type(data, header, path)
    header.setdefault("report_type", report_type)

    return ParsedReport(
        path=path,
        data=data,
        header=header,
        report_type=report_type,
    )


def extract_header(data: Any, path: Path) -> Dict[str, Any]:
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object at the top level")

    raw_header = data.get("header")
    if raw_header is None:
        header: Dict[str, Any] = {}
    elif isinstance(raw_header, dict):
        header = dict(raw_header)
    else:
        logger.warning(f"{path} has invalid header format; ignoring header")
        header = {}

    for field in REQUIRED_HEADER_FIELDS:
        if field not in header and field in data:
            header[field] = data[field]

    missing = [field for field in REQUIRED_HEADER_FIELDS if field not in header]
    if missing and missing != ["report_type"]:
        fields = ", ".join(sorted(missing))
        raise ValueError(f"{path} header missing required field(s): {fields}")

    return header


def determine_report_type(
    data: Dict[str, Any], header: Dict[str, Any], path: Path
) -> str:
    if "report_type" in header and header["report_type"]:
        return str(header["report_type"])
    if "report_type" in data and data["report_type"]:
        return str(data["report_type"])
    if "sample_stats" in data:
        return "aggregate"
    if "locus_results" in data:
        return "sample"
    raise ValueError(f"{path} missing 'report_type' and cannot infer report type")


def sample_name_from_report(report: ParsedReport) -> str:
    sample_name = report.data.get("sample_name")
    if not sample_name:
        logger.warning(f"{report.path} missing 'sample_name'; using filename stem")
        return report.path.stem
    return str(sample_name)


def aggregate_stats_by_sample(report: ParsedReport) -> Dict[str, Dict[str, Any]]:
    sample_stats = report.data.get("sample_stats")
    if not isinstance(sample_stats, list):
        raise ValueError(f"{report.path} missing 'sample_stats' for aggregate report")

    by_sample: Dict[str, Dict[str, Any]] = {}
    for entry in sample_stats:
        if not isinstance(entry, dict):
            logger.warning(
                f"{report.path} contains non-dict sample_stats entry; skipping"
            )
            continue
        sample_name = str(entry.get("sample_name", "")).strip()
        if not sample_name:
            logger.warning(f"{report.path} contains sample_stats entry without name")
            continue
        if sample_name in by_sample:
            logger.warning(
                f"{report.path} contains duplicate sample_stats for {sample_name}"
            )
            continue
        by_sample[sample_name] = entry
    return by_sample


def sample_stats_by_sample(
    sample_reports: List[ParsedReport],
) -> Dict[str, Dict[str, Any]]:
    stats_by_sample: Dict[str, Dict[str, Any]] = {}
    for report in sample_reports:
        sample_name = sample_name_from_report(report)
        stats = report.data.get("stats")
        stats_by_sample[sample_name] = stats if isinstance(stats, dict) else {}
    return stats_by_sample


def get_dict_value(container: Any, key: Any) -> Any:
    if not isinstance(container, dict):
        return None
    if key in container:
        return container[key]
    key_str = str(key)
    if key_str in container:
        return container.get(key_str)
    try:
        key_int = int(key_str)
    except (TypeError, ValueError):
        return None
    return container.get(key_int)


def table_from_aggregate_and_sample_reports(
    aggregate_report: ParsedReport,
    sample_reports: List[ParsedReport],
    sample_order: List[str],
    sample_display: Callable[[str], str],
    missing: MissingPolicy,
) -> CsvTable:
    aggregate_by_sample = aggregate_stats_by_sample(aggregate_report)
    sample_stats = sample_stats_by_sample(sample_reports)

    all_samples = sorted(set(aggregate_by_sample) | set(sample_stats))
    sample_order_set = set(sample_order)
    ordered_samples = [s for s in sample_order if s in all_samples] + [
        s for s in all_samples if s not in sample_order_set
    ]

    rows: List[List[str]] = []
    for sample_name in ordered_samples:
        summary = aggregate_by_sample.get(sample_name, {})
        stats = sample_stats.get(sample_name, {})

        coverage = summary.get("coverage") or {}
        read_length_stats = summary.get("read_length_stats") or {}
        read_quality_stats = summary.get("read_quality_stats") or {}

        on_target = stats.get("on_target") or {}
        off_target = stats.get("off_target") or {}
        unmapped = stats.get("unmapped") or {}
        loci = stats.get("loci") or {}
        loci_coverages = loci.get("coverages") or {}
        loci_coverage_stats = loci_coverages.get("stats") or {}
        sex_ratios = stats.get("sex_ratios") or {}

        total_reads = summary.get("total_reads")
        if total_reads is None:
            total_reads = stats.get("total_reads")

        rows.append(
            [
                sample_display(sample_name),
                format_int(total_reads, missing.numeric),
                format_int(read_length_stats.get("median"), missing.numeric),
                format_quality(read_quality_stats.get("median"), missing.string),
                format_float_trim(coverage.get("mean_coverage"), missing.numeric),
                format_float_trim(loci_coverage_stats.get("median"), missing.numeric),
                format_percentage(coverage.get("pct_targets_ge_10"), missing.numeric),
                format_percentage(coverage.get("pct_targets_ge_20"), missing.numeric),
                format_percentage(coverage.get("pct_targets_ge_30"), missing.numeric),
                format_percentage(coverage.get("pct_targets_lt_5"), missing.numeric),
                format_percentage(summary.get("on_target_fraction"), missing.numeric),
                format_percentage(summary.get("duplicate_fraction"), missing.numeric),
                format_int(on_target.get("total"), missing.numeric),
                format_int(on_target.get("hifi"), missing.numeric),
                format_int(on_target.get("fail"), missing.numeric),
                format_int(off_target.get("total"), missing.numeric),
                format_int(off_target.get("hifi"), missing.numeric),
                format_int(off_target.get("fail"), missing.numeric),
                format_int(unmapped.get("total"), missing.numeric),
                format_int(unmapped.get("hifi"), missing.numeric),
                format_int(unmapped.get("fail"), missing.numeric),
                format_int(loci.get("count"), missing.numeric),
                format_float_trim(sex_ratios.get("ratio_x_auto"), missing.numeric),
                format_float_trim(sex_ratios.get("ratio_y_auto"), missing.numeric),
            ]
        )

    return CsvTable(columns=AGGREGATE_OUTPUT_COLUMNS, rows=rows)


def table_from_trgt_reports(
    sample_reports: List[ParsedReport],
    sample_order: List[str],
    sample_display: Callable[[str], str],
    missing: MissingPolicy,
) -> CsvTable:
    trgt_loci: set[str] = set()
    repeat_units: Dict[str, str] = {}
    values: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]] = {}

    for report in sample_reports:
        sample_name = sample_name_from_report(report)
        locus_results = report.data.get("locus_results")
        if not isinstance(locus_results, dict):
            raise ValueError(f"{report.path} missing 'locus_results' for sample report")

        for _group, locus in locus_results.items():
            if not isinstance(locus, dict):
                continue
            trgt_results = locus.get("trgt_results")
            if not isinstance(trgt_results, dict):
                continue
            for trgt_name, trgt in trgt_results.items():
                if not isinstance(trgt, dict):
                    continue
                trgt_loci.add(trgt_name)

                repeat_unit = trgt.get("repeat_unit")
                if repeat_unit is not None and trgt_name not in repeat_units:
                    repeat_units[trgt_name] = str(repeat_unit)

                alleles = trgt.get("alleles") or {}
                read_stats = trgt.get("read_stats") or {}

                sample_store = values.setdefault(sample_name, {})
                locus_store = sample_store.setdefault(trgt_name, {})

                for allele_key in ("0", "1"):
                    allele = get_dict_value(alleles, allele_key) or {}
                    read_stat = get_dict_value(read_stats, allele_key) or {}

                    cov_value = get_nested(read_stat, "coverage")
                    if cov_value is None:
                        cov_value = get_nested(allele, "coverage")

                    locus_store[allele_key] = {
                        "read_count": cov_value,
                        "consensus_size": get_nested(allele, "size"),
                        "min_size": get_nested(get_nested(allele, "ci", {}), "min"),
                        "max_size": get_nested(get_nested(allele, "ci", {}), "max"),
                        "repeat_unit": repeat_unit,
                        "motif_counts": get_nested(allele, "motif_count"),
                        "motif_spans": get_nested(allele, "motif_spans"),
                    }

    rows: List[List[str]] = []
    for locus_name in sorted(trgt_loci):
        global_repeat = repeat_units.get(locus_name)
        for sample_name in sample_order:
            per_sample = values.get(sample_name, {})
            per_locus = per_sample.get(locus_name, {})
            for allele_key in ("0", "1"):
                allele_values = per_locus.get(allele_key, {})
                repeat_unit = allele_values.get("repeat_unit")
                if repeat_unit is None:
                    repeat_unit = global_repeat

                rows.append(
                    [
                        sample_display(sample_name),
                        locus_name,
                        allele_key,
                        format_int(allele_values.get("read_count"), missing.numeric),
                        format_int(
                            allele_values.get("consensus_size"), missing.numeric
                        ),
                        format_int(allele_values.get("min_size"), missing.numeric),
                        format_int(allele_values.get("max_size"), missing.numeric),
                        format_string(repeat_unit, missing.string),
                        format_string(
                            allele_values.get("motif_counts"), missing.string
                        ),
                        format_string(allele_values.get("motif_spans"), missing.string),
                    ]
                )

    return CsvTable(columns=TRGT_OUTPUT_COLUMNS, rows=rows)


def format_sv_breakpoints(sv_breakpoints: Any) -> str:
    if not isinstance(sv_breakpoints, list):
        return ""
    formatted: List[str] = []
    for breakpoint in sv_breakpoints:
        if not isinstance(breakpoint, dict):
            continue
        chrom = breakpoint.get("chrom", "")
        start = breakpoint.get("start", "")
        end = breakpoint.get("end", "")
        svtype = breakpoint.get("svtype", "")

        annotation = breakpoint.get("annotation")
        if isinstance(annotation, dict) and "name" in annotation:
            annotation_name = annotation.get("name", "")
            formatted.append(f"{chrom}:{start}:{end}:{svtype}:{annotation_name}")
        else:
            formatted.append(f"{chrom}:{start}:{end}:{svtype}")

    return ",".join(formatted)


def table_from_paraphase_reports(
    sample_reports: List[ParsedReport],
    sample_order: List[str],
    sample_display: Callable[[str], str],
    missing: MissingPolicy,
) -> CsvTable:
    regions: set[str] = set()
    values: Dict[str, Dict[str, Dict[str, Any]]] = {}

    for report in sample_reports:
        sample_name = sample_name_from_report(report)
        locus_results = report.data.get("locus_results")
        if not isinstance(locus_results, dict):
            raise ValueError(f"{report.path} missing 'locus_results' for sample report")

        for _group, locus in locus_results.items():
            if not isinstance(locus, dict):
                continue
            paraphase_results = locus.get("paraphase_results")
            if not isinstance(paraphase_results, dict):
                continue
            sawfish_results = locus.get("sawfish_results")

            homology_adjustment = paraphase_results.get("homology_adjustment")
            if not isinstance(homology_adjustment, dict):
                homology_adjustment = {}

            for paraphase_name, paraphase in paraphase_results.items():
                if paraphase_name == "homology_adjustment":
                    continue
                if not isinstance(paraphase, dict):
                    continue

                smn_info = paraphase.get("smn_info")
                if isinstance(smn_info, dict):
                    for region in ("smn1", "smn2"):
                        regions.add(region)
                        sample_store = values.setdefault(sample_name, {})
                        sample_store[region] = {
                            "CN": smn_info.get(f"{region}_cn"),
                            "CN_adjusted": homology_adjustment.get(
                                f"adjusted_{region}_cn"
                            ),
                        }
                    continue

                region = str(paraphase_name)
                regions.add(region)

                genotype_adjusted = paraphase.get("genotype_adjusted")
                genotype = (
                    genotype_adjusted
                    if genotype_adjusted is not None
                    else paraphase.get("genotype")
                )

                sv_calls = paraphase.get("sv_calls")
                sv_calls_value = (
                    ",".join(sorted(sv_calls.keys()))
                    if isinstance(sv_calls, dict)
                    else ""
                )

                inversion = ""
                f8_info = paraphase.get("f8_info")
                if isinstance(f8_info, dict) and "has_inversion" in f8_info:
                    inversion = format_bool(
                        f8_info.get("has_inversion"), missing.string
                    )

                sawfish_sv_calls = ""
                if isinstance(sawfish_results, dict):
                    sawfish_entry = sawfish_results.get(paraphase_name)
                    if isinstance(sawfish_entry, dict):
                        sawfish_sv_calls = format_sv_breakpoints(
                            sawfish_entry.get("sv_breakpoints")
                        )

                sample_store = values.setdefault(sample_name, {})
                sample_store[region] = {
                    "CN": paraphase.get("total_cn"),
                    "CN_adjusted": None,
                    "inversion": inversion,
                    "genotype_adjusted": genotype,
                    "paraphase_sv_calls": sv_calls_value,
                    "sawfish_sv_calls": sawfish_sv_calls,
                }

    rows: List[List[str]] = []
    for region in sorted(regions):
        for sample_name in sample_order:
            entry = values.get(sample_name, {}).get(region, {})
            rows.append(
                [
                    sample_display(sample_name),
                    region,
                    format_int(entry.get("CN"), missing.numeric),
                    format_int(entry.get("CN_adjusted"), missing.numeric),
                    format_string(entry.get("inversion"), missing.string),
                    format_string(entry.get("genotype_adjusted"), missing.string),
                    format_string(entry.get("paraphase_sv_calls"), missing.string),
                    format_string(entry.get("sawfish_sv_calls"), missing.string),
                ]
            )

    return CsvTable(columns=PARAPHASE_OUTPUT_COLUMNS, rows=rows)


def read_stats_to_fields(reads: Any) -> Dict[str, Any]:
    if not isinstance(reads, dict):
        return {}
    stats = reads.get("stats")
    if not isinstance(stats, dict):
        stats = {}
    read_quality = stats.get("read_quality")
    read_passes = stats.get("read_passes")
    read_length = stats.get("read_length")
    return {
        "total_reads": reads.get("n"),
        "hifi_reads": reads.get("hifi"),
        "median_read_quality": get_nested(read_quality, "median"),
        "median_read_passes": get_nested(read_passes, "median"),
        "median_read_length": get_nested(read_length, "median"),
    }


def table_from_coverage_reports(
    sample_reports: List[ParsedReport],
    sample_order: List[str],
    sample_display: Callable[[str], str],
    missing: MissingPolicy,
) -> CsvTable:
    keys: set[Tuple[str, str]] = set()
    values: Dict[str, Dict[Tuple[str, str], Dict[str, Any]]] = {}

    for report in sample_reports:
        sample_name = sample_name_from_report(report)
        locus_results = report.data.get("locus_results")
        if not isinstance(locus_results, dict):
            raise ValueError(f"{report.path} missing 'locus_results' for sample report")

        for region, locus in locus_results.items():
            if not isinstance(locus, dict):
                continue
            subgroups = locus.get("subgroups")
            has_subgroups = isinstance(subgroups, dict) and bool(subgroups)
            gene_total = f"{region}_total" if has_subgroups else str(region)
            key_total = (str(region), gene_total)
            keys.add(key_total)

            sample_values = values.setdefault(sample_name, {})
            sample_values[key_total] = read_stats_to_fields(locus.get("reads"))

            if isinstance(subgroups, dict):
                for gene, subgroup in subgroups.items():
                    key = (str(region), str(gene))
                    keys.add(key)
                    subgroup_reads = (
                        subgroup.get("reads") if isinstance(subgroup, dict) else None
                    )
                    sample_values[key] = read_stats_to_fields(subgroup_reads)

    rows: List[List[str]] = []
    for region, gene in sorted(keys):
        for sample_name in sample_order:
            entry = values.get(sample_name, {}).get((region, gene), {})
            rows.append(
                [
                    sample_display(sample_name),
                    region,
                    gene,
                    format_int(entry.get("total_reads"), missing.numeric),
                    format_int(entry.get("hifi_reads"), missing.numeric),
                    format_float_trim_nonnegative(
                        entry.get("median_read_quality"), missing.numeric
                    ),
                    format_float_trim_nonnegative(
                        entry.get("median_read_passes"), missing.numeric
                    ),
                    format_float_trim_nonnegative(
                        entry.get("median_read_length"), missing.numeric
                    ),
                ]
            )

    return CsvTable(columns=COVERAGE_OUTPUT_COLUMNS, rows=rows)


def parse_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def format_int(value: Any, missing: str) -> str:
    numeric = parse_float(value)
    if numeric is None or numeric < 0:
        return missing
    return str(int(math.floor(numeric + 0.5)))


def format_float_trim(value: Any, missing: str, decimals: int = 2) -> str:
    numeric = parse_float(value)
    if numeric is None:
        return missing
    formatted = f"{numeric:.{decimals}f}"
    return formatted.rstrip("0").rstrip(".")


def format_float_trim_nonnegative(value: Any, missing: str, decimals: int = 2) -> str:
    numeric = parse_float(value)
    if numeric is None or numeric < 0:
        return missing
    return format_float_trim(numeric, missing, decimals=decimals)


def format_quality(value: Any, missing: str = "") -> str:
    numeric = parse_float(value)
    if numeric is None or numeric < 0:
        return missing
    try:
        q_value = accuracy_to_q(numeric / 100.0)
    except ValueError:
        return missing
    if q_value is None:
        return missing
    capped_q = min(q_value, 60.0)
    return f"Q{capped_q:.0f}"


def format_percentage(value: Any, missing: str = "") -> str:
    numeric = parse_float(value)
    if numeric is None:
        return missing
    return f"{numeric * 100:.2f}%"


def format_bool(value: Any, missing: str = "") -> str:
    if value is None:
        return missing
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "false"}:
            return lowered.upper()
        return value
    numeric = parse_float(value)
    if numeric is None:
        return missing
    if numeric == 0:
        return "FALSE"
    if numeric == 1:
        return "TRUE"
    return missing


def format_string(value: Any, missing: str = "") -> str:
    if value is None:
        return missing
    return str(value)

def get_nested(container: Any, key: str, default: Any = None) -> Any:
    if not isinstance(container, dict):
        return default
    return container.get(key, default)


def run_parser(args: argparse.Namespace) -> None:
    if "all" in args.output_types:
        output_types = ["trgt", "paraphase", "coverage"]
    else:
        output_types = args.output_types

    input_files: list[Path] = list(args.inputs) if args.inputs else []

    if args.input_path:
        if not args.input_path.is_dir():
            logger.error(f"Input path is not a directory: {args.input_path}")
            sys.exit(1)
        json_files = list(args.input_path.glob("*.json"))
        if not json_files:
            logger.warning(f"No JSON files found in {args.input_path}")
        input_files.extend(json_files)

    input_files = sorted(input_files, key=lambda path: str(path))

    if not input_files:
        logger.error("No input files provided. Use --inputs or --input-path.")
        sys.exit(1)

    try:
        missing_paths = [path for path in input_files if not path.exists()]
        if missing_paths:
            missing_list = "\n  ".join(str(p) for p in missing_paths)
            logger.error(f"The following input file(s) do not exist:\n  {missing_list}")
            sys.exit(1)

        logger.info(f"Processing {len(input_files)} input file(s)")

        reports = []
        for path in input_files:
            logger.debug(f"Loading {path}")
            try:
                report = load_report(path)
                reports.append(report)
            except Exception as e:
                logger.error(f"Failed to load {path}: {e}")
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug(traceback.format_exc())
                sys.exit(1)

        if not reports:
            logger.error("No inputs provided")
            sys.exit(1)

        sample_reports = [r for r in reports if r.report_type == "sample"]
        aggregate_reports = [r for r in reports if r.report_type == "aggregate"]

        logger.info(
            f"Processing {len(sample_reports)} sample report(s) and "
            f"{len(aggregate_reports)} aggregate report(s)"
        )

        missing = (
            MissingPolicy(numeric="NA", string="NA")
            if args.missing_as_na
            else MissingPolicy()
        )

        all_samples: set[str] = set()
        all_samples.update(sample_name_from_report(r) for r in sample_reports)
        for aggregate_report in aggregate_reports:
            all_samples.update(aggregate_stats_by_sample(aggregate_report))

        sample_order = sorted(all_samples)
        sample_map = (
            {sample: f"Sample_{i + 1}" for i, sample in enumerate(sorted(all_samples))}
            if args.anon
            else {}
        )

        def sample_display(name: str) -> str:
            return sample_map.get(name, name)

        tables: Dict[str, CsvTable] = {}

        try:
            if sample_reports:
                logger.debug(f"Processing {len(sample_reports)} sample report(s)")
                for output_type in output_types:
                    if output_type == "trgt":
                        tables["trgt"] = table_from_trgt_reports(
                            sample_reports, sample_order, sample_display, missing
                        )
                    elif output_type == "paraphase":
                        tables["paraphase"] = table_from_paraphase_reports(
                            sample_reports, sample_order, sample_display, missing
                        )
                    elif output_type == "coverage":
                        tables["coverage"] = table_from_coverage_reports(
                            sample_reports, sample_order, sample_display, missing
                        )

            if aggregate_reports:
                logger.debug(f"Processing {len(aggregate_reports)} aggregate report(s)")
                for idx, aggregate_report in enumerate(aggregate_reports):
                    if len(aggregate_reports) == 1:
                        key = "aggregate"
                    else:
                        key = f"aggregate_{idx + 1}"
                    tables[key] = table_from_aggregate_and_sample_reports(
                        aggregate_report,
                        sample_reports,
                        sample_order,
                        sample_display,
                        missing,
                    )
        except Exception as e:
            logger.error(f"Failed to process reports: {e}")
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug(traceback.format_exc())
            sys.exit(1)

        try:
            output_stem = args.output.stem
            output_parent = args.output.parent
            output_parent.mkdir(parents=True, exist_ok=True)

            if args.anon:
                logger.info(f"Anonymized {len(sample_map)} sample(s)")

            for data_type, table in tables.items():
                if data_type.startswith("aggregate"):
                    output_path = output_parent / f"{output_stem}.{data_type}.csv"
                elif len(sample_reports) > 0 and len(tables) > len(aggregate_reports):
                    output_path = output_parent / f"{output_stem}.{data_type}.csv"
                else:
                    output_path = args.output

                logger.info(f"Writing {data_type} CSV output to {output_path}")

                with output_path.open("w", newline="", encoding="utf-8") as sink:
                    writer = csv.writer(sink)
                    writer.writerow(table.columns)
                    writer.writerows(table.rows)
                    logger.info(
                        f"Successfully wrote {len(table.rows)} row(s) for {data_type}"
                    )
        except IOError as e:
            logger.error(f"Failed to write output: {e}")
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug(traceback.format_exc())
            sys.exit(1)

    except KeyboardInterrupt:
        logger.error("Interrupted by user")
        sys.exit(2)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(traceback.format_exc())
        sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert ptcp-qc aggregate or sample JSON into CSV rows.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Convert a single aggregate JSON to CSV
  # Output: output.aggregate.csv
  %(prog)s -i aggregate.json -o output.csv

  # Convert multiple sample JSONs to CSV
  # Outputs: output.trgt.csv, output.paraphase.csv, output.coverage.csv
  %(prog)s -i sample1.json sample2.json sample3.json -o output.csv

  # Extract only TRGT results from sample JSONs using --inputs
  # Output: results.trgt.csv
  %(prog)s --inputs sample*.json -o results.csv -t trgt

  # Extract TRGT and paraphase results with verbose output
  # Outputs: results.trgt.csv, results.paraphase.csv
  %(prog)s -i sample*.json -o results.csv -t trgt paraphase -v

  # Process a mix of aggregate and sample JSONs
  # Outputs: combined.trgt.csv, combined.paraphase.csv, combined.coverage.csv, combined.aggregate.csv
  %(prog)s -i aggregate.json sample1.json sample2.json -o combined.csv

  # Use --input-path for directories with many JSON files (avoids shell expansion limits)
  # Outputs: results.trgt.csv, results.paraphase.csv, results.coverage.csv
  %(prog)s --input-path /path/to/json/directory -o results.csv
""",
    )
    parser.add_argument(
        "-i",
        "--inputs",
        nargs="+",
        type=Path,
        default=[],
        help="Path(s) to aggregate or sample JSON files. Accepts multiple files or shell globs.",
    )
    parser.add_argument(
        "-p",
        "--input-path",
        type=Path,
        help="Directory path containing JSON files. Uses Python glob to find *.json files (avoids shell expansion limits).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output CSV path (required).",
    )
    parser.add_argument(
        "-t",
        "--output-types",
        nargs="+",
        choices=["trgt", "paraphase", "coverage", "all"],
        default=["all"],
        metavar="TYPE",
        help="Data types to extract from sample JSONs: trgt, paraphase, coverage, all (default: all). Only applies to sample reports.",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity (-v for INFO, -vv for DEBUG).",
    )
    parser.add_argument(
        "--missing-as-na",
        action="store_true",
        default=False,
        help="Render missing values as 'NA' (default: empty).",
    )
    parser.add_argument(
        "--anon",
        action="store_true",
        default=False,
        help="Anonymize sample names in outputs (replaces them with Sample_1, Sample_2, ...).",
    )

    args = parser.parse_args()
    setup_logging(args.verbose)
    run_parser(args)


if __name__ == "__main__":
    main()
