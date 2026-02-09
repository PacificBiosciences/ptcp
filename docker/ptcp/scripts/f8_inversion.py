#!/usr/bin/env python3

__version__ = "1.0.0"

import argparse
import dataclasses
import datetime
import json
import logging
import sys
from dataclasses import asdict, dataclass
from collections import Counter
from pathlib import Path
from typing import Literal, Optional, TextIO, Dict, List

import pysam


@dataclass
class InversionRead:
    inversion_label: Optional[str]
    haplotype: Optional[str]
    ref_start: Optional[int]
    ref_end: Optional[int]
    qname: Optional[str] = None

    def __hash__(self):
        return hash(
            f"{self.inversion_label}_{self.haplotype}_{self.ref_start}_{self.ref_end}"
        )

    def __eq__(self, rhs):
        return (
            self.inversion_label == rhs.inversion_label
            and self.haplotype == rhs.haplotype
            and self.ref_start == rhs.ref_start
            and self.ref_end == rhs.ref_end
        )


@dataclass
class InversionBreakpoint:
    inversion_label: Optional[str]
    haplotype: Optional[str]
    region: Optional[str]
    count: Optional[int]


@dataclass
class InversionCall:
    has_inversion: Optional[bool]
    inversion_genotype: Optional[str]
    read_counts: Dict[str, int]
    spans: Optional[Dict[str, List[InversionBreakpoint]]]
    reads: Optional[Dict[str, InversionRead]] = None


@dataclass
class InversionConfig:
    name: str
    outer_start: int
    outer_end: int
    inner_start: int
    inner_end: int
    clip_left: int
    clip_right: int
    BP1: int
    BP2: int
    chromosome: str = "chrX"
    padding_pos: int = 10
    padding_clip: int = 20


inv1_config_hg38 = InversionConfig(
    name="inv1",
    outer_start=155005414,
    outer_end=155010478,  # TODO: Is this off-by-two (should it be 155010480)?
    inner_start=155006103,
    inner_end=155007155,
    clip_left=964,
    clip_right=2281,
    BP1=155007148,
    BP2=155148826,
    chromosome="chrX",
    padding_pos=10,
    padding_clip=20,
)
inv22_config_hg38 = InversionConfig(
    name="inv22",
    outer_start=154878269,
    outer_end=154892191,
    inner_start=154880815,
    inner_end=154890327,
    clip_left=622,
    clip_right=925,
    BP1=154890327,
    BP2=155454650,
    chromosome="chrX",
    padding_pos=10,
    padding_clip=20,
)

inv1_config_hg37 = InversionConfig(
    name="inv1",
    outer_start=154233689,
    outer_end=154238753,
    inner_start=154234378,
    inner_end=154235430,
    clip_left=964,
    clip_right=2281,
    BP1=154235423,
    BP2=154377101,
    chromosome="X",
    padding_pos=10,
    padding_clip=20,
)
inv22_config_hg37 = InversionConfig(
    name="inv22",
    outer_start=154106544,
    outer_end=154120466,
    inner_start=154109090,
    inner_end=154118602,
    clip_left=622,
    clip_right=925,
    BP1=154118602,
    BP2=154684311,
    chromosome="X",
    padding_pos=10,
    padding_clip=20,
)


def _get_label_for_read(
    read: pysam.AlignedSegment, config: InversionConfig
) -> InversionRead:
    """
    Determine the label for a single read based on alignment coordinates
    and clipping patterns relative to an inversion configuration.

    Labels correspond to:
      - '00': Read spanning the outer region.
      - '11': Read from the inner region, mapped to the outer reference (expect clips on both ends)
      - '01'/'10': Inversion reads spanning one of the breakpoints (expect asymmetric clipping).
      - None: Read does not fit any defined pattern.

    Args:
        read: A pysam AlignedSegment object.
        config: The InversionConfig dataclass instance for the inversion.

    Returns:
        The label ('00', '01', '10', '11', 'x').
    """

    def is_near(pos1: int, pos2: int, tolerance: int) -> bool:
        """Check if two positions are within a given tolerance."""
        return abs(pos1 - pos2) <= tolerance

    def clip_matches(
        observed_clip: Optional[int], expected_clip: int, tolerance: int
    ) -> bool:
        """Check if observed soft clipping matches expected length within tolerance."""
        return (
            observed_clip is not None
            and abs(observed_clip - expected_clip) <= tolerance
        )

    read_start = read.reference_start
    read_end = read.reference_end
    cigartuples = read.cigartuples
    haplotype = None
    if read.has_tag("HP"):
        haplotype = read.get_tag("HP")

    clip_len_5p = None
    if cigartuples and cigartuples[0][0] == pysam.CSOFT_CLIP:
        clip_len_5p = cigartuples[0][1]

    clip_len_3p = None
    if cigartuples and cigartuples[-1][0] == pysam.CSOFT_CLIP:
        clip_len_3p = cigartuples[-1][1]

    starts_near_outer = is_near(read_start, config.outer_start, config.padding_pos)
    ends_near_outer = is_near(read_end, config.outer_end, config.padding_pos)
    starts_near_inner = is_near(read_start, config.inner_start, config.padding_pos)
    ends_near_inner = is_near(read_end, config.inner_end, config.padding_pos)

    left_clip_matches = clip_matches(clip_len_5p, config.clip_left, config.padding_clip)
    right_clip_matches = clip_matches(
        clip_len_3p, config.clip_right, config.padding_clip
    )

    label = None
    if starts_near_outer and ends_near_outer:
        label = "00"
    elif starts_near_outer and ends_near_inner and right_clip_matches:
        label = "01"
    elif starts_near_inner and left_clip_matches and ends_near_outer:
        label = "10"
    elif (
        starts_near_inner
        and left_clip_matches
        and ends_near_inner
        and right_clip_matches
    ):
        label = "11"
    else:
        label = "x"

    return InversionRead(label, haplotype, read_start, read_end, read.qname)


def label_inversion_reads(
    bam_handle: pysam.AlignmentFile, config: InversionConfig
) -> Dict[str, InversionRead]:
    """
    Fetch reads from a BAM file within the inversion region and assign labels
    based on alignment coordinates and clipping patterns.

    Args:
        bam_handle: An open pysam AlignmentFile handle.
        config: The InversionConfig dataclass instance for the inversion.

    Returns:
        A dictionary mapping read names to their assigned InversionRead objects
        with labels ('00', '01', '10', '11', 'x').
    """
    read_info = {}

    for read in bam_handle.fetch(
        config.chromosome,
        config.outer_start - config.padding_pos,
        config.outer_end + config.padding_pos,
    ):
        if read.has_tag("rq") and read.get_tag("rq") < 0.99:
            continue
        if read.is_supplementary or read.is_secondary:
            continue
        label = _get_label_for_read(read, config)
        read_info[read.qname] = label

    return read_info


def call_inversion(
    read_tags: List[InversionRead],
    min_read: int = 2,
    sample_sex: Optional[Literal["M", "F"]] = None,
) -> InversionCall:
    """Call inversion from labelled reads"""
    has_inversion = None
    genotype = None

    def count_reads_with_inversion_label(
        reads: List[InversionRead], inversion_label: str
    ):
        return len([r for r in reads if r.inversion_label == inversion_label])

    count00 = count_reads_with_inversion_label(read_tags, "00")
    count01 = count_reads_with_inversion_label(read_tags, "01")
    count10 = count_reads_with_inversion_label(read_tags, "10")
    count11 = count_reads_with_inversion_label(read_tags, "11")
    countX = count_reads_with_inversion_label(read_tags, "x")

    if count01 + count10 + count00 + count11 >= min_read:
        if count01 + count10 >= min_read:
            has_inversion = True
            if count00 + count11 >= min_read:
                if sample_sex is None:
                    genotype = "0/1"
                elif sample_sex == "F":
                    genotype = "0/1"
            elif count00 + count11 == 0:
                if sample_sex is None:
                    # default to female
                    genotype = "1/1"
                elif sample_sex == "F":
                    genotype = "1/1"
                elif sample_sex == "M":
                    genotype = "1"
        elif count01 + count10 == 0:
            has_inversion = False
            if sample_sex is None:
                # default to female
                genotype = "0/0"
            elif sample_sex == "F":
                genotype = "0/0"
            elif sample_sex == "M":
                genotype = "0"

    return InversionCall(
        has_inversion=has_inversion,
        inversion_genotype=genotype,
        spans=None,
        read_counts={
            "00": count00,
            "01": count01,
            "10": count10,
            "11": count11,
            "x": countX,
        },
    )


def get_sample_id_from_header(bam_path: Path):
    """Get sample ID from RG SM from the bam header"""
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam_handle:
            header = bam_handle.header.to_dict()
        sample_ids = []
        rg_lines = header.get("RG")
        if rg_lines is not None:
            sample_ids = [a.get("SM") for a in rg_lines if "SM" in a]
        sample_ids = [a for a in sample_ids if a is not None]
        if len(set(sample_ids)) == 1:
            return "_".join(sample_ids[0].split())
        else:
            logging.warning(f"Multiple sample IDs found in {bam_path.name} header")
            return None
    except Exception as e:
        logging.error(f"Error reading BAM header from {bam_path.name}: {str(e)}")
        return None


def write_vcf_header(fout: TextIO, prog_cmd: str, sample_id: str, genome_version: int):
    """Write VCF header"""

    chromosome_lengths = {38: ("chrX", 156040895), 37: ("X", 155270560)}
    chrom, chrom_length = chromosome_lengths[genome_version]

    fout.write("##fileformat=VCFv4.2\n")
    fout.write('##ALT=<ID=INV,Description="Inversion">\n')
    fout.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    fout.write('##FILTER=<ID=LowDP,Description="Low depth">\n')
    fout.write(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">\n'
    )
    fout.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV">\n')
    fout.write(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n'
    )
    fout.write(
        '##INFO=<ID=EVIDENCE,Number=5,Type=Integer,Description="Number of reads supporting WT breakpoint 1, WT breakpoint 2, inversion breakpoint 1, inversion breakpoint 2, and unclassified respectively">\n'
    )
    fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fout.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
    fout.write(
        '##FORMAT=<ID=SUPP,Number=1,Type=Integer,Description="Read depth supporting the inversion">\n'
    )

    fout.write(f"##contig=<ID={chrom},length={chrom_length}>\n")
    fout.write(f"##command={prog_cmd}\n")
    header = [
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        sample_id,
    ]
    fout.write("\t".join(header) + "\n")


def inversion_vcf_line(
    config: InversionConfig,
    inv_call: InversionCall,
    inv_name: str,
) -> str:
    """Prepare the inversion record for VCF"""
    var_type = "INV"
    if inv_call.has_inversion is None:
        variant_filter = "LowDP"
        gt = "."
    else:
        variant_filter = "PASS"
        gt = inv_call.inversion_genotype or "."

    dp = sum(inv_call.read_counts.values())
    supp = inv_call.read_counts["01"] + inv_call.read_counts["10"]
    evidence = ",".join(
        str(inv_call.read_counts[key]) for key in ["00", "11", "01", "10", "x"]
    )
    info_field = f"SVTYPE={var_type};END={config.BP2};EVIDENCE={evidence}"
    vcf_line = [
        config.chromosome,
        config.BP1,
        inv_name,
        "N",
        f"<{var_type}>",
        ".",  # qual
        variant_filter,
        info_field,
        "GT:DP:SUPP",
    ] + [":".join([str(a) for a in [gt, dp, supp]])]
    return "\t".join([str(a) for a in vcf_line]) + "\n"


def configure_logging(verbose: bool):
    level = logging.DEBUG if verbose else logging.INFO
    fmt = "%(asctime)s %(levelname)-5s %(message)s"
    logging.basicConfig(level=level, format=fmt)


def get_spans(
    chrom: str,
    reads: List[InversionRead],
    min_fraction: float = 0.1,
    min_absolute: int = 30,
) -> List[InversionBreakpoint]:
    if not reads:
        return []
    frequent_reads = []
    threshold = int(max(len(reads) * min_fraction, min_absolute))
    cnt = Counter(reads)
    for pos, pos_count in cnt.items():
        if pos_count >= threshold:
            frequent_reads.append((pos, pos_count))

    if not frequent_reads:
        return []

    def overlaps(ia: InversionRead, ib: InversionRead):
        return (
            ia.inversion_label == ib.inversion_label
            and ia.haplotype == ib.haplotype
            and ia.ref_start <= ib.ref_end
            and ia.ref_end >= ib.ref_start
        )

    def merge_intervals(ia: InversionRead, ca: int, ib: InversionRead, cb: int):
        assert ia.inversion_label == ib.inversion_label
        assert ia.haplotype == ib.haplotype
        merged = InversionRead(
            ia.inversion_label,
            ia.haplotype,
            min(ia.ref_start, ib.ref_start),
            max(ia.ref_end, ib.ref_end),
            ia.qname,  # Use the first read's qname for the merged interval
        )
        return (merged, ca + cb)

    frequent_reads = sorted(
        frequent_reads,
        key=lambda x: (
            x[0].inversion_label,
            x[0].haplotype,
            x[0].ref_start,
            x[0].ref_end,
        ),
    )
    spans = []
    cur_read, cur_count = frequent_reads[0]

    for r, c in frequent_reads[1:]:
        if overlaps(cur_read, r):
            cur_read, cur_count = merge_intervals(cur_read, cur_count, r, c)
        else:
            spans.append(
                InversionBreakpoint(
                    cur_read.inversion_label,
                    cur_read.haplotype,
                    f"{chrom}:{cur_read.ref_start}-{cur_read.ref_end}",
                    cur_count,
                )
            )
            cur_read = r
            cur_count = c
    spans.append(
        InversionBreakpoint(
            cur_read.inversion_label,
            cur_read.haplotype,
            f"{chrom}:{cur_read.ref_start}-{cur_read.ref_end}",
            cur_count,
        )
    )
    return spans


def get_args():
    parser = argparse.ArgumentParser(
        description="Call F8 Intron 1 and Intron 22 inversions from the PureTarget panel data",
    )
    inputp = parser.add_argument_group("Input Options")
    outputp = parser.add_argument_group("Output Options")
    inputp.add_argument(
        "-b",
        "--bam",
        help="Input bam, phased by Paraphase",
        type=Path,
        required=True,
    )
    inputp.add_argument(
        "--genome",
        help="Genome version",
        choices=[38, 37],
        default=38,
        type=int,
        required=False,
    )
    outputp.add_argument(
        "-o",
        "--out",
        help="Output directory",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "-n",
        "--nread",
        help="Minimum number of reads required to call an allele, default is 2",
        required=False,
        type=int,
        default=2,
    )
    parser.add_argument(
        "-s",
        "--sex",
        help="Sample sex, optional.",
        choices=["F", "M", None],
        required=False,
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="Prefix to output files, optional. Dafault is to infer from bam header",
        required=False,
    )
    parser.add_argument(
        "--json",
        help="Optional. If specified, write a separate JSON output file",
        required=False,
        action="store_true",
    )
    parser.add_argument("--verbose", action="store_true", help="More verbose logging")
    return parser.parse_args()


def run(args):
    bam = args.bam
    nread = args.nread
    sample_sex = args.sex

    genome_version = args.genome

    if genome_version == 37:
        inv1_cfg, inv22_cfg = inv1_config_hg37, inv22_config_hg37
        logging.debug(f"Using hg37 coordinates (genome specified: {args.genome}).")
    elif genome_version == 38:
        inv1_cfg, inv22_cfg = inv1_config_hg38, inv22_config_hg38
        logging.debug(f"Using hg38 coordinates (genome specified: {args.genome}).")
    else:
        # Should be unreachable given argparse
        raise ValueError(f"Internal error: Unhandled genome version: {args.genome!r}")

    if sample_sex is not None:
        if sample_sex not in set(["M", "F"]):
            sample_sex = None

    outdir = args.out
    outdir.mkdir(parents=True, exist_ok=True)

    if args.prefix is not None:
        sample_id = args.prefix
    else:
        sample_id = get_sample_id_from_header(bam)
    if sample_id is None or sample_id == "":
        sample_id = bam.name.split(".", 1)[0]

    logging.info(f"Started analysis for {sample_id} at {datetime.datetime.now()}...")
    vcf_out = outdir / f"{sample_id}.f8inversion.vcf"

    with pysam.AlignmentFile(bam, "r") as bam_in:
        inv22_reads = label_inversion_reads(bam_handle=bam_in, config=inv22_cfg)
        inv1_reads = label_inversion_reads(bam_handle=bam_in, config=inv1_cfg)

    inv22_call = call_inversion(
        list(inv22_reads.values()), min_read=nread, sample_sex=sample_sex
    )
    logging.debug(f"inv22 call: {inv22_call}")
    inv1_call = call_inversion(
        list(inv1_reads.values()), min_read=nread, sample_sex=sample_sex
    )
    logging.debug(f"inv1 call: {inv1_call}")

    with open(vcf_out, "w") as vcf_handle:
        prog_cmd = " ".join(sys.argv)
        write_vcf_header(vcf_handle, prog_cmd, sample_id, genome_version)
        inv22_vcf_line = inversion_vcf_line(inv22_cfg, inv22_call, "f8inv22")
        vcf_handle.write(inv22_vcf_line)
        inv1_vcf_line = inversion_vcf_line(inv1_cfg, inv1_call, "f8inv1")
        vcf_handle.write(inv1_vcf_line)

    if args.json:
        json_out = outdir / f"{sample_id}.f8inversion.json"
        inv1_call = dataclasses.replace(
            inv1_call,
            reads=inv1_reads,
            spans=get_spans(inv1_cfg.chromosome, list(inv1_reads.values())),
        )
        inv22_call = dataclasses.replace(
            inv22_call,
            reads=inv22_reads,
            spans=get_spans(inv22_cfg.chromosome, list(inv22_reads.values())),
        )
        sample_call = {"f8inv1": asdict(inv1_call), "f8inv22": asdict(inv22_call)}
        with open(json_out, "w") as fout:
            json.dump(sample_call, fout, indent=4)

    logging.info(f"Completed analysis for {sample_id}")
    return 0


def main():
    bam_path = None
    try:
        args = get_args()
        bam_path = args.bam
        configure_logging(args.verbose)
        return run(args)
    except ValueError as ve:
        logging.error("Invalid input for BAM file %s: %s", bam_path or "unknown", ve)
        sys.exit(1)
    except Exception:
        logging.exception(
            "Unhandled exception while processing BAM file %s, aborting.",
            bam_path or "unknown",
        )
        sys.exit(2)


if __name__ == "__main__":
    sys.exit(main())
