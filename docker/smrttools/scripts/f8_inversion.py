#!/usr/bin/env python3

import argparse
import os
import pysam
import re
import traceback
import logging
import datetime
import sys
import json


INV1_BP1 = 155007148
INV1_BP2 = 155148826
INV22_BP1 = 154890327
INV22_BP2 = 155454650


def inv_reads(
    bam_handle,
    outer_start,
    outer_end,
    inner_start,
    inner_end,
    clip_left,
    clip_right,
):
    """Assign reads to four classes of molecules"""
    read_info = {}
    nchr = "chrX"
    padding_pos = 10
    padding_clip = 20
    clip_5p = r"^\d+S"
    clip_3p = r"\d+S$"
    for read in bam_handle.fetch(
        nchr, outer_start - padding_pos, outer_end + padding_pos
    ):
        read_qc = True
        if read.has_tag("rq"):
            read_rq = read.get_tag("rq")
            if read_rq < 0.99:
                read_qc = False
        if read_qc:
            read_name = read.qname
            read_info.setdefault(read_name, None)
            read_start = read.reference_start
            read_end = read.reference_end
            read_cigar = read.cigarstring
            clip_len_5p = None
            find_clip_5p = re.findall(clip_5p, read_cigar)
            if find_clip_5p != []:
                clip_len_5p = int(read_cigar.split("S")[0])
            clip_len_3p = None
            find_clip_3p = re.findall(clip_3p, read_cigar)
            if find_clip_3p != []:
                clip_len_3p = int(find_clip_3p[0][:-1])
            if abs(read_start - outer_start) <= padding_pos:
                if abs(read_end - outer_end) <= padding_pos:
                    read_info[read_name] = "00"
                elif (
                    abs(read_end - inner_end) <= padding_pos
                    and clip_len_3p is not None
                    and abs(clip_len_3p - clip_right) <= padding_clip
                ):
                    read_info[read_name] = "01"
            elif (
                abs(read_start - inner_start) <= padding_pos
                and clip_len_5p is not None
                and abs(clip_len_5p - clip_left) <= padding_clip
            ):
                if abs(read_end - outer_end) <= padding_pos:
                    read_info[read_name] = "10"
                elif (
                    abs(read_end - inner_end) <= padding_pos
                    and clip_len_3p is not None
                    and abs(clip_len_3p - clip_right) <= padding_clip
                ):
                    read_info[read_name] = "11"
    return read_info


def call_inv(read_tags, min_read=2, sample_sex=None):
    """Call inversion from labelled reads"""
    has_inversion = None
    genotype = None
    count00 = read_tags.count("00")
    count01 = read_tags.count("01")
    count10 = read_tags.count("10")
    count11 = read_tags.count("11")
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
    return {
        "has_inversion": has_inversion,
        "inversion_genotype": genotype,
        "read_counts": {
            "00": count00,
            "01": count01,
            "10": count10,
            "11": count11,
        },
    }


def get_sample_id_from_header(bam):
    """Get sample ID from RG SM from the bam header"""
    bam_handle = pysam.AlignmentFile(bam, "rb")
    header = bam_handle.header
    header = header.to_dict()
    sample_ids = []
    rg_lines = header.get("RG")
    if rg_lines is not None:
        sample_ids = [a.get("SM") for a in rg_lines if "SM" in a]
    bam_handle.close()
    sample_ids = [a for a in sample_ids if a is not None]
    if len(set(sample_ids)) == 1:
        return "_".join(sample_ids[0].split())
    else:
        return None


def write_vcf_header(fout, prog_cmd, sample_id):
    """Write VCF header"""
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
        '##INFO=<ID=EVIDENCE,Number=4,Type=Integer,Description="Number of reads supporting WT breakpoint 1, WT breakpoint 2, inversion breakpoint 1 and inversion breakpoint 2, respectively">\n'
    )
    fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fout.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
    fout.write(
        '##FORMAT=<ID=SUPP,Number=1,Type=Integer,Description="Read depth supporting the inversion">\n'
    )
    fout.write("##contig=<ID=chrX,length=156040895>\n")
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


def inversion_vcf_line(nstart, nend, inv_call, inv_name):
    """Prepare the inversion record for VCF"""
    var_type = "INV"
    if inv_call["has_inversion"] is None:
        variant_filter = "LowDP"
        gt = "."
    else:
        variant_filter = "PASS"
        inversion_genotype = inv_call["inversion_genotype"]
        if inversion_genotype is None:
            gt = "."
        else:
            gt = inversion_genotype
    dp = sum(list(inv_call["read_counts"].values()))
    supp = inv_call["read_counts"]["01"] + inv_call["read_counts"]["10"]
    evidence = ",".join(
        [
            str(a)
            for a in [
                inv_call["read_counts"]["00"],
                inv_call["read_counts"]["11"],
                inv_call["read_counts"]["01"],
                inv_call["read_counts"]["10"],
            ]
        ]
    )
    info_field = f"SVTYPE={var_type};END={nend};EVIDENCE={evidence}"
    vcf_line = [
        "chrX",
        nstart,
        inv_name,
        "N",
        f"<{var_type}>",
        ".",  # qual
        variant_filter,
        info_field,
        "GT:DP:SUPP",
    ] + [":".join([str(a) for a in [gt, dp, supp]])]
    return "\t".join([str(a) for a in vcf_line]) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description="Call F8 Intron 1 and Intron 22 inversions from Puretarget panel data",
    )
    inputp = parser.add_argument_group("Input Options")
    outputp = parser.add_argument_group("Output Options")
    inputp.add_argument(
        "-b",
        "--bam",
        help="Input bam, realigned to a masked genome (hg38))",
        required=True,
    )
    outputp.add_argument(
        "-o",
        "--out",
        help="Output directory",
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
        help="Sample sex, optional. Accepted values are F or M",
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

    logging.basicConfig(level=logging.INFO)
    try:
        # parse iput parameters
        args = parser.parse_args()
        bam = args.bam
        nread = args.nread
        sample_sex = args.sex
        if sample_sex is not None:
            if sample_sex not in ["M", "F"]:
                sample_sex = None
        outdir = args.out
        os.makedirs(outdir, exist_ok=True)
        if args.prefix is not None:
            sample_id = args.prefix
        else:
            sample_id = get_sample_id_from_header(bam)
        if sample_id is None or sample_id == "":
            sample_id = bam.split("/")[-1].split(".")[0]
        logging.info(
            f"Started analysis for {sample_id} at {datetime.datetime.now()}..."
        )
        vcf_out = os.path.join(outdir, f"{sample_id}.f8inversion.vcf")
        if args.json:
            json_out = os.path.join(outdir, f"{sample_id}.f8inversion.json")
        bam_handle = pysam.AlignmentFile(bam)

        # label reads
        inv1_reads = inv_reads(
            bam_handle=bam_handle,
            outer_start=155004820,
            outer_end=155008983,
            inner_start=155006107,
            inner_end=155007147,
            clip_left=969,
            clip_right=2290,
        )
        inv22_reads = inv_reads(
            bam_handle=bam_handle,
            outer_start=154879365,
            outer_end=154892105,
            inner_start=154880815,
            inner_end=154890318,
            clip_left=622,
            clip_right=935,
        )
        bam_handle.close()

        # call inversions and write to vcf
        inv1_call = call_inv(
            list(inv1_reads.values()), min_read=nread, sample_sex=sample_sex
        )
        inv22_call = call_inv(
            list(inv22_reads.values()), min_read=nread, sample_sex=sample_sex
        )
        vcf_handle = open(vcf_out, "w")
        prog_cmd = " ".join(sys.argv)
        write_vcf_header(vcf_handle, prog_cmd, sample_id)
        inv22_vcf_line = inversion_vcf_line(INV22_BP1, INV22_BP2, inv22_call, "inv22")
        vcf_handle.write(inv22_vcf_line)
        inv1_vcf_line = inversion_vcf_line(INV1_BP1, INV1_BP2, inv1_call, "inv1")
        vcf_handle.write(inv1_vcf_line)
        vcf_handle.close()

        # write to json
        if args.json:
            inv1_call.setdefault("reads", inv1_reads)
            inv22_call.setdefault("reads", inv22_reads)
            sample_call = {"inv1": inv1_call, "inv22": inv22_call}
            with open(json_out, "w") as fout:
                json.dump(sample_call, fout, indent=4)

    except Exception:
        logging.error("Error running the program...See error message below")
        traceback.print_exc()
    finally:
        logging.info(f"Completed analysis at {datetime.datetime.now()}...")


if __name__ == "__main__":
    main()
