#!/usr/bin/env python3
"""
Parse bed file containing pbaa guide intervals and extract sequences
from reference fasta.
"""

import argparse
import pysam
from collections import defaultdict, namedtuple


def parse_bed(bed_file):
    """Group intervals by `group name`.

    Args:
        bed_file (str): Path to bed file containing pbaa guide intervals.

    Returns:
        dict: group name -> list of Interval objects
    """
    Interval = namedtuple("Interval", ["chr", "start", "end", "name"])
    group_dict = defaultdict(list)
    with open(bed_file, "r") as bed:
        for line in bed:
            interval = Interval._make(line.strip().split("\t"))
            # name is in the format of "sequence name|group name"
            group_dict[interval.name.split("|")[1]].append(interval)
    return group_dict


def extract_sequences(fasta_file, intervals):
    """Extract sequence for interval from fasta_file.

    Args:
        fasta_file (str): Path to reference fasta.
        intervals (list): List of Interval objects.

    Returns:
        dict: name -> sequence
    """
    sequences = {}
    infasta = pysam.FastaFile(fasta_file)
    for chr, start, end, name in intervals:
        sequences[name] = infasta.fetch(chr, int(start), int(end))
    return sequences


def write_fasta(sequences, group):
    """Write list of sequences to fasta file.

    Args:
        sequences (dict): name -> sequence
        group (str): group name
    """
    with open(f"{group}.fasta", "w") as out_fasta:
        for name, seq in sequences.items():
            out_fasta.write(f">{name}\n{seq}\n")


def main(bed_file, fasta_file):
    group_dict = parse_bed(bed_file)
    for group, intervals in group_dict.items():
        sequences = extract_sequences(fasta_file, intervals)
        write_fasta(sequences, group)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bed_file", help="Bed file containing pbaa guide intervals.")
    parser.add_argument("fasta_file", help="Reference fasta.")
    args = parser.parse_args()
    main(args.bed_file, args.fasta_file)
