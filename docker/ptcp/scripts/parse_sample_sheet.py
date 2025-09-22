#!/usr/bin/env python3
"""Parse a CSV sample sheet with bam_name, bam_id and sex columns.
bam_name is the full bam name without stripped of .hifi_reads/.fail_reads/.reads.
bam_id must be unique and match the BAM filename with .hifi_reads/.fail_reads/.reads stripped.
Search for a bam_id and return the row as JSON.
Optionally return a specific value from the row.
"""

import csv
import json
import sys
import argparse
from typing import Dict, Optional


def parse_tsv(csv_file: str) -> Dict[str, Dict[str, str]]:
    "Parse the CSV file and return a dictionary with bam_id as key and row dict as value."
    data_dict = {}
    bam_ids = set()
    with open(csv_file, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["bam_id"] in bam_ids:
                raise ValueError(f"Duplicate bam_id found: {row['bam_id']}")
            bam_ids.add(row["bam_id"])
            data_dict[row["bam_id"]] = row
    return data_dict


def search_data(data_dict: Dict[str, Dict[str, str]], query: str) -> Dict[str, str]:
    "Lookup a bam_prefix in the data dictionary and return the row dict."
    if query in data_dict:
        return data_dict[query]
    else:
        available_ids = ", ".join(sorted(data_dict.keys()))
        raise ValueError(
            f"Error: bam_id '{query}' not found in sample sheet. Available bam_ids are: {available_ids}"
        )


def format_json(data: Dict[str, str], output: Optional[str] = None) -> str:
    "Format the data as JSON and return the output value if specified."
    if output:
        json_output = json.dumps(data[output], indent=4)
        return json_output.replace("\n", "").replace('"', "")
    return json.dumps(data, indent=4)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_file", help="The CSV file to parse.")

    query_group = parser.add_mutually_exclusive_group(required=True)
    query_group.add_argument(
        "query",
        help="The bam_prefix to search for",
        nargs="?",
    )
    query_group.add_argument(
        "--batch-input",
        "-b",
        action="store_true",
        help="Read multiple queries from stdin, one per line",
    )

    parser.add_argument(
        "--output", "-o", help="The value to return. Default is all data."
    )
    args = parser.parse_args()

    if args.output:
        args.output = args.output.strip()

    data_dict = parse_tsv(args.csv_file)
    errors = []

    if args.batch_input:
        for query in sys.stdin:
            query = query.strip()
            if query:  # Skip empty lines
                try:
                    data = search_data(data_dict, query)
                    print(format_json(data, args.output))
                except ValueError as e:
                    errors.append(str(e))
        if errors:
            sys.exit(
                f"Error(s) occurred processing {len(errors)} sample(s):\n"
                + "\n".join(errors)
            )
    else:
        try:
            data = search_data(data_dict, args.query)
            print(format_json(data, args.output))
        except ValueError as e:
            print(str(e), file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
