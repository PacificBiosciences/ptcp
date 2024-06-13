#!/usr/bin/env python3
"""Parse a TSV sample sheet with movie, barcode, sample_id, and sex columns.
barcode and sample_id must be unique keys within this file.
Search for a sample_id or barcode and return the row as JSON.
Optionally return a specific value from the row.
"""

import csv
import json
import sys
import argparse


def parse_tsv(tsv_file):
    "Parse the TSV file and return a dictionary with sample_id and barcode as keys and row dicts as values."
    data_dict = {}
    sample_ids = set()
    barcodes = set()
    with open(tsv_file, "r") as tsv:
        reader = csv.DictReader(tsv, delimiter="\t")
        for row in reader:
            if row["sample_id"] in sample_ids:
                sys.exit(f"Duplicate sample_id found: {row['sample_id']}")
            if row["barcode"] in barcodes:
                sys.exit(f"Duplicate barcode found: {row['barcode']}")
            sample_ids.add(row["sample_id"])
            barcodes.add(row["barcode"])
            data_dict[row["sample_id"]] = row
            data_dict[row["barcode"]] = row
    return data_dict


def search_data(data_dict, query):
    "Lookup a key (sample_id or barcode) in the data dictionary and return the row dict."
    if query in data_dict:
        return data_dict[query]
    else:
        sys.exit("Input not found.")


def format_json(data, output=None):
    "Format the data as JSON and return the output value if specified."
    if output:
        json_output = json.dumps(data[output], indent=4)
        return json_output.replace('\n', '').replace('"', '')
    return json.dumps(data, indent=4)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tsv_file", help="The TSV file to parse.")
    parser.add_argument("query", help="The sample_id or barcode to search for.")
    parser.add_argument(
        "--output", "-o", help="The value to return. Default is all data."
    )
    args = parser.parse_args()

    if args.output:
        args.output = args.output.strip()

    data_dict = parse_tsv(args.tsv_file)
    data = search_data(data_dict, args.query)
    print(format_json(data, args.output))

if __name__ == "__main__":
    main()
