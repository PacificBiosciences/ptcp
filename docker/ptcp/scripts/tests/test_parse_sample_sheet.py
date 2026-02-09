import csv
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import parse_sample_sheet  # noqa: E402
import utils  # noqa: E402


class TestParseSampleSheet(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.csv_path = os.path.join(self.temp_dir.name, "sample_sheet.csv")
        self.sample_rows = [
            {"bam_id": "sample1", "bam_name": "sample_1", "sex": "M"},
            {"bam_id": "sample2", "bam_name": "sample_2", "sex": "F"},
        ]
        with open(self.csv_path, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=self.sample_rows[0].keys())
            writer.writeheader()
            writer.writerows(self.sample_rows)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_parse_tsv_returns_expected_dict(self):
        data = parse_sample_sheet.parse_tsv(self.csv_path)
        self.assertEqual(set(data.keys()), {"sample1", "sample2"})
        self.assertEqual(data["sample1"]["bam_name"], "sample_1")
        self.assertEqual(data["sample2"]["sex"], "F")

    def test_parse_tsv_raises_on_duplicate_bam_id(self):
        duplicate_path = os.path.join(self.temp_dir.name, "duplicate.csv")
        with open(duplicate_path, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=["bam_id", "bam_name", "sex"])
            writer.writeheader()
            writer.writerow({"bam_id": "dup", "bam_name": "dup_name", "sex": "M"})
            writer.writerow({"bam_id": "dup", "bam_name": "dup_other", "sex": "F"})
        with self.assertRaisesRegex(ValueError, "Duplicate bam_id found: dup"):
            parse_sample_sheet.parse_tsv(duplicate_path)

    def test_search_data_returns_row(self):
        data = parse_sample_sheet.parse_tsv(self.csv_path)
        result = parse_sample_sheet.search_data(data, "sample2")
        self.assertEqual(result, self.sample_rows[1])

    def test_search_data_raises_when_missing(self):
        data = parse_sample_sheet.parse_tsv(self.csv_path)
        with self.assertRaisesRegex(
            ValueError,
            "Error: bam_id 'missing' not found in sample sheet. Available bam_ids are: sample1, sample2",
        ):
            parse_sample_sheet.search_data(data, "missing")

    def test_reads_bam_suffix_normalizes_to_bam_id(self):
        data = parse_sample_sheet.parse_tsv(self.csv_path)
        normalized = utils.normalize_bam_name(Path("sample1.reads.bam"))
        result = parse_sample_sheet.search_data(data, normalized)
        self.assertEqual(result, self.sample_rows[0])

    def test_format_json_outputs_full_row(self):
        data = parse_sample_sheet.parse_tsv(self.csv_path)
        entry = data["sample1"]
        json_output = parse_sample_sheet.format_json(entry)
        self.assertEqual(json.loads(json_output), entry)

    def test_format_json_single_value(self):
        data = parse_sample_sheet.parse_tsv(self.csv_path)
        entry = data["sample1"]
        sex_value = parse_sample_sheet.format_json(entry, "sex")
        self.assertEqual(sex_value, "M")

    def test_two_column_sheet_supported(self):
        two_col_path = os.path.join(self.temp_dir.name, "two_col.csv")
        with open(two_col_path, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=["bam_id", "sex"])
            writer.writeheader()
            writer.writerow({"bam_id": "sample3", "sex": "F"})
        data = parse_sample_sheet.parse_tsv(two_col_path)
        self.assertEqual(set(data.keys()), {"sample3"})
        self.assertEqual(data["sample3"]["bam_name"], "sample3")
        self.assertEqual(data["sample3"]["sex"], "F")


if __name__ == "__main__":
    unittest.main()
