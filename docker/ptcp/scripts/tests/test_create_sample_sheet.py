import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import create_sample_sheet  # noqa: E402


class TestCreateSampleSheet(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        tmp_path = Path(self.temp_dir.name)
        self.sample1 = tmp_path / "sample1.hifi_reads.bam"
        self.sample2_fail = tmp_path / "sample2.fail_reads.bam"
        self.sample2_hifi = tmp_path / "sample2.hifi_reads.bam"
        for path in [self.sample1, self.sample2_fail, self.sample2_hifi]:
            path.touch()
        self.bam_files = [self.sample1, self.sample2_fail, self.sample2_hifi]

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_two_column_output_default(self):
        lines, groups = create_sample_sheet.create_sample_sheet(
            self.bam_files, include_bam_name=False
        )

        self.assertEqual(lines[0], "bam_id,sex")
        self.assertEqual(set(lines[1:]), {"sample1,F", "sample2,F"})
        self.assertEqual(set(groups.keys()), {"sample1", "sample2"})

    def test_three_column_output_with_flag(self):
        lines, groups = create_sample_sheet.create_sample_sheet(
            self.bam_files, include_bam_name=True
        )
        self.assertEqual(lines[0], "bam_name,bam_id,sex")
        # sample2.hifi_reads.bam should be preferred over fail for bam_name
        self.assertEqual(
            set(lines[1:]),
            {"sample1.hifi_reads.bam,sample1,F", "sample2.hifi_reads.bam,sample2,F"},
        )
        self.assertEqual(set(groups.keys()), {"sample1", "sample2"})


if __name__ == "__main__":
    unittest.main()
