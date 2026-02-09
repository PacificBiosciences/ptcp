import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import utils  # noqa: E402


class TestNormalizeBamName(unittest.TestCase):
    def test_strips_hifi_reads_suffix(self):
        result = utils.normalize_bam_name(Path("sample1.hifi_reads.bam"))
        self.assertEqual(result, "sample1")

    def test_strips_fail_reads_suffix(self):
        result = utils.normalize_bam_name(Path("sample1.fail_reads.bam"))
        self.assertEqual(result, "sample1")

    def test_strips_reads_suffix_without_bam(self):
        result = utils.normalize_bam_name(Path("sample1.reads"))
        self.assertEqual(result, "sample1")

    def test_strips_reads_bam_suffix(self):
        result = utils.normalize_bam_name(Path("sample1.reads.bam"))
        self.assertEqual(result, "sample1")

    def test_plain_bam_extension(self):
        result = utils.normalize_bam_name(Path("sample1.bam"))
        self.assertEqual(result, "sample1")

    def test_preserves_sample_ending_in_bam_chars(self):
        result = utils.normalize_bam_name(Path("alabama.bam"))
        self.assertEqual(result, "alabama")


class TestFileHelpers(unittest.TestCase):
    def test_find_bam_files_respects_max_depth(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            bam_root = root / "root.bam"
            bam_root.touch()

            depth1 = root / "depth1"
            depth1.mkdir()
            bam_depth1 = depth1 / "depth1.bam"
            bam_depth1.touch()

            depth2 = depth1 / "depth2"
            depth2.mkdir()
            bam_depth2 = depth2 / "depth2.bam"
            bam_depth2.touch()

            files_depth0 = utils.find_bam_files(root, max_depth=0)
            self.assertEqual(set(files_depth0), {bam_root})

            files_depth1 = utils.find_bam_files(root, max_depth=1)
            self.assertEqual(set(files_depth1), {bam_root, bam_depth1})

            files_depth2 = utils.find_bam_files(root, max_depth=2)
            self.assertEqual(set(files_depth2), {bam_root, bam_depth1, bam_depth2})

    def test_should_skip_file(self):
        self.assertTrue(utils.should_skip_file(Path("unassigned_sample.bam")))
        self.assertTrue(utils.should_skip_file(Path("Sequencing_Control.bam")))
        self.assertFalse(utils.should_skip_file(Path("sample.bam")))

    def test_group_bam_files_skips_patterns(self):
        bam_files = [
            Path("sample1.hifi_reads.bam"),
            Path("sample1.fail_reads.bam"),
            Path("unassigned.hifi_reads.bam"),
        ]
        groups = utils.group_bam_files(bam_files)
        self.assertEqual(set(groups.keys()), {"sample1"})
        self.assertEqual(set(groups["sample1"]), {bam_files[0], bam_files[1]})

    def test_categorize_bam_files(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            hifi = root / "sample1.hifi_reads.bam"
            fail = root / "sample1.fail_reads.bam"
            reads = root / "sample1.reads.bam"
            for path in (hifi, fail, reads):
                path.touch()

            hifi_reads, fail_reads = utils.categorize_bam_files([hifi, fail, reads])
            self.assertEqual(
                set(hifi_reads),
                {str(hifi.resolve()), str(reads.resolve())},
            )
            self.assertEqual(set(fail_reads), {str(fail.resolve())})


if __name__ == "__main__":
    unittest.main()
