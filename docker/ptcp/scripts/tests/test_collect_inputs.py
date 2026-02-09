import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))


import collect_inputs  # noqa: E402


class TestParseBamArg(unittest.TestCase):
    def test_splits_semicolon_separated(self):
        result = collect_inputs.parse_bam_arg("a.bam;b.bam;c.bam")
        self.assertEqual(result, ["a.bam", "b.bam", "c.bam"])

    def test_filters_empty_strings(self):
        result = collect_inputs.parse_bam_arg("a.bam;;b.bam;")
        self.assertEqual(result, ["a.bam", "b.bam"])

    def test_empty_input(self):
        result = collect_inputs.parse_bam_arg("")
        self.assertEqual(result, [])

    def test_single_bam(self):
        result = collect_inputs.parse_bam_arg("single.bam")
        self.assertEqual(result, ["single.bam"])


class FakeBamReader:
    def __init__(self, bam_path):
        self.path = bam_path
        self.readGroupTable = [
            SimpleNamespace(StringID="rg/1/ccs", SampleName="sample_from_rg")
        ]
        self.peer = SimpleNamespace(header={"PG": [{"ID": "jasmine"}]})

    def hasBaseFeature(self, feature):
        return True


class TestCollectInputs(unittest.TestCase):
    def setUp(self):
        self.orig_cwd = os.getcwd()
        self.temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self.temp_dir.name)

        # Fake dependencies to avoid external commands
        self.original_validate = collect_inputs.validate_bam_file
        self.original_bamreader = collect_inputs.BamReader
        self.original_check_call = collect_inputs.subprocess.check_call

        collect_inputs.validate_bam_file = lambda path: None
        collect_inputs.BamReader = FakeBamReader
        collect_inputs.subprocess.check_call = self._fake_check_call
        self.check_calls = []

    def tearDown(self):
        collect_inputs.validate_bam_file = self.original_validate
        collect_inputs.BamReader = self.original_bamreader
        collect_inputs.subprocess.check_call = self.original_check_call
        os.chdir(self.orig_cwd)
        self.temp_dir.cleanup()

    def _fake_check_call(self, args):
        self.check_calls.append(args)
        # Simulate output creation
        out_idx = args.index("-o") + 1
        Path(args[out_idx]).touch()

    def test_single_bam_copies(self):
        bam_path = Path("sample1.hifi_reads.bam")
        bam_path.touch()

        samples = collect_inputs.collect_bam_files(
            hifi_bams=[str(bam_path)],
            fail_bams=[],
            nproc=2,
            biosamples_out="biosamples.txt",
            chunk_bams_out="chunk_bams.txt",
        )

        self.assertEqual(samples, ["sample1"])
        self.assertTrue(Path("sample0001_sample1.unmapped.bam").exists())
        biosamples = Path("biosamples.txt").read_text().strip().splitlines()
        self.assertEqual(biosamples, ["sample1"])
        chunk_bams = Path("chunk_bams.txt").read_text().strip().splitlines()
        self.assertEqual(chunk_bams, ["sample0001_sample1.unmapped.bam"])
        self.assertEqual(self.check_calls, [])

    def test_multiple_bams_merge(self):
        bam_a = Path("sample2.hifi_reads.bam")
        bam_b = Path("sample2.fail_reads.bam")
        for path in (bam_a, bam_b):
            path.touch()

        samples = collect_inputs.collect_bam_files(
            hifi_bams=[str(bam_a)],
            fail_bams=[str(bam_b)],
            nproc=4,
            biosamples_out="biosamples.txt",
            chunk_bams_out="chunk_bams.txt",
        )

        self.assertEqual(samples, ["sample2"])
        self.assertTrue(Path("sample0001_sample2.unmapped.bam").exists())
        self.assertEqual(len(self.check_calls), 1)
        self.assertEqual(
            self.check_calls[0][:5],
            ["pbmerge", "-j", "4", "-o", "sample0001_sample2.unmapped.bam"],
        )

    def test_chunk_bams_preserve_sample_order(self):
        bam_first = Path("beta.hifi_reads.bam")
        bam_second = Path("alpha.hifi_reads.bam")
        for path in (bam_first, bam_second):
            path.touch()

        samples = collect_inputs.collect_bam_files(
            hifi_bams=[str(bam_first), str(bam_second)],
            fail_bams=[],
            nproc=1,
            biosamples_out="biosamples.txt",
            chunk_bams_out="chunk_bams.txt",
        )

        self.assertEqual(samples, ["beta", "alpha"])
        chunk_bams = Path("chunk_bams.txt").read_text().strip().splitlines()
        self.assertEqual(
            chunk_bams,
            [
                "sample0001_beta.unmapped.bam",
                "sample0002_alpha.unmapped.bam",
            ],
        )

    def test_empty_bam_strings_filtered(self):
        """Empty strings in BAM lists should be filtered out."""
        bam_path = Path("sample1.hifi_reads.bam")
        bam_path.touch()

        samples = collect_inputs.collect_bam_files(
            hifi_bams=[str(bam_path), "", ""],
            fail_bams=["", ""],
            nproc=1,
            biosamples_out="biosamples.txt",
            chunk_bams_out="chunk_bams.txt",
        )

        self.assertEqual(samples, ["sample1"])
        self.assertEqual(len(self.check_calls), 0)


if __name__ == "__main__":
    unittest.main()
