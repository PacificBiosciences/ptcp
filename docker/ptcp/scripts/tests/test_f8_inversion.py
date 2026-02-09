import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from f8_inversion import (  # noqa: E402
    _get_label_for_read,
    call_inversion,
    get_spans,
    inv1_config_hg38,
    inv22_config_hg38,
    InversionRead,
)

try:
    import pysam

    CSOFT_CLIP = pysam.CSOFT_CLIP
    CMATCH = pysam.CMATCH
except ImportError:
    CSOFT_CLIP = 4  # S operation in SAM
    CMATCH = 0  # M operation in SAM


def create_mock_read(
    start: int,
    end: int,
    cigartuples: list = None,
    haplotype: int = None,
    qname: str = "test_read",
):
    read = SimpleNamespace()
    read.reference_start = start
    read.reference_end = end
    read.cigartuples = cigartuples
    read.qname = qname
    read.has_tag = lambda tag: tag == "HP" and haplotype is not None
    read.get_tag = lambda tag: haplotype if tag == "HP" else None
    return read


class TestGetLabelForRead(unittest.TestCase):
    def test_inv1_label_00(self):
        config = inv1_config_hg38
        read = create_mock_read(
            config.outer_start,
            config.outer_end,
            [(CMATCH, config.outer_end - config.outer_start)],
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "00")
        self.assertIsNone(result.haplotype)
        # Within padding
        read = create_mock_read(
            config.outer_start + config.padding_pos,
            config.outer_end - config.padding_pos,
            [(CMATCH, 100)],
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "00")
        self.assertIsNone(result.haplotype)
        # With small clips
        read = create_mock_read(
            config.outer_start,
            config.outer_end,
            [(CSOFT_CLIP, 5), (CMATCH, 100), (CSOFT_CLIP, 5)],
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "00")
        self.assertIsNone(result.haplotype)

    def test_inv1_label_11(self):
        config = inv1_config_hg38
        match_len = config.inner_end - config.inner_start
        # Exact match
        cigar = [
            (CSOFT_CLIP, config.clip_left),
            (CMATCH, match_len),
            (CSOFT_CLIP, config.clip_right),
        ]
        read = create_mock_read(config.inner_start, config.inner_end, cigar)
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "11")
        self.assertIsNone(result.haplotype)
        # Within padding (coordinates and clips)
        cigar_padded = [
            (CSOFT_CLIP, config.clip_left - config.padding_clip),
            (CMATCH, match_len),
            (CSOFT_CLIP, config.clip_right + config.padding_clip),
        ]
        read = create_mock_read(
            config.inner_start + config.padding_pos,
            config.inner_end - config.padding_pos,
            cigar_padded,
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "11")
        self.assertIsNone(result.haplotype)

    def test_inv1_label_01(self):
        config = inv1_config_hg38
        match_len = config.inner_end - config.outer_start  # Approx match len
        # Exact match
        cigar = [(CMATCH, match_len), (CSOFT_CLIP, config.clip_right)]
        read = create_mock_read(config.outer_start, config.inner_end, cigar)
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "01")
        self.assertIsNone(result.haplotype)
        # Within padding (coordinates and clip)
        cigar_padded = [
            (CMATCH, match_len),
            (CSOFT_CLIP, config.clip_right - config.padding_clip),
        ]
        read = create_mock_read(
            config.outer_start + config.padding_pos,
            config.inner_end - config.padding_pos,
            cigar_padded,
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "01")
        self.assertIsNone(result.haplotype)

    def test_inv1_label_10(self):
        config = inv1_config_hg38
        match_len = config.outer_end - config.inner_start  # Approx match len
        # Exact match
        cigar = [(CSOFT_CLIP, config.clip_left), (CMATCH, match_len)]
        read = create_mock_read(config.inner_start, config.outer_end, cigar)
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "10")
        self.assertIsNone(result.haplotype)
        # Within padding (coordinates and clip)
        cigar_padded = [
            (CSOFT_CLIP, config.clip_left + config.padding_clip),
            (CMATCH, match_len),
        ]
        read = create_mock_read(
            config.inner_start - config.padding_pos,
            config.outer_end + config.padding_pos,
            cigar_padded,
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "10")
        self.assertIsNone(result.haplotype)

    def test_inv1_label_none(self):
        config = inv1_config_hg38
        # Start coord outside padding
        read = create_mock_read(
            config.outer_start - config.padding_pos - 5,  # Too far left
            config.outer_end,
            [(CMATCH, 100)],
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "x")
        # End coord outside padding
        read = create_mock_read(
            config.outer_start,
            config.outer_end + config.padding_pos + 5,  # Too far right
            [(CMATCH, 100)],
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "x")
        # Inner start coord outside padding
        read = create_mock_read(
            config.inner_start - config.padding_pos - 5,  # Too far left for inner start
            config.outer_end,
            [(CSOFT_CLIP, config.clip_left), (CMATCH, 100)],  # Correct clip for '10'
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "x")
        # Inner end coord outside padding
        read = create_mock_read(
            config.outer_start,
            config.inner_end + config.padding_pos + 5,  # Too far right for inner end
            [(CMATCH, 100), (CSOFT_CLIP, config.clip_right)],  # Correct clip for '01'
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "x")

    def test_inv22_label_00(self):
        config = inv22_config_hg38
        read = create_mock_read(
            config.outer_start,
            config.outer_end,
            [(CMATCH, config.outer_end - config.outer_start)],
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "00")
        self.assertIsNone(result.haplotype)
        read = create_mock_read(
            config.outer_start + config.padding_pos,
            config.outer_end - config.padding_pos,
            [(CMATCH, 100)],
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "00")
        self.assertIsNone(result.haplotype)

    def test_inv22_label_11(self):
        config = inv22_config_hg38
        match_len = config.inner_end - config.inner_start
        cigar = [
            (CSOFT_CLIP, config.clip_left),
            (CMATCH, match_len),
            (CSOFT_CLIP, config.clip_right),
        ]
        read = create_mock_read(config.inner_start, config.inner_end, cigar)
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "11")
        self.assertIsNone(result.haplotype)
        cigar_padded = [
            (CSOFT_CLIP, config.clip_left - config.padding_clip),
            (CMATCH, match_len),
            (CSOFT_CLIP, config.clip_right + config.padding_clip),
        ]
        read = create_mock_read(
            config.inner_start + config.padding_pos,
            config.inner_end - config.padding_pos,
            cigar_padded,
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "11")
        self.assertIsNone(result.haplotype)

    def test_inv22_label_01(self):
        config = inv22_config_hg38
        match_len = config.inner_end - config.outer_start
        cigar = [(CMATCH, match_len), (CSOFT_CLIP, config.clip_right)]
        read = create_mock_read(config.outer_start, config.inner_end, cigar)
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "01")
        self.assertIsNone(result.haplotype)
        cigar_padded = [
            (CMATCH, match_len),
            (CSOFT_CLIP, config.clip_right - config.padding_clip),
        ]
        read = create_mock_read(
            config.outer_start + config.padding_pos,
            config.inner_end - config.padding_pos,
            cigar_padded,
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "01")
        self.assertIsNone(result.haplotype)

    def test_inv22_label_10(self):
        config = inv22_config_hg38
        match_len = config.outer_end - config.inner_start
        cigar = [(CSOFT_CLIP, config.clip_left), (CMATCH, match_len)]
        read = create_mock_read(config.inner_start, config.outer_end, cigar)
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "10")
        self.assertIsNone(result.haplotype)
        cigar_padded = [
            (CSOFT_CLIP, config.clip_left + config.padding_clip),
            (CMATCH, match_len),
        ]
        read = create_mock_read(
            config.inner_start - config.padding_pos,
            config.outer_end + config.padding_pos,
            cigar_padded,
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "10")
        self.assertIsNone(result.haplotype)

    def test_inv22_label_none(self):
        config = inv22_config_hg38
        # Start coord outside padding
        read = create_mock_read(
            config.outer_start - config.padding_pos - 5,  # Too far left
            config.outer_end,
            [(CMATCH, 100)],
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "x")
        # End coord outside padding
        read = create_mock_read(
            config.outer_start,
            config.outer_end + config.padding_pos + 5,  # Too far right
            [(CMATCH, 100)],
        )
        result = _get_label_for_read(read, config)
        self.assertEqual(result.inversion_label, "x")


class TestCallInversion(unittest.TestCase):
    def test_homozygous_reference_female(self):
        reads = [InversionRead("00", None, 100, 200, "read_00")] * 5 + [
            InversionRead("11", None, 300, 400, "read_11")
        ] * 5
        result = call_inversion(reads, min_read=2, sample_sex="F")
        self.assertFalse(result.has_inversion)
        self.assertEqual(result.inversion_genotype, "0/0")
        self.assertEqual(
            result.read_counts, {"00": 5, "01": 0, "10": 0, "11": 5, "x": 0}
        )

    def test_homozygous_reference_male(self):
        reads = [InversionRead("00", None, 100, 200, "read_00")] * 5 + [
            InversionRead("11", None, 300, 400, "read_11")
        ] * 5
        result = call_inversion(reads, min_read=2, sample_sex="M")
        self.assertFalse(result.has_inversion)
        self.assertEqual(result.inversion_genotype, "0")
        self.assertEqual(
            result.read_counts, {"00": 5, "01": 0, "10": 0, "11": 5, "x": 0}
        )

    def test_homozygous_reference_sex_unknown(self):
        reads = [InversionRead("00", None, 100, 200, "read_00")] * 5 + [
            InversionRead("11", None, 300, 400, "read_11")
        ] * 5
        result = call_inversion(reads, min_read=2, sample_sex=None)
        self.assertFalse(result.has_inversion)
        self.assertEqual(result.inversion_genotype, "0/0")
        self.assertEqual(
            result.read_counts, {"00": 5, "01": 0, "10": 0, "11": 5, "x": 0}
        )

    def test_heterozygous_inversion(self):
        reads = (
            [InversionRead("00", None, 100, 200, "read_00")] * 5
            + [InversionRead("11", None, 300, 400, "read_11")] * 5
            + [InversionRead("01", None, 500, 600, "read_01")] * 3
            + [InversionRead("10", None, 700, 800, "read_10")] * 3
        )
        result = call_inversion(reads, min_read=2, sample_sex="F")
        self.assertTrue(result.has_inversion)
        self.assertEqual(result.inversion_genotype, "0/1")
        self.assertEqual(
            result.read_counts, {"00": 5, "01": 3, "10": 3, "11": 5, "x": 0}
        )

    def test_homozygous_inversion_female(self):
        reads = [InversionRead("01", None, 100, 200, "read_01")] * 5 + [
            InversionRead("10", None, 300, 400, "read_10")
        ] * 5
        result = call_inversion(reads, min_read=2, sample_sex="F")
        self.assertTrue(result.has_inversion)
        self.assertEqual(result.inversion_genotype, "1/1")
        self.assertEqual(
            result.read_counts, {"00": 0, "01": 5, "10": 5, "11": 0, "x": 0}
        )

    def test_homozygous_inversion_male(self):
        reads = [InversionRead("01", None, 100, 200, "read_01")] * 5 + [
            InversionRead("10", None, 300, 400, "read_10")
        ] * 5
        result = call_inversion(reads, min_read=2, sample_sex="M")
        self.assertTrue(result.has_inversion)
        self.assertEqual(result.inversion_genotype, "1")
        self.assertEqual(
            result.read_counts, {"00": 0, "01": 5, "10": 5, "11": 0, "x": 0}
        )

    def test_homozygous_inversion_sex_unknown(self):
        reads = [InversionRead("01", None, 100, 200, "read_01")] * 5 + [
            InversionRead("10", None, 300, 400, "read_10")
        ] * 5
        result = call_inversion(reads, min_read=2, sample_sex=None)
        self.assertTrue(result.has_inversion)
        self.assertEqual(result.inversion_genotype, "1/1")
        self.assertEqual(
            result.read_counts, {"00": 0, "01": 5, "10": 5, "11": 0, "x": 0}
        )

    def test_low_depth_no_call(self):
        reads = [InversionRead("00", None, 100, 200, "read_00")] * 1
        result = call_inversion(reads, min_read=2)
        self.assertIsNone(result.has_inversion)
        self.assertIsNone(result.inversion_genotype)
        self.assertEqual(
            result.read_counts, {"00": 1, "01": 0, "10": 0, "11": 0, "x": 0}
        )

    def test_low_depth_inversion_reads_no_call(self):
        reads = (
            [InversionRead("00", None, 100, 200, "read_00")] * 5
            + [InversionRead("11", None, 300, 400, "read_11")] * 5
            + [InversionRead("01", None, 500, 600, "read_01")] * 1
        )
        result = call_inversion(reads, min_read=2)
        self.assertIsNone(result.has_inversion)
        self.assertIsNone(result.inversion_genotype)
        self.assertEqual(
            result.read_counts, {"00": 5, "01": 1, "10": 0, "11": 5, "x": 0}
        )

    def test_low_depth_reference_reads_hom_inv_call_female(self):
        reads = [InversionRead("01", None, 100, 200, "read_01")] * 3 + [
            InversionRead("10", None, 300, 400, "read_10")
        ] * 3
        result = call_inversion(reads, min_read=2, sample_sex="F")
        self.assertTrue(result.has_inversion)
        self.assertEqual(result.inversion_genotype, "1/1")
        self.assertEqual(
            result.read_counts, {"00": 0, "01": 3, "10": 3, "11": 0, "x": 0}
        )

    def test_low_depth_reference_reads_hom_inv_call_male(self):
        reads = [InversionRead("01", None, 100, 200, "read_01")] * 3 + [
            InversionRead("10", None, 300, 400, "read_10")
        ] * 3
        result = call_inversion(reads, min_read=2, sample_sex="M")
        self.assertTrue(result.has_inversion)
        self.assertEqual(result.inversion_genotype, "1")
        self.assertEqual(
            result.read_counts, {"00": 0, "01": 3, "10": 3, "11": 0, "x": 0}
        )

    def test_with_none_reads(self):
        reads = (
            [InversionRead("00", None, 100, 200, "read_00")] * 5
            + [InversionRead("11", None, 300, 400, "read_11")] * 5
            + [InversionRead("x", None, 300, 400, "read_x")] * 10
        )
        result = call_inversion(reads, min_read=2, sample_sex="F")
        self.assertFalse(result.has_inversion)
        self.assertEqual(result.inversion_genotype, "0/0")
        self.assertEqual(
            result.read_counts, {"00": 5, "01": 0, "10": 0, "11": 5, "x": 10}
        )


class TestGetBreakpoints(unittest.TestCase):
    def test_empty_reads_list(self):
        result = get_spans("chrX", [])
        self.assertEqual(result, [])

    def test_insufficient_reads_below_threshold(self):
        reads = [
            InversionRead("00", None, 100, 200, "read_00")
        ] * 5  # Only 5 reads, below default threshold of 40,
        result = get_spans("chrX", reads, min_fraction=0.1, min_absolute=40)
        self.assertEqual(result, [])

    def test_sufficient_reads_by_fraction(self):
        # 100 reads, 10% = 10 reads threshold
        reads = [InversionRead("00", None, 100, 200, "read_00")] * 50 + [
            InversionRead("01", None, 300, 400, "read_01")
        ] * 50
        result = get_spans("chrX", reads, min_fraction=0.1, min_absolute=5)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0].count, 50)
        self.assertEqual(result[1].count, 50)

    def test_sufficient_reads_by_absolute(self):
        reads = [
            InversionRead("00", None, 100, 200, "read_00")
        ] * 45  # 45 reads meets absolute threshold of 40
        result = get_spans("chrX", reads, min_fraction=0.9, min_absolute=40)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].count, 45)
        self.assertEqual(result[0].region, "chrX:100-200")

    def test_overlapping_reads_same_label_haplotype(self):
        reads = (
            [InversionRead("00", 1, 100, 150, "read_00")] * 50  # Overlaps with next
            + [InversionRead("00", 1, 140, 200, "read_00")]
            * 50  # Overlaps with previous
        )
        result = get_spans("chrX", reads, min_fraction=0.1, min_absolute=40)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].count, 100)  # Merged count
        self.assertEqual(result[0].region, "chrX:100-200")  # Merged range
        self.assertEqual(result[0].inversion_label, "00")
        self.assertEqual(result[0].haplotype, 1)

    def test_overlapping_reads_different_labels_not_merged(self):
        reads = (
            [InversionRead("00", None, 100, 150, "read_00")] * 50
            + [InversionRead("01", None, 140, 200, "read_01")]
            * 50  # Different label, should not merge
        )
        result = get_spans("chrX", reads, min_fraction=0.1, min_absolute=40)
        self.assertEqual(len(result), 2)
        # Should be sorted by start position
        self.assertEqual(result[0].inversion_label, "00")
        self.assertEqual(result[0].region, "chrX:100-150")
        self.assertEqual(result[1].inversion_label, "01")
        self.assertEqual(result[1].region, "chrX:140-200")

    def test_overlapping_reads_different_haplotypes_not_merged(self):
        reads = (
            [InversionRead("00", 1, 100, 150, "read_00")] * 50
            + [InversionRead("00", 2, 140, 200, "read_00")]
            * 50  # Different haplotype, should not merge
        )
        result = get_spans("chrX", reads, min_fraction=0.1, min_absolute=40)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0].haplotype, 1)
        self.assertEqual(result[1].haplotype, 2)

    def test_non_overlapping_reads_not_merged(self):
        reads = (
            [InversionRead("00", None, 100, 150, "read_00")] * 50
            + [InversionRead("00", None, 200, 250, "read_00")] * 50  # No overlap
        )
        result = get_spans("chrX", reads, min_fraction=0.1, min_absolute=40)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0].region, "chrX:100-150")
        self.assertEqual(result[1].region, "chrX:200-250")

    def test_multiple_overlapping_groups(self):
        reads = (
            [InversionRead("00", None, 100, 150, "read_00")] * 50  # Group 1
            + [InversionRead("00", None, 140, 180, "read_00")]
            * 50  # Overlaps with group 1
            + [InversionRead("01", None, 300, 350, "read_01")]
            * 50  # Group 2 (different label)
            + [InversionRead("01", None, 340, 380, "read_01")]
            * 50  # Overlaps with group 2
        )
        result = get_spans("chrX", reads, min_fraction=0.1, min_absolute=40)
        self.assertEqual(len(result), 2)
        # First group (label "00")
        self.assertEqual(result[0].inversion_label, "00")
        self.assertEqual(result[0].region, "chrX:100-180")  # Merged range
        self.assertEqual(result[0].count, 100)
        # Second group (label "01")
        self.assertEqual(result[1].inversion_label, "01")
        self.assertEqual(result[1].region, "chrX:300-380")  # Merged range
        self.assertEqual(result[1].count, 100)

    def test_chain_of_overlapping_reads(self):
        reads = (
            [InversionRead("00", None, 100, 150, "read_00")] * 50  # Overlaps with next
            + [InversionRead("00", None, 140, 190, "read_00")]
            * 50  # Overlaps with prev and next
            + [InversionRead("00", None, 180, 230, "read_00")]
            * 50  # Overlaps with prev
        )
        result = get_spans("chrX", reads, min_fraction=0.1, min_absolute=40)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].region, "chrX:100-230")  # Full merged range
        self.assertEqual(result[0].count, 150)  # All counts combined

    def test_sorting_by_genomic_position(self):
        reads = (
            [InversionRead("01", None, 300, 400, "read_01")] * 50  # Later position
            + [InversionRead("00", None, 100, 200, "read_00")] * 50  # Earlier position
        )
        result = get_spans("chrX", reads, min_fraction=0.1, min_absolute=40)
        self.assertEqual(len(result), 2)
        # Should be sorted by start position
        self.assertEqual(result[0].region, "chrX:100-200")  # Earlier position first
        self.assertEqual(result[1].region, "chrX:300-400")  # Later position second

    def test_edge_case_exact_threshold(self):
        reads = [
            InversionRead("00", None, 100, 200, "read_00")
        ] * 40  # Exactly meets min_absolute=40
        result = get_spans("chrX", reads, min_fraction=0.9, min_absolute=40)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].count, 40)


if __name__ == "__main__":
    unittest.main()
