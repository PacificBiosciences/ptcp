import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List

SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import create_input_json  # noqa: E402


def parse_json_stream(text: str) -> List[Dict]:
    decoder = json.JSONDecoder()
    idx = 0
    objs = []
    while idx < len(text):
        while idx < len(text) and text[idx].isspace():
            idx += 1
        if idx >= len(text):
            break
        obj, end = decoder.raw_decode(text, idx)
        objs.append(obj)
        idx = end
    return objs


class TestCreateInputJson(unittest.TestCase):
    def test_fill_template_chunk_respects_max_files(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            template_path = root / "template.json"
            template_path.write_text("{}")
            sample_sheet = root / "sample_sheet.csv"
            sample_sheet.write_text("bam_id,sex\nsample1,F\nsample2,F\n")

            bam1 = root / "sample1.hifi_reads.bam"
            bam2 = root / "sample1.fail_reads.bam"
            bam3 = root / "sample2.hifi_reads.bam"
            for path in (bam1, bam2, bam3):
                path.touch()

            bam_groups = {"sample1": [bam1, bam2], "sample2": [bam3]}
            template, next_idx = create_input_json.fill_template_chunk(
                template_path,
                sample_sheet,
                bam_groups,
                max_files=2,
                start_sample_idx=0,
            )

            self.assertEqual(next_idx, 1)
            self.assertEqual(template["ptcp.sample_sheet"], str(sample_sheet))
            self.assertEqual(template["ptcp.hifi_reads"], [str(bam1.resolve())])
            self.assertEqual(template["ptcp.fail_reads"], [str(bam2.resolve())])
            self.assertNotIn(str(bam3.resolve()), template["ptcp.hifi_reads"])

    def test_process_chunks_outputs_multiple_jsons(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            template_path = root / "template.json"
            template_path.write_text("{}")
            sample_sheet = root / "sample_sheet.csv"
            sample_sheet.write_text("bam_id,sex\na,F\nb,F\n")

            bam_a = root / "a.hifi_reads.bam"
            bam_b = root / "b.hifi_reads.bam"
            bam_a.touch()
            bam_b.touch()

            bam_groups = {"b": [bam_b], "a": [bam_a]}

            buf = io.StringIO()
            with redirect_stdout(buf):
                create_input_json.process_chunks(
                    template_path,
                    sample_sheet,
                    bam_groups,
                    max_files=1,
                )

            objs = parse_json_stream(buf.getvalue())
            self.assertEqual(len(objs), 2)
            self.assertEqual(objs[0]["ptcp.hifi_reads"], [str(bam_a.resolve())])
            self.assertEqual(objs[0]["ptcp.fail_reads"], [])
            self.assertEqual(objs[1]["ptcp.hifi_reads"], [str(bam_b.resolve())])
            self.assertEqual(objs[1]["ptcp.fail_reads"], [])
            self.assertEqual(objs[0]["ptcp.sample_sheet"], str(sample_sheet))
            self.assertEqual(objs[1]["ptcp.sample_sheet"], str(sample_sheet))

    def test_validate_arguments_missing_paths(self):
        args = SimpleNamespace(
            data=Path("missing_data"),
            sample_sheet=Path("missing_sheet.csv"),
            template=Path("missing_template.json"),
        )
        with self.assertRaises(SystemExit):
            create_input_json.validate_arguments(args)


if __name__ == "__main__":
    unittest.main()
