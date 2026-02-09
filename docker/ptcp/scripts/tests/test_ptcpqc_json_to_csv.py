import argparse
import json
import math
import sys
import unittest
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from tempfile import TemporaryDirectory


def load_ptcpqc_json_to_csv_module():
    script_path = Path(__file__).resolve().parents[1] / "ptcpqc_json_to_csv.py"
    spec = spec_from_file_location("ptcpqc_json_to_csv", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module from {script_path}")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


PTCPQC_JSON_TO_CSV = load_ptcpqc_json_to_csv_module()


def make_header(report_type=None):
    header = {
        "genome_version": "hg38",
        "targets_bed": "targets.bed",
        "ptcp-qc_version": "0.0.0",
        "timestamp": "2026-01-01T00:00:00Z",
    }
    if report_type is not None:
        header["report_type"] = report_type
    return header


def write_json(path, payload):
    path.write_text(json.dumps(payload), encoding="utf-8")


def make_sample_report(sample_name):
    return {
        "header": make_header("sample"),
        "sample_name": sample_name,
        "stats": {"total_reads": 10},
        "locus_results": {
            "L1": {
                "trgt_results": {
                    "HTT": {
                        "repeat_unit": "CAG",
                        "alleles": {
                            "1": {
                                "coverage": 5,
                                "size": 15,
                                "ci": {"min": 14, "max": 16},
                            }
                        },
                    }
                },
                "paraphase_results": {
                    "HTT": {
                        "total_reads": 10,
                        "total_cn": 2,
                        "genotype": "gt",
                    }
                },
                "reads": {"n": 10, "hifi": 8},
            }
        },
    }


def make_aggregate_report(sample_name):
    return {
        "header": make_header("aggregate"),
        "sample_stats": [
            {
                "sample_name": sample_name,
                "total_reads": 100,
                "coverage": {
                    "mean_coverage": 12.345,
                    "pct_targets_ge_10": 0.5,
                    "pct_targets_ge_20": 0.4,
                    "pct_targets_ge_30": 0.3,
                    "pct_targets_lt_5": 0.1,
                },
                "on_target_fraction": 0.125,
                "duplicate_fraction": 0.2,
                "read_length_stats": {"median": 1234.4},
                "read_quality_stats": {"median": 99},
            }
        ],
    }


class TestFormatters(unittest.TestCase):
    def test_accuracy_to_q(self):
        self.assertIsNone(PTCPQC_JSON_TO_CSV.accuracy_to_q(-1.0))
        self.assertTrue(math.isinf(PTCPQC_JSON_TO_CSV.accuracy_to_q(1.0)))
        self.assertAlmostEqual(PTCPQC_JSON_TO_CSV.accuracy_to_q(0.9), 10.0)
        with self.assertRaises(ValueError):
            PTCPQC_JSON_TO_CSV.accuracy_to_q(1.1)

    def test_formatters_basic(self):
        self.assertEqual(PTCPQC_JSON_TO_CSV.format_int(None, ""), "")
        self.assertEqual(PTCPQC_JSON_TO_CSV.format_int(5.7, ""), "6")
        self.assertEqual(PTCPQC_JSON_TO_CSV.format_int(4.5, ""), "5")
        self.assertEqual(PTCPQC_JSON_TO_CSV.format_int(1234.4, ""), "1234")
        self.assertEqual(PTCPQC_JSON_TO_CSV.format_quality(99), "Q20")
        self.assertEqual(PTCPQC_JSON_TO_CSV.format_quality(100), "Q60")
        self.assertEqual(PTCPQC_JSON_TO_CSV.format_float_trim(1.234, ""), "1.23")
        self.assertEqual(PTCPQC_JSON_TO_CSV.format_percentage(0.125), "12.50%")
        self.assertEqual(PTCPQC_JSON_TO_CSV.format_int(-3, ""), "")
        self.assertEqual(PTCPQC_JSON_TO_CSV.format_int(2.6, ""), "3")


class TestHeaderParsing(unittest.TestCase):
    def test_extract_header_requires_fields(self):
        data = {"header": {"report_type": "sample"}}
        with self.assertRaises(ValueError):
            PTCPQC_JSON_TO_CSV.extract_header(data, Path("sample.json"))

    def test_extract_header_allows_missing_report_type(self):
        header = PTCPQC_JSON_TO_CSV.extract_header(make_header(), Path("sample.json"))
        self.assertEqual(header["genome_version"], "hg38")

    def test_determine_report_type(self):
        path = Path("report.json")
        header = {"report_type": "sample"}
        data = {"locus_results": {}}
        self.assertEqual(
            PTCPQC_JSON_TO_CSV.determine_report_type(data, header, path), "sample"
        )

        header = {}
        data = {"report_type": "aggregate"}
        self.assertEqual(
            PTCPQC_JSON_TO_CSV.determine_report_type(data, header, path), "aggregate"
        )

        header = {}
        data = {"sample_stats": []}
        self.assertEqual(
            PTCPQC_JSON_TO_CSV.determine_report_type(data, header, path), "aggregate"
        )

        header = {}
        data = {"locus_results": {}}
        self.assertEqual(
            PTCPQC_JSON_TO_CSV.determine_report_type(data, header, path), "sample"
        )

        with self.assertRaises(ValueError):
            PTCPQC_JSON_TO_CSV.determine_report_type({}, {}, path)


class TestLoadReport(unittest.TestCase):
    def test_load_report_sets_report_type(self):
        with TemporaryDirectory() as tempdir:
            path = Path(tempdir) / "report.json"
            payload = {"header": make_header(), "sample_stats": []}
            write_json(path, payload)
            report = PTCPQC_JSON_TO_CSV.load_report(path)

        self.assertEqual(report.report_type, "aggregate")
        self.assertEqual(report.header["report_type"], "aggregate")


class TestAggregateTables(unittest.TestCase):
    def test_table_from_aggregate_reports(self):
        data = make_aggregate_report("sample1")
        report = PTCPQC_JSON_TO_CSV.ParsedReport(
            path=Path("aggregate.json"),
            data=data,
            header=data["header"],
            report_type="aggregate",
        )
        table = PTCPQC_JSON_TO_CSV.table_from_aggregate_and_sample_reports(
            report,
            [],
            ["sample1"],
            lambda sample_name: sample_name,
            PTCPQC_JSON_TO_CSV.MissingPolicy(),
        )

        self.assertEqual(table.columns, PTCPQC_JSON_TO_CSV.AGGREGATE_OUTPUT_COLUMNS)
        self.assertEqual(table.rows[0][0], "sample1")
        self.assertEqual(
            table.rows[0][table.columns.index("Median Mapped Read Length")], "1234"
        )
        self.assertEqual(
            table.rows[0][table.columns.index("Median Mapped Read Quality")], "Q20"
        )
        self.assertEqual(
            table.rows[0][table.columns.index("Mean Target Coverage")], "12.35"
        )
        self.assertEqual(
            table.rows[0][table.columns.index("Percent of Targets >=10-fold Coverage")],
            "50.00%",
        )
        self.assertEqual(
            table.rows[0][table.columns.index("Percent of On-Target Reads")], "12.50%"
        )
        self.assertEqual(
            table.rows[0][table.columns.index("Percent of Duplicate Reads")], "20.00%"
        )
        self.assertEqual(table.rows[0][table.columns.index("on_target_total")], "")


class TestTrgtTables(unittest.TestCase):
    def test_table_from_trgt_reports_coverage_precedence(self):
        data = {
            "sample_name": "sample1",
            "locus_results": {
                "group1": {
                    "trgt_results": {
                        "HTT": {
                            "repeat_unit": "CAG",
                            "alleles": {
                                "0": {
                                    "coverage": 4,
                                    "size": 20,
                                    "ci": {"min": 19, "max": 21},
                                    "motif_count": "mc0",
                                    "motif_spans": "ms0",
                                },
                            },
                            "read_stats": {"0": {"coverage": 9}},
                        }
                    }
                }
            }
        }

        report = PTCPQC_JSON_TO_CSV.ParsedReport(
            path=Path("sample.json"),
            data=data,
            header=make_header("sample"),
            report_type="sample",
        )
        table = PTCPQC_JSON_TO_CSV.table_from_trgt_reports(
            [report],
            ["sample1"],
            lambda sample_name: sample_name,
            PTCPQC_JSON_TO_CSV.MissingPolicy(),
        )

        self.assertEqual(table.columns, PTCPQC_JSON_TO_CSV.TRGT_OUTPUT_COLUMNS)
        self.assertEqual(
            table.rows[0],
            ["sample1", "HTT", "0", "9", "20", "19", "21", "CAG", "mc0", "ms0"],
        )
        self.assertEqual(
            table.rows[1],
            ["sample1", "HTT", "1", "", "", "", "", "CAG", "", ""],
        )

class TestCoverageTables(unittest.TestCase):
    def test_table_from_coverage_reports_subgroups(self):
        data = {
            "sample_name": "sample1",
            "locus_results": {
                "L1": {
                    "reads": {
                        "n": 5,
                        "hifi": 4,
                        "stats": {
                            "read_quality": {"median": 99},
                            "read_passes": {"median": 12.4},
                            "read_length": {"median": 1000.6},
                        },
                    },
                    "subgroups": {"subA": {"reads": {"n": 2, "hifi": 1}}},
                }
            },
        }

        report = PTCPQC_JSON_TO_CSV.ParsedReport(
            path=Path("sample.json"),
            data=data,
            header=make_header("sample"),
            report_type="sample",
        )
        table = PTCPQC_JSON_TO_CSV.table_from_coverage_reports(
            [report],
            ["sample1"],
            lambda sample_name: sample_name,
            PTCPQC_JSON_TO_CSV.MissingPolicy(),
        )

        self.assertEqual(table.columns, PTCPQC_JSON_TO_CSV.COVERAGE_OUTPUT_COLUMNS)
        self.assertEqual(
            table.rows[0],
            ["sample1", "L1", "L1_total", "5", "4", "99", "12.4", "1000.6"],
        )
        self.assertEqual(
            table.rows[1],
            ["sample1", "L1", "subA", "2", "1", "", "", ""],
        )


class TestSampleTables(unittest.TestCase):
    def test_table_from_paraphase_reports_alignment_and_sorting(self):
        data1 = {
            "sample_name": "sample1",
            "locus_results": {
                "group": {"paraphase_results": {"B": {"total_reads": 2, "total_cn": 2}}}
            },
        }
        data2 = {
            "sample_name": "sample2",
            "locus_results": {
                "group": {"paraphase_results": {"A": {"total_reads": 3, "total_cn": 3}}}
            },
        }

        report1 = PTCPQC_JSON_TO_CSV.ParsedReport(
            path=Path("sample1.json"),
            data=data1,
            header=make_header("sample"),
            report_type="sample",
        )
        report2 = PTCPQC_JSON_TO_CSV.ParsedReport(
            path=Path("sample2.json"),
            data=data2,
            header=make_header("sample"),
            report_type="sample",
        )

        table = PTCPQC_JSON_TO_CSV.table_from_paraphase_reports(
            [report1, report2],
            ["sample1", "sample2"],
            lambda sample_name: sample_name,
            PTCPQC_JSON_TO_CSV.MissingPolicy(),
        )
        self.assertEqual(table.columns, PTCPQC_JSON_TO_CSV.PARAPHASE_OUTPUT_COLUMNS)

        self.assertEqual(table.rows[0][0:2], ["sample1", "A"])
        self.assertEqual(table.rows[0][2], "")
        self.assertEqual(table.rows[1][0:2], ["sample2", "A"])
        self.assertEqual(table.rows[1][2], "3")
        self.assertEqual(table.rows[2][0:2], ["sample1", "B"])
        self.assertEqual(table.rows[2][2], "2")
        self.assertEqual(table.rows[3][0:2], ["sample2", "B"])
        self.assertEqual(table.rows[3][2], "")

    def test_table_from_paraphase_reports_formats_sv_calls_and_inversion(self):
        data = {
            "sample_name": "sample1",
            "locus_results": {
                "group1": {
                    "paraphase_results": {
                        "F8": {"total_cn": 2, "f8_info": {"has_inversion": True}},
                        "GENE": {
                            "total_cn": 2,
                            "genotype_adjusted": "g1_adj",
                            "sv_calls": {"dup": {}, "del": {}},
                        },
                    },
                    "sawfish_results": {
                        "GENE": {
                            "sv_breakpoints": [
                                {
                                    "chrom": "1",
                                    "start": 100,
                                    "end": 200,
                                    "svtype": "DEL",
                                    "annotation": {"name": "annot"},
                                },
                                {"chrom": "1", "start": 300, "end": 400, "svtype": "DUP"},
                            ]
                        }
                    },
                }
            },
        }

        report = PTCPQC_JSON_TO_CSV.ParsedReport(
            path=Path("sample.json"),
            data=data,
            header=make_header("sample"),
            report_type="sample",
        )
        table = PTCPQC_JSON_TO_CSV.table_from_paraphase_reports(
            [report],
            ["sample1"],
            lambda sample_name: sample_name,
            PTCPQC_JSON_TO_CSV.MissingPolicy(numeric="NA", string="NA"),
        )

        self.assertEqual(table.columns, PTCPQC_JSON_TO_CSV.PARAPHASE_OUTPUT_COLUMNS)
        self.assertEqual(
            table.rows[0],
            ["sample1", "F8", "2", "NA", "TRUE", "NA", "", ""],
        )
        self.assertEqual(
            table.rows[1],
            [
                "sample1",
                "GENE",
                "2",
                "NA",
                "",
                "g1_adj",
                "del,dup",
                "1:100:200:DEL:annot,1:300:400:DUP",
            ],
        )


class TestRunParser(unittest.TestCase):
    def test_run_parser_writes_sample_outputs(self):
        with TemporaryDirectory() as tempdir:
            temp_path = Path(tempdir)
            sample_path = temp_path / "sample.json"
            write_json(sample_path, make_sample_report("Sample1"))

            args = argparse.Namespace(
                inputs=[sample_path],
                input_path=None,
                output=temp_path / "out.csv",
                output_types=["all"],
                verbose=0,
                missing_as_na=False,
                anon=False,
            )
            PTCPQC_JSON_TO_CSV.run_parser(args)

            trgt_path = temp_path / "out.trgt.csv"
            paraphase_path = temp_path / "out.paraphase.csv"
            coverage_path = temp_path / "out.coverage.csv"

            self.assertTrue(trgt_path.exists())
            self.assertTrue(paraphase_path.exists())
            self.assertTrue(coverage_path.exists())
            self.assertEqual(
                trgt_path.read_text().splitlines()[0],
                ",".join(PTCPQC_JSON_TO_CSV.TRGT_OUTPUT_COLUMNS),
            )
            self.assertEqual(
                paraphase_path.read_text().splitlines()[0],
                ",".join(PTCPQC_JSON_TO_CSV.PARAPHASE_OUTPUT_COLUMNS),
            )
            self.assertEqual(
                coverage_path.read_text().splitlines()[0],
                ",".join(PTCPQC_JSON_TO_CSV.COVERAGE_OUTPUT_COLUMNS),
            )

    def test_run_parser_writes_aggregate_output(self):
        with TemporaryDirectory() as tempdir:
            temp_path = Path(tempdir)
            agg_path = temp_path / "aggregate.json"
            write_json(agg_path, make_aggregate_report("Sample1"))

            args = argparse.Namespace(
                inputs=[agg_path],
                input_path=None,
                output=temp_path / "out.csv",
                output_types=["all"],
                verbose=0,
                missing_as_na=False,
                anon=False,
            )
            PTCPQC_JSON_TO_CSV.run_parser(args)

            output_path = temp_path / "out.aggregate.csv"
            self.assertTrue(output_path.exists())
            header = output_path.read_text().splitlines()[0]
            self.assertEqual(
                header, ",".join(PTCPQC_JSON_TO_CSV.AGGREGATE_OUTPUT_COLUMNS)
            )

    def test_run_parser_anonymizes_sample_names(self):
        with TemporaryDirectory() as tempdir:
            temp_path = Path(tempdir)
            sample_a = temp_path / "sample_a.json"
            sample_b = temp_path / "sample_b.json"
            write_json(sample_a, make_sample_report("Alpha"))
            write_json(sample_b, make_sample_report("Beta"))

            args = argparse.Namespace(
                inputs=[sample_a, sample_b],
                input_path=None,
                output=temp_path / "out.csv",
                output_types=["trgt"],
                verbose=0,
                missing_as_na=False,
                anon=True,
            )
            PTCPQC_JSON_TO_CSV.run_parser(args)

            output_path = temp_path / "out.trgt.csv"
            contents = output_path.read_text()
            self.assertNotIn("Alpha", contents)
            self.assertNotIn("Beta", contents)
            self.assertIn("Sample_1", contents)
            self.assertIn("Sample_2", contents)

    def test_run_parser_missing_input_exits(self):
        with TemporaryDirectory() as tempdir:
            temp_path = Path(tempdir)
            missing = temp_path / "missing.json"
            args = argparse.Namespace(
                inputs=[missing],
                input_path=None,
                output=temp_path / "out.csv",
                output_types=["all"],
                verbose=0,
                missing_as_na=False,
                anon=False,
            )
            with self.assertRaises(SystemExit) as context:
                PTCPQC_JSON_TO_CSV.run_parser(args)

        self.assertEqual(context.exception.code, 1)

    def test_missing_as_na_controls_numeric_missing(self):
        with TemporaryDirectory() as tempdir:
            temp_path = Path(tempdir)

            sample_a = temp_path / "sample_a.json"
            sample_b = temp_path / "sample_b.json"

            write_json(
                sample_a,
                {
                    "header": make_header("sample"),
                    "sample_name": "sample_a",
                    "locus_results": {"A": {"reads": {"n": 1, "hifi": 1}}},
                },
            )
            write_json(
                sample_b,
                {
                    "header": make_header("sample"),
                    "sample_name": "sample_b",
                    "locus_results": {"B": {"reads": {"n": 2, "hifi": 2}}},
                },
            )

            args = argparse.Namespace(
                inputs=[sample_a, sample_b],
                input_path=None,
                output=temp_path / "out.csv",
                output_types=["coverage"],
                verbose=0,
                missing_as_na=False,
                anon=False,
            )
            PTCPQC_JSON_TO_CSV.run_parser(args)

            coverage_path = temp_path / "out.coverage.csv"
            contents = coverage_path.read_text()
            self.assertNotIn("NA", contents)

            args = argparse.Namespace(
                inputs=[sample_a, sample_b],
                input_path=None,
                output=temp_path / "out_na.csv",
                output_types=["coverage"],
                verbose=0,
                missing_as_na=True,
                anon=False,
            )
            PTCPQC_JSON_TO_CSV.run_parser(args)

            coverage_path = temp_path / "out_na.coverage.csv"
            contents = coverage_path.read_text()
            self.assertIn("NA", contents)


if __name__ == "__main__":
    unittest.main()
