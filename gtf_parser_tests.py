#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

import gtf_parser

GTF_FP = Path(sys.argv[1])

if not GTF_FP.exists():
    raise ValueError(f"Could not find {GTF_FP} in this dir")
    sys.exit()

genes_of_interest = ["FBLIM1", "RERE", "SMARCAD1", "ESYT2", "NGLY1"]

features = gtf_parser.parse_gtf(GTF_FP, genes_of_interest)

position_lookups = {}

for feat in features:
    gtf_parser.process_transcript(feat, features, position_lookups)


class GTFParserTests(unittest.TestCase):
    def test_correct_output_0(self):
        position = 473
        expected = sorted([("chr1", 8365842), ("chr1", 8365841), ("chr1", 8365840)])

        self.assertEqual(sorted(position_lookups["RERE"][position]), expected)

    def test_correct_output_1(self):
        position = 272
        expected = sorted([("chr1", 15774720), ("chr1", 15774721), ("chr1", 15774722)])

        self.assertEqual(sorted(position_lookups["FBLIM1"][position]), expected)

    def test_correct_output_2(self):
        position = 247
        expected = sorted([("chr4", 94249687), ("chr4", 94249688), ("chr4", 94249689)])

        self.assertEqual(sorted(position_lookups["SMARCAD1"][position]), expected)

    def test_correct_output_3(self):
        position = 857
        expected = sorted(
            [
                ("chr7", 158734239),
                ("chr7", 158734238),
                ("chr7", 158734237),
                ("chr7", 158734245),
                ("chr7", 158734244),
                ("chr7", 158734243),
            ]
        )

        self.assertEqual(sorted(position_lookups["ESYT2"][position]), expected)

    def test_correct_output_4(self):
        position = 121
        expected = sorted(
            [
                ("chr3", 25764197),
                ("chr3", 25764196),
                ("chr3", 25764195),
                ("chr3", 25764071),
                ("chr3", 25764070),
                ("chr3", 25764069),
            ]
        )

        self.assertEqual(sorted(position_lookups["NGLY1"][position]), expected)


if __name__ == "__main__":
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
