import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from run_featurecounts import _project_blocks_to_gene_body, _pysam_blocks_1based_half_open, load_gene_models


class GeneBodyCoverageTests(unittest.TestCase):
    def test_transcript_body_projection_keeps_strand_orientation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = os.path.join(tmpdir, "genes.gtf")
            with open(gtf, "w") as f:
                f.write('chr1\tT\texon\t1\t100\t.\t+\t.\tgene_id "plus"; transcript_id "plus_t1"; gene_name "plus";\n')
                f.write('chr1\tT\texon\t1\t100\t.\t-\t.\tgene_id "minus"; transcript_id "minus_t1"; gene_name "minus";\n')

            models = load_gene_models(gtf)

            plus = models["plus"]
            minus = models["minus"]

            plus_overlap, plus_bins = _project_blocks_to_gene_body(plus, [(1, 11)])
            minus_overlap, minus_bins = _project_blocks_to_gene_body(minus, [(91, 101)])

            self.assertEqual(10, plus_overlap)
            self.assertEqual(10, minus_overlap)
            self.assertEqual([1.0] * 10, plus_bins[:10])
            self.assertEqual([1.0] * 10, minus_bins[:10])
            self.assertEqual(0.0, plus_bins[10])
            self.assertEqual(0.0, minus_bins[10])

    def test_longest_transcript_is_selected_per_gene(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = os.path.join(tmpdir, "genes.gtf")
            with open(gtf, "w") as f:
                f.write('chr1\tT\texon\t1\t100\t.\t+\t.\tgene_id "gene1"; transcript_id "short";\n')
                f.write('chr1\tT\texon\t1\t150\t.\t+\t.\tgene_id "gene1"; transcript_id "long";\n')

            models = load_gene_models(gtf)
            self.assertEqual("long", models["gene1"]["transcript_id"])
            self.assertEqual(150, models["gene1"]["length"])

    def test_pysam_blocks_are_converted_to_1based_half_open(self):
        class DummyRead:
            def get_blocks(self):
                return [(0, 100), (150, 175)]

        self.assertEqual(
            [(1, 101), (151, 176)],
            _pysam_blocks_1based_half_open(DummyRead()),
        )


if __name__ == "__main__":
    unittest.main()
