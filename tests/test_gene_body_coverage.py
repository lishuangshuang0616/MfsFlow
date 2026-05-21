import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from run_featurecounts import load_gene_models


class GeneBodyCoverageTests(unittest.TestCase):
    def test_percentile_bins_keep_strand_orientation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = os.path.join(tmpdir, "genes.gtf")
            with open(gtf, "w") as f:
                f.write('chr1\tT\texon\t1\t100\t.\t+\t.\tgene_id "plus"; gene_name "plus";\n')
                f.write('chr1\tT\texon\t1\t100\t.\t-\t.\tgene_id "minus"; gene_name "minus";\n')

            models = load_gene_models(gtf)

            plus = models["plus"]
            minus = models["minus"]

            self.assertEqual(plus["percentile_coords"][0], 1)
            self.assertEqual(plus["percentile_bins"][0], 0)
            self.assertEqual(plus["percentile_coords"][-1], 100)
            self.assertEqual(plus["percentile_bins"][-1], 99)

            self.assertEqual(minus["percentile_coords"][0], 1)
            self.assertEqual(minus["percentile_bins"][0], 99)
            self.assertEqual(minus["percentile_coords"][-1], 100)
            self.assertEqual(minus["percentile_bins"][-1], 0)


if __name__ == "__main__":
    unittest.main()
