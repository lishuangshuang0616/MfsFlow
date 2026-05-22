import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from run_featurecounts import parse_gtf_and_create_saf


class SafBuilderTests(unittest.TestCase):
    def test_introns_are_clipped_away_from_overlapping_gene_bodies(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = os.path.join(tmpdir, "genes.gtf")
            with open(gtf, "w") as f:
                f.write('chr1\tT\texon\t100\t120\t.\t+\t.\tgene_id "geneA"; gene_name "A";\n')
                f.write('chr1\tT\texon\t200\t220\t.\t+\t.\tgene_id "geneA"; gene_name "A";\n')
                f.write('chr1\tT\texon\t150\t160\t.\t+\t.\tgene_id "geneB"; gene_name "B";\n')
                f.write('chr1\tT\texon\t400\t420\t.\t+\t.\tgene_id "geneB"; gene_name "B";\n')

            saf, gene_map = parse_gtf_and_create_saf(gtf, os.path.join(tmpdir, "sample"))
            self.assertEqual(gene_map["geneA"], "A")
            self.assertTrue(os.path.exists(os.path.join(tmpdir, "sample.exon.saf")))
            self.assertTrue(os.path.exists(os.path.join(tmpdir, "sample.intron.saf")))

            with open(saf) as f:
                lines = [line.strip().split("\t") for line in f if not line.startswith("GeneID")]

            introns = [row for row in lines if row[0] == "geneA__INTRON__"]
            self.assertEqual(
                introns,
                [
                    ["geneA__INTRON__", "chr1", "121", "149", "+"],
                ],
            )

            self.assertNotIn(["geneA__INTRON__", "chr1", "150", "160", "+"], introns)
            self.assertNotIn(["geneA__INTRON__", "chr1", "161", "199", "+"], introns)

    def test_introns_respect_r_length_limits(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = os.path.join(tmpdir, "genes.gtf")
            with open(gtf, "w") as f:
                f.write('chr1\tT\texon\t100\t120\t.\t+\t.\tgene_id "short"; gene_name "short";\n')
                f.write('chr1\tT\texon\t126\t140\t.\t+\t.\tgene_id "short"; gene_name "short";\n')
                f.write('chr1\tT\texon\t1000\t1020\t.\t+\t.\tgene_id "long"; gene_name "long";\n')
                f.write('chr1\tT\texon\t101050\t101070\t.\t+\t.\tgene_id "long"; gene_name "long";\n')

            saf, _gene_map = parse_gtf_and_create_saf(gtf, os.path.join(tmpdir, "sample"))
            with open(saf) as f:
                content = f.read()

            self.assertNotIn("short__INTRON__", content)
            self.assertNotIn("long__INTRON__", content)


if __name__ == "__main__":
    unittest.main()
