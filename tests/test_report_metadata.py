import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from report import _infer_transcriptome_label, export_deliverables_to_outs


class ReportMetadataTests(unittest.TestCase):
    def test_transcriptome_label_uses_parent_for_star_index_dir(self):
        label = _infer_transcriptome_label({"STAR_index": "/path/to/reference/star"})
        self.assertEqual(label, "reference")

    def test_transcriptome_label_uses_index_basename_when_specific(self):
        label = _infer_transcriptome_label({"STAR_index": "/path/to/GRCh38_2024"})
        self.assertEqual(label, "GRCh38_2024")

    def test_transcriptome_label_prefers_explicit_config(self):
        label = _infer_transcriptome_label({
            "STAR_index": "/path/to/reference/star",
            "transcriptome_name": "custom_v1",
        })
        self.assertEqual(label, "custom_v1")

    def test_export_deliverables_to_outs_moves_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            outdir = root / "XPRESS_PROCESSING"
            outs = root / "outs"
            expr = outdir / "expression"
            stats = outdir / "stats"
            config = outdir / "config"
            mex = expr / "Sample01.inex.umi"
            for path in (mex, stats, config):
                path.mkdir(parents=True)

            (expr / "Sample01.h5ad").write_text("h5ad")
            (mex / "matrix.mtx.gz").write_text("matrix")
            (mex / "features.tsv.gz").write_text("features")
            (mex / "barcodes.tsv.gz").write_text("barcodes")
            (stats / "Sample01.stats.tsv").write_text("stats")
            (stats / "Sample01.read_stats.json").write_text("{}")
            (outdir / "Sample01.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam").write_text("bam")
            (outdir / "Sample01.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam.bai").write_text("bai")
            (config / "run_config.yaml").write_text("project: Sample01")
            (config / "expect_id_barcode.tsv").write_text("wellID\tumi_barcodes\tinternal_barcodes\n")

            export_deliverables_to_outs(outdir, outs, "Sample01")

            self.assertTrue((outs / "expression" / "Sample01.h5ad").exists())
            self.assertTrue((outs / "expression" / "Sample01.inex.umi" / "matrix.mtx.gz").exists())
            self.assertTrue((outs / "stats" / "Sample01.stats.tsv").exists())
            self.assertTrue((outs / "stats" / "Sample01.read_stats.json").exists())
            self.assertTrue((outs / "bam" / "Sample01.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam").exists())
            self.assertTrue((outs / "bam" / "Sample01.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam.bai").exists())
            self.assertTrue((outs / "config" / "run_config.yaml").exists())
            self.assertTrue((outs / "config" / "expect_id_barcode.tsv").exists())
            self.assertFalse((expr / "Sample01.h5ad").exists())
            self.assertFalse(mex.exists())
            self.assertFalse((stats / "Sample01.stats.tsv").exists())
            self.assertFalse((outdir / "Sample01.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam").exists())
            self.assertFalse((config / "run_config.yaml").exists())


if __name__ == "__main__":
    unittest.main()
