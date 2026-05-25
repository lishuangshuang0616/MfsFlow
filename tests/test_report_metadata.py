import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from report import _infer_transcriptome_label


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


if __name__ == "__main__":
    unittest.main()
