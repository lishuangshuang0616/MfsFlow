"""
Standalone scripts executed as subprocesses by the pipeline.

Each script runs in its own Python process via ``subprocess``. They share
utility modules within this package (``path_layout``, ``umi_utils``, etc.)
but do NOT import from the parent ``mfsflow`` package directly.
"""