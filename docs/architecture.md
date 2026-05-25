# MfsFlow Architecture

MfsFlow is a package-oriented pipeline. All code lives in the `mfsflow/`
package. After `pip install .`, the package is self-contained and requires
no external script directories.

## Entry Points

- `mfsflow.cli`: canonical Python CLI implementation.
- `mfsflow` console script: installed entry point defined in `pyproject.toml`.
- `run_analysis_pipeline.py`: compatibility launcher (deprecated).

## Package Layout

```
mfsflow/
├── cli.py                  # CLI entry point
├── runtime.py              # Runtime context, path derivation, timing
├── bootstrap.py            # Pre-run setup, barcode table creation
├── pipeline_config.py      # YAML config building and barcode resolution
├── run_config.py          # Run configuration writer
├── path_layout.py         # Output directory layout constants
├── constant.py            # Shared constants
├── report.py              # HTML report generation
├── config/                # Input validation and typed configuration models
├── pipeline/              # Stage orchestration
├── stages/                # Stage implementations (filtering, mapping, counting, statistics)
├── scripts/               # Executable stage tools (fqfilter, mapping_analysis, etc.)
├── software/              # Bundled bioinformatics tools (STAR, samtools, etc.)
├── yaml/                  # Example config and barcode list YAML files
└── report_assets/         # HTML report templates and Plotly.js
```

## Key Modules

- `mfsflow.config`: input validation and typed configuration models.
- `mfsflow.bootstrap`: pre-run setup, output directory creation, barcode table creation.
- `mfsflow.runtime`: runtime context, path derivation, timing, subprocess execution helpers.
- `mfsflow.pipeline`: stage orchestration.
- `mfsflow.stages`: stage-level orchestration for filtering, mapping, counting, and statistics.
- `mfsflow.scripts`: executable stage tools (previously in `src/`).

## Bundled Tools

Bioinformatics tools are bundled in `mfsflow/software/` and automatically
added to `PATH` at runtime (Linux only). No manual installation required:

- `STAR` — splice-aware aligner
- `samtools` — BAM manipulation
- `featureCounts` — read quantification
- `pigz` — parallel gzip
- `seqkit` — FASTQ processing

## Refactoring Rules

Keep each refactor behavior-preserving unless explicitly changing analysis
logic. Move code by responsibility first, then optimize internals:

1. CLI and configuration.
2. Runtime and stage orchestration.
3. Stage implementation modules.
4. Data model and artifact contracts.
5. Algorithm-level performance work.

All stage-level changes should keep `mfsflow --help`,
`python3 -m unittest discover -s tests`, and `git diff --check` passing.
