# MfsFlow Architecture

MfsFlow is moving from a script-oriented pipeline into a package-oriented
pipeline. The legacy scripts in `src/` are still used as executable stage tools,
while orchestration and shared runtime code now live in `mfsflow/`.

## Entry Points

- `run_analysis_pipeline.py`: compatibility launcher for existing command lines.
- `mfsflow.cli`: canonical Python CLI implementation.
- `mfsflow` console script: installed entry point defined in `pyproject.toml`.

## Package Layout

- `mfsflow.config`: input validation and future typed configuration models.
- `mfsflow.bootstrap`: pre-run setup, output directory creation, barcode table creation.
- `mfsflow.runtime`: runtime context, path derivation, timing, subprocess execution helpers.
- `mfsflow.pipeline`: stage orchestration.
- `mfsflow.stages`: stage-level orchestration for filtering, mapping, counting, and statistics.

## Legacy Tool Modules

The existing `src/` modules still contain the heavy algorithmic implementations:

- `fqfilter.py`
- `mapping_analysis.py`
- `run_featurecounts.py`
- `dge_analysis.py`
- `generate_stats.py`
- barcode and report helpers

These modules should be migrated gradually into package modules after behavior is
locked down by tests. Until then, the pipeline runner resolves and executes them
through the same paths used by the previous command-line workflow.

## Refactoring Rules

Keep each refactor behavior-preserving unless explicitly changing analysis
logic. Move code by responsibility first, then optimize internals:

1. CLI and configuration.
2. Runtime and stage orchestration.
3. Stage implementation modules.
4. Data model and artifact contracts.
5. Algorithm-level performance work.

All stage-level changes should keep `run_analysis_pipeline.py --help`,
`python3 -m unittest discover -s tests`, and `git diff --check` passing.
