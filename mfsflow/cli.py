import argparse
import os
import sys

_SRC_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "src")
if os.path.isdir(_SRC_DIR) and _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

FILTERING = "Filtering"
STAGE_ORDER = ("Filtering", "Mapping", "Counting", "Summarising")


def build_parser():
    parser = argparse.ArgumentParser(description="MfsFlow Data Analysis Pipeline")
    parser.add_argument("--fastqs", required=True, help="Directory containing input R1/R2 FASTQ files")
    parser.add_argument("--samplesheet", help="CSV samplesheet for equal-length R1/R2 data")
    parser.add_argument("--genomeDir", required=True, help="Reference directory containing star/ and genes/genes.gtf or genes.gtf.gz")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--outdir", help="Output directory (default: ./<sample_name>)")
    parser.add_argument("--threads", type=int, default=20, help="Number of threads")
    parser.add_argument("--tmpRoot", help="Temporary root for chunk BAM/FASTQ files, e.g. /dev/shm")
    parser.add_argument("--stage", choices=STAGE_ORDER, default=FILTERING, help="Analysis stage to start from")

    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument("--manual", help='Manual sample IDs (comma separated, e.g. "20,21"). Sets sample_type=manual.')
    mode_group.add_argument("--plate", help='Plate ID (e.g. "1"). Sets sample_type=auto.')
    mode_group.add_argument("--expectBarcode", help="Path to custom barcode file. Sets sample_type=custom.")
    mode_group.add_argument("--discoverBarcodes", action="store_true", help="Infer plate/manual barcode set from observed reads before barcode detection.")
    return parser


def generate_report(config):
    try:
        from datetime import datetime
        from pathlib import Path

        from path_layout import logs_dir
        from mfsflow.runtime import PipelineTimer

        if str(Path(__file__).resolve().parents[1] / "src") not in sys.path:
            sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
        import report as report

        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Generating HTML Report...', flush=True)
        timing_path = os.path.join(logs_dir(config["out_dir"]), "pipeline_timing.tsv")
        report_timer = PipelineTimer(timing_path, config["project"])
        with report_timer.section("Report: HTML generation"):
            report.generate_multi_report(config["project"], config["out_dir"], config)
    except Exception as exc:
        print(f"Error generating HTML report: {exc}", file=sys.stderr)


def main(argv=None):
    args = build_parser().parse_args(argv)

    import time
    from datetime import datetime
    from pathlib import Path

    from pipeline_config import build_base_config, resolve_samplesheet_barcodes
    from run_config import write_run_config
    from mfsflow.bootstrap import create_barcode_tables, create_output_dirs, process_fastq_inputs
    from mfsflow.config.validation import require_supported_python, validate_input_files
    from mfsflow.pipeline.runner import run_pipeline_stages
    from mfsflow.runtime import format_duration

    require_supported_python()

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Start analysis for {args.sample}.', flush=True)

    script_dir = str(Path(__file__).resolve().parents[1])
    config, samplesheet_records = build_base_config(args, script_dir)

    create_output_dirs(config)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Directories created.', flush=True)

    validate_input_files(config)

    process_fastq_inputs(config)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Fastq processed.', flush=True)

    create_barcode_tables(config)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Barcode files created.', flush=True)

    if config.get("barcode_source") == "samplesheet_barcode":
        expect_id_file = os.path.join(config["out_dir"], "config", "expect_id_barcode.tsv")
        config["fastq_groups"] = resolve_samplesheet_barcodes(samplesheet_records, expect_id_file)

    final_yaml_path = write_run_config(config)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Config generated: {final_yaml_path}', flush=True)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Starting Pipeline...', flush=True)

    pipeline_start = time.perf_counter()
    run_pipeline_stages(final_yaml_path)
    pipeline_duration = time.perf_counter() - pipeline_start

    print(
        f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} All analysis finished '
        f"({format_duration(pipeline_duration)}).",
        flush=True,
    )

    generate_report(config)


if __name__ == "__main__":
    main()
