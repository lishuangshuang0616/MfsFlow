import os
import shutil
import sys

_SRC_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "src")
if os.path.isdir(_SRC_DIR) and _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

from barcode_discovery import build_expected_records, discover_barcodes, write_expected_tables
from path_layout import barcode_dir, config_dir, ensure_layout, outputs_dir


def create_output_dirs(config):
    out_path = config["out_dir"]
    outs_path = outputs_dir(out_path)

    if os.path.exists(out_path):
        print(f"Warning: Processing directory '{out_path}' already exists. Resuming/Overwriting analysis.")

    ensure_layout(out_path)

    print(f"Directory 'XPRESS_PROCESSING' (out_dir) created/verified at: {out_path}")
    print(f"Directory 'outs' created/verified at: {outs_path}")


def create_barcode_tables(config):
    sample_type = config["sample"]["sample_type"].lower()
    out_path = config["out_dir"]
    script_path = config.get("toolkit_directory")

    if sample_type in ("custom", "external"):
        provided_bc = config["barcodes"]["barcode_file"]
        print(f"Using custom barcode file: {provided_bc}")
        dest_summary = os.path.join(config_dir(out_path), "expect_id_barcode.tsv")
        dest_pipe = os.path.join(config_dir(out_path), "expect_barcode.tsv")

        shutil.copy(provided_bc, dest_summary)

        with open(provided_bc, "r") as infile, open(dest_pipe, "w") as outfile:
            for line in infile:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    if parts[0].lower() == "wellid":
                        continue
                    for barcode in parts[1].split(",") + parts[2].split(","):
                        barcode = barcode.strip()
                        if barcode:
                            outfile.write(barcode + "\n")

        config["barcodes"]["barcode_file"] = dest_pipe
        return

    if sample_type == "discover":
        records = build_expected_records(script_path)
        pipe_path, _summary_path = write_expected_tables(records, config_dir(out_path))
        config["barcodes"]["barcode_file"] = pipe_path
        return

    sample_ids = [s.strip() for s in str(config["sample"]["sample_id"]).split(",")]

    with open(os.path.join(config_dir(out_path), "expect_barcode.tsv"), "w") as pipe_file, \
         open(os.path.join(config_dir(out_path), "expect_id_barcode.tsv"), "w") as summary_file:

        print("\t".join(["wellID", "umi_barcodes", "internal_barcodes"]), file=summary_file)

        if sample_type == "manual":
            records = build_expected_records(script_path, "manual", sample_ids)
        else:
            records = build_expected_records(script_path, "auto", sample_ids)

        grouped = {}
        for rec in records:
            grouped.setdefault(rec["wellID"], {"umi": [], "internal": []})
            grouped[rec["wellID"]][rec["barcode_type"]].append(rec["barcode"])

        for well_id in sorted(grouped):
            umi_str = ",".join(grouped[well_id]["umi"])
            int_str = ",".join(grouped[well_id]["internal"])
            print(f"{well_id}\t{umi_str}\t{int_str}", file=summary_file)
            for barcode in grouped[well_id]["umi"] + grouped[well_id]["internal"]:
                print(barcode, file=pipe_file)

    config["barcodes"]["barcode_file"] = os.path.join(config_dir(out_path), "expect_barcode.tsv")


def run_barcode_discovery(config, project, analysis_dir):
    records = build_expected_records(config.get("toolkit_directory", "."))
    bcstats_file = os.path.join(analysis_dir, f"{project}.BCstats.txt")
    report_file = os.path.join(barcode_dir(analysis_dir), f"{project}.barcode_discovery.tsv")
    selected, selected_records = discover_barcodes(bcstats_file, records, report_file)
    pipe_path, summary_path = write_expected_tables(selected_records, config_dir(analysis_dir))
    config["barcodes"]["barcode_file"] = pipe_path
    if selected:
        config.setdefault("sample", {})["discovered_sample_type"] = selected[0]["candidate_type"]
        config.setdefault("sample", {})["discovered_sample_ids"] = ",".join(
            str(row["candidate_id"]) for row in selected
        )

    selected_label = ", ".join(
        f"{row['candidate_type']}:{row['candidate_id']}({row['matched_reads']} reads/{row['matched_expected_barcodes']} BCs)"
        for row in selected
    )
    print(f">>> Barcode discovery selected: {selected_label}")
    print(f">>> Barcode discovery report: {report_file}")
    print(f">>> Barcode tables updated: {summary_path}")


def process_fastq_inputs(_config):
    # The pipeline uses original FASTQ locations directly.
    return None
