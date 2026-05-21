#!/usr/bin/env python3
#-*-coding:utf-8-*- 

import os
import argparse
import yaml
import sys
import subprocess
import shutil
from datetime import datetime
import multiprocessing
import gzip
import math
import glob

# Import constants
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))
from pipeline_config import build_base_config, resolve_samplesheet_barcodes
from path_layout import barcode_dir, config_dir, ensure_layout, logs_dir, outputs_dir, tmp_merge_dir
from run_config import write_run_config
from barcode_discovery import build_expected_records, discover_barcodes, write_expected_tables

MIN_PYTHON = (3, 8)


def require_supported_python():
    if sys.version_info < MIN_PYTHON:
        required = ".".join(map(str, MIN_PYTHON))
        current = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
        raise RuntimeError(f"Python {required}+ is required. Current interpreter: Python {current}")


def make_dir(data):
    out_path = data['out_dir']
    outs_path = outputs_dir(out_path)
    
    if os.path.exists(out_path):
        print(f"Warning: Processing directory '{out_path}' already exists. Resuming/Overwriting analysis.")
    
    ensure_layout(out_path)
    
    print(f"Directory 'XPRESS_PROCESSING' (out_dir) created/verified at: {out_path}")
    print(f"Directory 'outs' created/verified at: {outs_path}")

def create_barcode(data):
    sample_type = data['sample']['sample_type'].lower()
    out_path = data['out_dir'] # XPRESS_PROCESSING
    script_path = data.get('toolkit_directory')
    
    # Custom/External Mode
    if sample_type == 'custom' or sample_type == 'external':
        provided_bc = data['barcodes']['barcode_file']
        print(f"Using custom barcode file: {provided_bc}")
        dest_summary = os.path.join(config_dir(out_path), 'expect_id_barcode.tsv')
        dest_pipe = os.path.join(config_dir(out_path), 'expect_barcode.tsv')
        
        shutil.copy(provided_bc, dest_summary)
        
        # Extract annotated barcodes for cell-barcode selection.
        with open(provided_bc, 'r') as infile, open(dest_pipe, 'w') as outfile:
            for line in infile:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    if parts[0].lower() == 'wellid':
                        continue
                    for barcode in parts[1].split(',') + parts[2].split(','):
                        barcode = barcode.strip()
                        if barcode:
                            outfile.write(barcode + '\n')
        
        # Update barcode file path in config
        data['barcodes']['barcode_file'] = dest_pipe
        return

    if sample_type == 'discover':
        records = build_expected_records(script_path)
        pipe_path, _summary_path = write_expected_tables(records, config_dir(out_path))
        data['barcodes']['barcode_file'] = pipe_path
        return

    sample_id_str = str(data['sample']['sample_id'])
    sample_ids = [s.strip() for s in sample_id_str.split(',')]

    with open(os.path.join(config_dir(out_path), 'expect_barcode.tsv'), 'w') as pipe_file, \
         open(os.path.join(config_dir(out_path), 'expect_id_barcode.tsv'), 'w') as summary_file:
        
        print('\t'.join(['wellID','umi_barcodes','internal_barcodes']), file=summary_file)

        if sample_type == 'manual':
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
    
    # Update barcode file path in config
    data['barcodes']['barcode_file'] = os.path.join(config_dir(out_path), 'expect_barcode.tsv')


def run_barcode_discovery(config, project, analysis_dir):
    records = build_expected_records(config.get('toolkit_directory', '.'))
    bcstats_file = os.path.join(analysis_dir, f"{project}.BCstats.txt")
    report_file = os.path.join(barcode_dir(analysis_dir), f"{project}.barcode_discovery.tsv")
    selected, selected_records = discover_barcodes(bcstats_file, records, report_file)
    pipe_path, summary_path = write_expected_tables(selected_records, config_dir(analysis_dir))
    config['barcodes']['barcode_file'] = pipe_path

    selected_label = ", ".join(
        f"{row['candidate_type']}:{row['candidate_id']}({row['matched_reads']} reads/{row['matched_expected_barcodes']} BCs)"
        for row in selected
    )
    print(f">>> Barcode discovery selected: {selected_label}")
    print(f">>> Barcode discovery report: {report_file}")
    print(f">>> Barcode tables updated: {summary_path}")

def check_file_exists(data):
    files_to_check=[data['reference']['STAR_index'], data['reference']['GTF_file']]
    for f in files_to_check:
        if not f or not os.path.exists(f):
            raise FileNotFoundError(f'Reference file not found: {f}')
    
    # Check sequence files
    fq1 = data['sequence_files']['file1']['name']
    fq2 = data['sequence_files']['file2']['name']
    
    for f in [fq1, fq2]:
        if f:
             for subf in f.split(','):
                 if not os.path.exists(subf.strip()):
                     raise FileNotFoundError(f'Fastq file not found: {subf}')

def process_fq(data):
    # User requested to use original data location, skipping symlinks.
    pass

def run_pipeline_stages(yaml_file):
    """
    Orchestrates the pipeline stages.
    """
    print(f"Loading config from {yaml_file}...")
    with open(yaml_file, 'r') as f:
        config = yaml.safe_load(f)

    project = config['project']
    out_dir = config['out_dir']
    num_threads = int(config['num_threads'])
    which_stage = config['which_Stage']
    python_exec = sys.executable or "python3"
    
    # Executables
    samtools = config.get('samtools_exec', 'samtools')
    pigz = config.get('pigz_exec', 'pigz')
    seqkit = config.get('seqkit_exec', 'seqkit')
    toolkit_dir = config.get('toolkit_directory', '.')

    exec_env = os.environ.copy()
    software_dir = os.path.join(toolkit_dir, 'software')
    if sys.platform.startswith('linux') and os.path.isdir(software_dir):
        exec_env['PATH'] = software_dir + os.pathsep + exec_env.get('PATH', '')
    toolkit_src_dir = os.path.join(toolkit_dir, 'src')
    if os.path.isdir(toolkit_src_dir) and toolkit_src_dir not in sys.path:
        sys.path.insert(0, toolkit_src_dir)
    import pipeline_modules

    def resolve_script(script_name):
        direct = os.path.join(toolkit_dir, script_name)
        if os.path.exists(direct):
            return direct
        in_src = os.path.join(toolkit_dir, 'src', script_name)
        if os.path.exists(in_src):
            return in_src
        raise FileNotFoundError(f"Script not found: {script_name}. Tried: {direct}, {in_src}")
    
    class Tee:
        def __init__(self, *streams):
            self._streams = streams

        def write(self, s):
            for stream in self._streams:
                stream.write(s)

        def flush(self):
            for stream in self._streams:
                stream.flush()

    def estimate_avg_line_len(path, sample_lines=1000):
        opener = gzip.open if path.endswith('.gz') else open
        total = 0
        n = 0
        with opener(path, 'rb') as fh:
            for _ in range(sample_lines):
                line = fh.readline()
                if not line:
                    break
                total += len(line)
                n += 1
        return (total / n) if n else 0.0

    
    # config['out_dir'] is already set to XPRESS_PROCESSING in main()
    analysis_dir = out_dir
    log_path = os.path.join(logs_dir(analysis_dir), 'pipeline.log')
    tmp_merge_path = tmp_merge_dir(analysis_dir)
    
    # Ensure log directory exists
    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    original_stdout = sys.stdout
    with open(log_path, 'a') as run_log:
        sys.stdout = Tee(original_stdout, run_log)
        try:
            print(f"Starting Pipeline for project: {project}")
            print(f"Stage: {which_stage}")

            chunk_suffixes = []

            def run_stage_cmd(cmd, stage_name, shell=False):
                if isinstance(cmd, list) and not shell:
                    cmd_str = " ".join(cmd)
                else:
                    cmd_str = str(cmd)
                res = subprocess.run(cmd, stdout=run_log, stderr=subprocess.STDOUT, shell=shell, env=exec_env)
                if res.returncode != 0:
                    run_log.flush()
                    try:
                        with open(log_path, 'r') as lr:
                            print(f"\n[ERROR] {stage_name} failed (rc={res.returncode}). Last 30 lines of log ({log_path}):\n", file=sys.stderr)
                            print("".join(lr.readlines()[-30:]), file=sys.stderr)
                    except Exception:
                        pass
                    raise RuntimeError(f"{stage_name} failed with exit code {res.returncode}.")

            def remove_path(path):
                if not path:
                    return
                if os.path.exists(path):
                    os.remove(path)
                bai = path + ".bai"
                if os.path.exists(bai):
                    os.remove(bai)

            if which_stage == "Filtering":
                print(">>> Starting Filtering Stage")
                
                f1_str = config.get('sequence_files', {}).get('file1', {}).get('name', '')
                f2_str = config.get('sequence_files', {}).get('file2', {}).get('name', '')
                
                fq1_files = [f.strip() for f in f1_str.split(',')] if f1_str else []
                fq2_files = [f.strip() for f in f2_str.split(',')] if f2_str else []

                if not fq1_files:
                    raise ValueError("No file1 found in YAML configuration.")

                total_size_bytes = 0
                first_fq = fq1_files[0]
                
                for f in fq1_files:
                    if f.endswith('.gz'):
                        total_size_bytes += os.path.getsize(f) * 3
                    else:
                        total_size_bytes += os.path.getsize(f)

                avg_line_len = estimate_avg_line_len(first_fq, sample_lines=1000)
                if avg_line_len <= 0:
                    raise ValueError(f"Failed to estimate average line length for {first_fq}")

                total_lines_est = total_size_bytes / avg_line_len
                
                lines_per_chunk = int(math.ceil(total_lines_est / num_threads))
                rem = lines_per_chunk % 4
                if rem != 0: lines_per_chunk += (4 - rem)
                if lines_per_chunk < 4000: lines_per_chunk = 4000
                
                print(f"Total input estimation: {int(total_lines_est/4)} reads.")
                print(f"Split config: {lines_per_chunk} lines per chunk.")

                pool = multiprocessing.Pool(processes=min(2, num_threads)) 
                results = []
                
                if fq2_files:
                    res = pool.apply_async(
                        pipeline_modules.split_fastq,
                        (fq1_files, num_threads, lines_per_chunk, tmp_merge_path, project, pigz, seqkit, fq2_files),
                    )
                    results.append(res)
                else:
                    res = pool.apply_async(
                        pipeline_modules.split_fastq,
                        (fq1_files, num_threads, lines_per_chunk, tmp_merge_path, project, pigz, seqkit),
                    )
                    results.append(res)

                pool.close()
                pool.join()

                chunk_suffixes = results[0].get() 

                print(">>> Running fqfilter.py on chunks")
                max_reads = config.get('counting_opts', {}).get('max_reads', 0)
                
                processes = []
                # Each fqfilter worker can fan out to pigz/samtools, so cap the
                # number of concurrent chunks to avoid oversubscribing CPUs.
                max_filter_jobs = max(1, min(len(chunk_suffixes), max(1, num_threads // 3)))

                def wait_for_filter_process(proc):
                    proc.wait()
                    if proc.returncode != 0:
                        raise RuntimeError(f"fqfilter failed (rc={proc.returncode}). Check {log_path} for details.")

                for suffix in chunk_suffixes:
                    cmd = [python_exec, resolve_script('fqfilter.py'), yaml_file, samtools, pigz, toolkit_dir, suffix]
                    
                    if max_reads and int(max_reads) > 0:
                        chunk_limit = int(int(max_reads) / len(chunk_suffixes))
                        if chunk_limit < 1: chunk_limit = 1
                        cmd.extend(['--limit', str(chunk_limit)])
                        
                    processes.append(subprocess.Popen(cmd, stdout=run_log, stderr=subprocess.STDOUT, env=exec_env))
                    if len(processes) >= max_filter_jobs:
                        wait_for_filter_process(processes.pop(0))

                for p in processes:
                    wait_for_filter_process(p)

                print(">>> Cleaning up temporary FASTQ chunks...")
                import glob
                cleanup_candidates = (
                    glob.glob(os.path.join(tmp_merge_path, "*.part_*"))
                    + glob.glob(os.path.join(tmp_merge_path, "*.part_*.gz"))
                    + glob.glob(os.path.join(tmp_merge_path, "*.fq.part_*"))
                    + glob.glob(os.path.join(tmp_merge_path, "*.fq.part_*.gz"))
                    + glob.glob(os.path.join(tmp_merge_path, "*.fastq.part_*"))
                    + glob.glob(os.path.join(tmp_merge_path, "*.fastq.part_*.gz"))
                )
                for f in cleanup_candidates:
                    if not os.path.exists(f): continue
                    base = os.path.basename(f)
                    if base.endswith(".bam") or base.endswith(".bai") or base.endswith(".txt"): continue
                    if ".raw.tagged." in base or ".filtered.tagged." in base: continue
                    os.remove(f)

                print(">>> Merging BAM Stats")
                pipeline_modules.merge_bam_stats(tmp_merge_path, project, analysis_dir, yaml_file, samtools)

                bc_bin_table = os.path.join(barcode_dir(analysis_dir), f"{project}.BCbinning.txt")
                expect_id_barcode_file = os.path.join(config_dir(out_dir), 'expect_id_barcode.tsv')

                if config.get('barcode_source') == 'samplesheet_barcode':
                    print(">>> Samplesheet barcode mode: skipping barcode detection/binning")
                    bc_bin_for_correction = os.devnull
                else:
                    if config.get('sample', {}).get('sample_type', '').lower() == 'discover':
                        print(">>> Running Barcode Discovery")
                        run_barcode_discovery(config, project, analysis_dir)
                    print(">>> Running Barcode Detection")
                    run_stage_cmd([python_exec, resolve_script("barcode_detection.py"), yaml_file], "BCdetection")
                    bc_bin_for_correction = bc_bin_table

                if config.get('barcode_source') == 'samplesheet_barcode' or os.path.exists(bc_bin_table):
                    print(">>> Correcting BC Tags")
                    correct_processes = []
                    umi_chunks = []
                    int_chunks = []

                    for suffix in chunk_suffixes:
                        raw_bam = os.path.join(tmp_merge_path, f"{project}{suffix}.raw.tagged.bam")
                        fixed_bam_umi = os.path.join(tmp_merge_path, f"{project}{suffix}.filtered.tagged.umi.bam")
                        fixed_bam_int = os.path.join(tmp_merge_path, f"{project}{suffix}.filtered.tagged.internal.bam")
                        
                        umi_chunks.append(fixed_bam_umi)
                        int_chunks.append(fixed_bam_int)

                        if os.path.exists(fixed_bam_umi): os.remove(fixed_bam_umi) 
                        if os.path.exists(fixed_bam_int): os.remove(fixed_bam_int)
                        
                        cmd_args = [python_exec, resolve_script('correct_BCtag.py'), raw_bam, fixed_bam_umi, fixed_bam_int, bc_bin_for_correction, expect_id_barcode_file]
                        correct_processes.append(subprocess.Popen(cmd_args, stdout=run_log, stderr=subprocess.STDOUT, env=exec_env))

                    for p in correct_processes:
                        p.wait()
                        if p.returncode != 0:
                            raise RuntimeError(f"correct_BCtag failed (rc={p.returncode}). Check {log_path} for details.")

                    for suffix in chunk_suffixes:
                        raw_bam = os.path.join(tmp_merge_path, f"{project}{suffix}.raw.tagged.bam")
                        if os.path.exists(raw_bam): os.remove(raw_bam)
                            
                    print(">>> Skipping physical merge of chunks (will stream to STAR)...")
                    
            if which_stage in ["Filtering", "Mapping"]:
                print(">>> Starting Mapping Stage")
                
                umi_arg = ""
                int_arg = ""
                
                # Helper to find chunks if not in memory (e.g. restarting from Mapping)
                def find_chunks(suffix_pattern):
                    found = glob.glob(os.path.join(tmp_merge_path, suffix_pattern))
                    return sorted(found)

                # Resolve UMI inputs
                if 'umi_chunks' in locals() and umi_chunks:
                    umi_arg = ",".join(umi_chunks)
                else:
                    # Try to find chunks on disk
                    disk_umi_chunks = find_chunks(f"{project}*.filtered.tagged.umi.bam")
                    if disk_umi_chunks:
                        print(f"Found {len(disk_umi_chunks)} UMI chunks on disk.")
                        umi_chunks = disk_umi_chunks # Update local var for cleanup later
                        umi_arg = ",".join(disk_umi_chunks)
                    else:
                        legacy_umi = os.path.join(analysis_dir, f"{project}.filtered.tagged.umi.unmapped.bam")
                        if os.path.exists(legacy_umi): 
                            umi_arg = legacy_umi
                        else:
                            # If we are starting at Mapping, we expect inputs.
                            if which_stage == "Mapping":
                                raise FileNotFoundError(f"Could not find input BAMs for Mapping stage. Checked for chunks in {tmp_merge_path} and merged file {legacy_umi}")

                # Resolve Internal inputs
                if 'int_chunks' in locals() and int_chunks:
                    int_arg = ",".join(int_chunks)
                else:
                    disk_int_chunks = find_chunks(f"{project}*.filtered.tagged.internal.bam")
                    if disk_int_chunks:
                        print(f"Found {len(disk_int_chunks)} Internal chunks on disk.")
                        int_chunks = disk_int_chunks # Update local var for cleanup later
                        int_arg = ",".join(disk_int_chunks)
                    else:
                        legacy_int = os.path.join(analysis_dir, f"{project}.filtered.tagged.internal.unmapped.bam")
                        if os.path.exists(legacy_int): 
                            int_arg = legacy_int

                map_cmd = [python_exec, resolve_script('mapping_analysis.py'), yaml_file, '--umi_bam', umi_arg, '--internal_bam', int_arg]
                expect_id_file = os.path.join(out_dir, 'config', 'expect_id_barcode.tsv')
                map_cmd.extend(['--expect_id_file', expect_id_file])
                run_stage_cmd(map_cmd, "mapping_analysis.py")
                                
                if 'umi_chunks' in locals() and umi_chunks:
                     for f in umi_chunks:
                         if os.path.exists(f): os.remove(f)
                
                if 'int_chunks' in locals() and int_chunks:
                     for f in int_chunks:
                         if os.path.exists(f): os.remove(f)

            if which_stage in ["Filtering", "Mapping", "Counting"]:
                print(">>> Starting Counting Stage")
                
                umi_aligned = os.path.join(analysis_dir, f"{project}.filtered.tagged.umi.Aligned.out.bam")
                int_aligned = os.path.join(analysis_dir, f"{project}.filtered.tagged.internal.Aligned.out.bam")
                umi_to_tx = os.path.join(analysis_dir, f"{project}.filtered.tagged.umi.Aligned.toTranscriptome.out.bam")
                int_to_tx = os.path.join(analysis_dir, f"{project}.filtered.tagged.internal.Aligned.toTranscriptome.out.bam")
                
                featurecounts_cmd = [python_exec, resolve_script('run_featurecounts.py'), yaml_file, '--umi_bam', umi_aligned, '--internal_bam', int_aligned]
                run_stage_cmd(featurecounts_cmd, "FeatureCounts (Python)")

                remove_path(umi_aligned)
                remove_path(int_aligned)
                remove_path(umi_to_tx)
                remove_path(int_to_tx)

                print(">>> Starting DGE Analysis (Python)")
                dge_cmd = [python_exec, resolve_script('dge_analysis.py'), yaml_file, samtools]
                run_stage_cmd(dge_cmd, "dge_analysis.py")

                gene_tagged_bam = os.path.join(analysis_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
                stats_enabled = str(config.get('make_stats', 'yes')).lower() in ['yes', 'true']
                if not stats_enabled:
                    remove_path(gene_tagged_bam)

            if which_stage in ["Filtering", "Mapping", "Counting", "Summarising"]:
                if str(config.get('make_stats', 'yes')).lower() in ['yes', 'true']:
                    print(">>> Starting Statistics Stage")
                    stats_cmd = [python_exec, resolve_script('generate_stats.py'), yaml_file]
                    run_stage_cmd(stats_cmd, "Stats (Python)")
                    gene_tagged_bam = os.path.join(analysis_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
                    remove_path(gene_tagged_bam)

            print("Pipeline Finished Successfully.")
        finally:
            sys.stdout = original_stdout

def main():
    require_supported_python()

    parser = argparse.ArgumentParser(description='Mhsflt Data Analysis Pipeline')
    parser.add_argument('--fastqs', required=True, help='Directory containing input R1/R2 FASTQ files')
    parser.add_argument('--samplesheet', help='CSV samplesheet for equal-length R1/R2 data')
    parser.add_argument('--genomeDir', required=True, help='Reference directory containing star/ and genes/genes.gtf')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--outdir', help='Output directory (default: ./<sample_name>)')
    parser.add_argument('--threads', type=int, default=20, help='Number of threads')
    parser.add_argument('--stage', choices=['Filtering', 'Mapping', 'Counting', 'Summarising'], default='Filtering', help='Analysis stage to start from')
    
    # Mutually exclusive group for sample mode
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--manual', help='Manual sample IDs (comma separated, e.g. "20,21"). Sets sample_type=manual.')
    mode_group.add_argument('--plate', help='Plate ID (e.g. "1"). Sets sample_type=auto.')
    mode_group.add_argument('--expectBarcode', help='Path to custom barcode file. Sets sample_type=custom.')
    mode_group.add_argument('--discoverBarcodes', action='store_true', help='Infer plate/manual barcode set from observed reads before barcode detection.')

    args = parser.parse_args()

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Start analysis for {args.sample}.', flush=True)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    config, samplesheet_records = build_base_config(args, script_dir)

    # Setup Directories
    make_dir(config)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Directories created.', flush=True)

    # Validate Files
    check_file_exists(config)

    # Process Fastq (Link to data dir)
    process_fq(config)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Fastq processed.', flush=True)

    # Create Barcode File
    create_barcode(config)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Barcode files created.', flush=True)

    if config.get('barcode_source') == "samplesheet_barcode":
        expect_id_file = os.path.join(config['out_dir'], 'config', 'expect_id_barcode.tsv')
        config['fastq_groups'] = resolve_samplesheet_barcodes(samplesheet_records, expect_id_file)

    final_yaml_path = write_run_config(config)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Config generated: {final_yaml_path}', flush=True)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Starting Pipeline...', flush=True)

    run_pipeline_stages(final_yaml_path)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} All analysis finished.', flush=True)

    # Generate HTML Report
    try:
        from pathlib import Path
        import sys as _sys
        _sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
        import report as report
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Generating HTML Report...', flush=True)
        report.generate_multi_report(config['project'], config['out_dir'], config)
    except Exception as e:
        print(f"Error generating HTML report: {e}", file=sys.stderr)

if __name__ == '__main__':
    main()
