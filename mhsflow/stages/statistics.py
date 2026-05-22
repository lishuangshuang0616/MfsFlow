import os

from mhsflow.runtime import remove_path


def run_statistics_stage(runtime, run_stage_cmd):
    config = runtime.config
    if str(config.get("make_stats", "yes")).lower() not in ["yes", "true"]:
        return

    print(">>> Starting Statistics Stage")
    stats_cmd = [runtime.python_exec, runtime.resolve_script("generate_stats.py"), runtime.yaml_file]
    run_stage_cmd(stats_cmd, "Stats (Python)")

    gene_tagged_bam = os.path.join(runtime.analysis_dir, f"{runtime.project}.filtered.Aligned.GeneTagged.bam")
    remove_path(gene_tagged_bam)
