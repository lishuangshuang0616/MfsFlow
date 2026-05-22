import os

from mhsflow.runtime import remove_path


def run_counting_stage(runtime, run_stage_cmd):
    project = runtime.project
    analysis_dir = runtime.analysis_dir
    yaml_file = runtime.yaml_file
    python_exec = runtime.python_exec
    samtools = runtime.tools.samtools
    resolve_script = runtime.resolve_script
    config = runtime.config

    print(">>> Starting Counting Stage")

    umi_aligned = os.path.join(analysis_dir, f"{project}.filtered.tagged.umi.Aligned.out.bam")
    int_aligned = os.path.join(analysis_dir, f"{project}.filtered.tagged.internal.Aligned.out.bam")
    umi_to_tx = os.path.join(analysis_dir, f"{project}.filtered.tagged.umi.Aligned.toTranscriptome.out.bam")
    int_to_tx = os.path.join(analysis_dir, f"{project}.filtered.tagged.internal.Aligned.toTranscriptome.out.bam")

    featurecounts_cmd = [
        python_exec,
        resolve_script("run_featurecounts.py"),
        yaml_file,
        "--umi_bam",
        umi_aligned,
        "--internal_bam",
        int_aligned,
    ]
    run_stage_cmd(featurecounts_cmd, "FeatureCounts (Python)")

    remove_path(umi_aligned)
    remove_path(int_aligned)
    remove_path(umi_to_tx)
    remove_path(int_to_tx)

    print(">>> Starting DGE Analysis (Python)")
    dge_cmd = [python_exec, resolve_script("dge_analysis.py"), yaml_file, samtools]
    run_stage_cmd(dge_cmd, "dge_analysis.py")

    gene_tagged_bam = os.path.join(analysis_dir, f"{project}.filtered.Aligned.GeneTagged.bam")
    stats_enabled = str(config.get("make_stats", "yes")).lower() in ["yes", "true"]
    if not stats_enabled:
        remove_path(gene_tagged_bam)
