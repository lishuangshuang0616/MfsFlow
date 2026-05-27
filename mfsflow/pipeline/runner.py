"""
Pipeline stage orchestration: loads configuration and executes stages sequentially.

This module provides the main pipeline runner that loads the YAML
configuration, initializes the runtime environment, and executes
the pipeline stages (Filtering, Mapping, Counting, Summarising)
in the correct order based on the user-specified starting stage.
"""

import os
import sys

import yaml

from mfsflow.runtime import PipelineRuntime, PipelineTimer, Tee, run_stage_cmd as run_timed_stage_cmd, log_info
from mfsflow.stages.counting import run_counting_stage
from mfsflow.stages.filtering import run_filtering_stage
from mfsflow.stages.mapping import run_mapping_stage
from mfsflow.stages.statistics import run_statistics_stage


def run_pipeline_stages(yaml_file):
    """Orchestrate pipeline stages from a generated run_config.yaml.
    
    Loads the configuration, initializes the runtime environment, and
    executes pipeline stages sequentially starting from the specified stage.
    
    Args:
        yaml_file (str): Path to the run configuration YAML file.
    """
    log_info(f"Loading config from {yaml_file}...")
    with open(yaml_file, "r") as handle:
        config = yaml.safe_load(handle)

    runtime = PipelineRuntime.from_config(config, yaml_file)

    which_stage = runtime.which_stage
    exec_env = runtime.exec_env
    log_path = runtime.log_path
    timing_path = runtime.timing_path

    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    original_stdout = sys.stdout
    with open(log_path, "a") as run_log:
        sys.stdout = Tee(original_stdout, run_log)
        try:
            log_info(f"Starting Pipeline for project: {runtime.project}")
            log_info(f"Stage: {which_stage}")
            timer = PipelineTimer(timing_path, runtime.project)
            log_info(f"Timing log: {timing_path}")
            log_info(f"Temporary chunk directory: {runtime.tmp_merge_path}")

            def run_stage_cmd(cmd, stage_name, shell=False):
                run_timed_stage_cmd(cmd, stage_name, run_log, exec_env, timer, log_path, shell=shell)

            if which_stage == "Filtering":
                umi_chunks, int_chunks = run_filtering_stage(runtime, timer, run_stage_cmd, run_log)

            if which_stage in ["Filtering", "Mapping"]:
                run_mapping_stage(
                    runtime,
                    run_stage_cmd,
                    umi_chunks=locals().get("umi_chunks"),
                    int_chunks=locals().get("int_chunks"),
                )

            if which_stage in ["Filtering", "Mapping", "Counting"]:
                run_counting_stage(runtime, run_stage_cmd)

            if which_stage in ["Filtering", "Mapping", "Counting", "Summarising"]:
                run_statistics_stage(runtime, run_stage_cmd)

            log_info("Pipeline Finished Successfully.")
        finally:
            sys.stdout = original_stdout
