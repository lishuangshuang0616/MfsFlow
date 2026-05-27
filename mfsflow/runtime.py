"""
Pipeline runtime environment: logging, timing, tool resolution, and stage execution.

This module provides utilities for managing the pipeline runtime, including
logging functions, execution timing, tool path resolution, and stage command
execution with proper error handling and logging.
"""

import logging
import os
import subprocess
import sys
import time
import hashlib
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime

# Configure module-level logger
logger = logging.getLogger(__name__)
logger.propagate = False  # Prevent duplicate logs via root logger

# Create a formatter that matches the original [timestamp] [LEVEL] message format
_formatter = logging.Formatter(
    fmt="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# Set up stdout handler for INFO level
_stdout_handler = logging.StreamHandler(sys.stdout)
_stdout_handler.setLevel(logging.INFO)
_stdout_handler.setFormatter(_formatter)
logger.addHandler(_stdout_handler)

# Set up stderr handler for ERROR level only
_stderr_handler = logging.StreamHandler(sys.stderr)
_stderr_handler.setLevel(logging.ERROR)
_stderr_handler.setFormatter(_formatter)
logger.addHandler(_stderr_handler)


class _StdoutFilter(logging.Filter):
    """Filter that blocks ERROR+ messages from reaching the stdout handler.

    This prevents duplicate output when both handlers are attached: ERROR
    messages go only to stderr (via _stderr_handler), while INFO/WARNING
    messages go to stdout (via _stdout_handler).
    """
    def filter(self, record):
        return record.levelno < logging.ERROR


_stdout_handler.addFilter(_StdoutFilter())

# Set logger level to INFO to capture all messages
logger.setLevel(logging.INFO)


def log_info(msg):
    """Log an informational message with timestamp to stdout.
    
    This function is kept for backward compatibility. Prefer using the
    module-level logger directly for new code.
    
    Args:
        msg (str): Message to log.
    """
    logger.info(msg)


def log_error(msg):
    """Log an error message with timestamp to stderr.
    
    This function is kept for backward compatibility. Prefer using the
    module-level logger directly for new code.
    
    Args:
        msg (str): Error message to log.
    """
    logger.error(msg)

from mfsflow.path_layout import logs_dir, tmp_merge_dir


def format_duration(seconds):
    """Format a duration in seconds to a human-readable string.
    
    Args:
        seconds (float): Duration in seconds.
        
    Returns:
        str: Formatted duration string (e.g., "1.23s", "2m30.00s", "1h15m30.00s").
    """
    seconds = float(seconds)
    if seconds < 60:
        return f"{seconds:.2f}s"
    minutes, sec = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)}m{sec:05.2f}s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)}h{int(minutes):02d}m{sec:05.2f}s"


class Tee:
    """Tee-like stream multiplexer that writes to multiple streams simultaneously.
    
    This class allows writing to multiple output streams (e.g., stdout and a log file)
    with a single write call.
    """
    def __init__(self, *streams):
        """Initialize the Tee with multiple output streams.
        
        Args:
            *streams: Variable number of file-like objects to write to.
        """
        self._streams = streams

    def write(self, s):
        """Write a string to all streams.
        
        Args:
            s (str): String to write.
        """
        for stream in self._streams:
            stream.write(s)

    def flush(self):
        """Flush all streams."""
        for stream in self._streams:
            stream.flush()


class PipelineTimer:
    """Timer for recording pipeline stage execution times.
    
    This class manages timing records for pipeline stages, writing them
    to a TSV file with timestamps, stage names, status, and durations.
    """
    def __init__(self, timing_path, project):
        """Initialize the pipeline timer.
        
        Args:
            timing_path (str): Path to the timing TSV output file.
            project (str): Project name for timing records.
        """
        self.timing_path = timing_path
        self.project = project
        self._ensure_header()

    def _ensure_header(self):
        """Ensure the timing file exists and has the proper header."""
        os.makedirs(os.path.dirname(self.timing_path), exist_ok=True)
        if not os.path.exists(self.timing_path) or os.path.getsize(self.timing_path) == 0:
            with open(self.timing_path, "w") as handle:
                handle.write("timestamp\tproject\tstage\tstatus\tduration_sec\tduration_human\tdetails\n")

    def record(self, stage, status, duration, details=""):
        """Record a timing entry for a pipeline stage.
        
        Args:
            stage (str): Name of the pipeline stage.
            status (str): Stage status ('ok' or 'failed').
            duration (float): Duration in seconds.
            details (str, optional): Additional details about the stage execution.
        """
        safe_details = str(details or "").replace("\t", " ").replace("\n", " ")
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(self.timing_path, "a") as handle:
            handle.write(
                f"{ts}\t{self.project}\t{stage}\t{status}\t{duration:.3f}\t"
                f"{format_duration(duration)}\t{safe_details}\n"
            )

    @contextmanager
    def section(self, stage, details=""):
        """Context manager for timing a pipeline stage section.
        
        Records timing for the enclosed block, handling both success and failure.
        
        Args:
            stage (str): Name of the pipeline stage.
            details (str, optional): Additional details about the stage execution.
            
        Yields:
            None: Context manager yield point.
        """
        start = time.perf_counter()
        try:
            yield
        except Exception:
            duration = time.perf_counter() - start
            self.record(stage, "failed", duration, details)
            log_error(f"Failed {stage} (Duration: {format_duration(duration)})")
            raise
        else:
            duration = time.perf_counter() - start
            self.record(stage, "ok", duration, details)
            log_info(f"Finished {stage} (Duration: {format_duration(duration)})")


@dataclass
class PipelineTools:
    """Data class holding executable paths for external tools."""
    samtools: str
    pigz: str
    seqkit: str


@dataclass
class PipelineRuntime:
    """Runtime configuration and environment for the pipeline.
    
    This data class encapsulates all runtime parameters needed for pipeline
    execution, including configuration, paths, tools, and environment variables.
    """
    config: dict
    yaml_file: str
    project: str
    out_dir: str
    num_threads: int
    which_stage: str
    python_exec: str
    toolkit_dir: str
    tools: PipelineTools
    exec_env: dict
    analysis_dir: str
    log_path: str
    timing_path: str
    tmp_merge_path: str

    @classmethod
    def from_config(cls, config, yaml_file):
        """Create a PipelineRuntime instance from configuration dictionary.
        
        Args:
            config (dict): Pipeline configuration dictionary.
            yaml_file (str): Path to the YAML configuration file.
            
        Returns:
            PipelineRuntime: Configured runtime instance.
        """
        toolkit_dir = config.get("toolkit_directory")
        if not toolkit_dir:
            # Default to mfsflow package location
            toolkit_dir = os.path.dirname(os.path.abspath(__file__))
        elif not os.path.isabs(toolkit_dir):
            # Relative path: resolve relative to YAML file's directory
            toolkit_dir = os.path.join(os.path.dirname(os.path.abspath(yaml_file)), toolkit_dir)
        exec_env = os.environ.copy()
        software_dir = os.path.join(toolkit_dir, "software")
        if sys.platform.startswith("linux") and os.path.isdir(software_dir):
            exec_env["PATH"] = software_dir + os.pathsep + exec_env.get("PATH", "")

        out_dir = config["out_dir"]
        tmp_merge_path = cls._resolve_tmp_merge_path(config, out_dir)
        os.makedirs(tmp_merge_path, exist_ok=True)
        return cls(
            config=config,
            yaml_file=yaml_file,
            project=config["project"],
            out_dir=out_dir,
            num_threads=int(config["num_threads"]),
            which_stage=config["which_Stage"],
            python_exec=sys.executable or "python3",
            toolkit_dir=toolkit_dir,
            tools=PipelineTools(
                samtools=config.get("samtools_exec", "samtools"),
                pigz=config.get("pigz_exec", "pigz"),
                seqkit=config.get("seqkit_exec", "seqkit"),
            ),
            exec_env=exec_env,
            analysis_dir=out_dir,
            log_path=os.path.join(logs_dir(out_dir), "pipeline.log"),
            timing_path=os.path.join(logs_dir(out_dir), "pipeline_timing.tsv"),
            tmp_merge_path=tmp_merge_path,
        )

    @staticmethod
    def _resolve_tmp_merge_path(config, out_dir):
        """Resolve the temporary merge directory path based on configuration.
        
        Args:
            config (dict): Pipeline configuration dictionary.
            out_dir (str): Output directory path.
            
        Returns:
            str: Resolved temporary merge directory path.
        """
        tmp_root = (config.get("performance_opts", {}) or {}).get("tmp_root")
        if not tmp_root:
            return tmp_merge_dir(out_dir)
        project = str(config.get("project", "sample"))
        safe_project = "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in project)
        out_hash = hashlib.sha1(os.path.abspath(out_dir).encode("utf-8")).hexdigest()[:10]
        return os.path.join(os.path.abspath(tmp_root), f"mfsflow_{safe_project}_{out_hash}", "tmp_merge")

    def resolve_script(self, script_name):
        """Resolve the path to a pipeline script.
        
        Searches for the script in multiple locations: toolkit directory,
        scripts subdirectory, and mfsflow package directory.
        
        Args:
            script_name (str): Name of the script to resolve.
            
        Returns:
            str: Full path to the script.
            
        Raises:
            FileNotFoundError: If the script cannot be found in any location.
        """
        # mfsflow package location (where runtime.py resides)
        mfsflow_pkg = os.path.dirname(os.path.abspath(__file__))
        candidates = [
            os.path.join(self.toolkit_dir, script_name),
            os.path.join(self.toolkit_dir, "scripts", script_name),
            os.path.join(mfsflow_pkg, "scripts", script_name),
            os.path.join(mfsflow_pkg, script_name),
        ]
        for candidate in candidates:
            if os.path.exists(candidate):
                return candidate
        raise FileNotFoundError(f"Script not found: {script_name}. Tried: {', '.join(candidates)}")


def remove_path(path):
    """Remove a file and its associated .bai index file if they exist.
    
    Args:
        path (str): Path to the file to remove.
    """
    if not path:
        return
    if os.path.exists(path):
        os.remove(path)
    bai = path + ".bai"
    if os.path.exists(bai):
        os.remove(bai)


def run_stage_cmd(cmd, stage_name, run_log, exec_env, timer, log_path, shell=False):
    """Execute a pipeline stage command with logging and error handling.
    
    Args:
        cmd (str or list): Command to execute.
        stage_name (str): Name of the pipeline stage for logging.
        run_log (file): File object for command output logging.
        exec_env (dict): Environment variables for command execution.
        timer (PipelineTimer): Timer for recording stage execution.
        log_path (str): Path to the pipeline log file.
        shell (bool, optional): Whether to execute command through shell. Defaults to False.
        
    Raises:
        RuntimeError: If the command fails with non-zero exit code.
    """
    if isinstance(cmd, list) and not shell:
        cmd_str = " ".join(map(str, cmd))
    else:
        cmd_str = str(cmd)
    start = time.perf_counter()
    res = subprocess.run(cmd, stdout=run_log, stderr=subprocess.STDOUT, shell=shell, env=exec_env)
    duration = time.perf_counter() - start
    if res.returncode != 0:
        timer.record(stage_name, "failed", duration, cmd_str)
        log_error(f"Failed {stage_name} (Duration: {format_duration(duration)})")
        run_log.flush()
        try:
            with open(log_path, "r") as lr:
                log_error(f"{stage_name} failed (rc={res.returncode}). Last 30 lines of log ({log_path}):\n" + "".join(lr.readlines()[-30:]))
        except Exception:
            pass
        raise RuntimeError(f"{stage_name} failed with exit code {res.returncode}.")
    timer.record(stage_name, "ok", duration, cmd_str)
    log_info(f"Finished {stage_name} (Duration: {format_duration(duration)})")
