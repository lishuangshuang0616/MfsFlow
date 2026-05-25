import os
import subprocess
import sys
import time
import hashlib
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime

_SRC_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "src")
if os.path.isdir(_SRC_DIR) and _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

from path_layout import logs_dir, tmp_merge_dir


def format_duration(seconds):
    seconds = float(seconds)
    if seconds < 60:
        return f"{seconds:.2f}s"
    minutes, sec = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)}m{sec:05.2f}s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)}h{int(minutes):02d}m{sec:05.2f}s"


class Tee:
    def __init__(self, *streams):
        self._streams = streams

    def write(self, s):
        for stream in self._streams:
            stream.write(s)

    def flush(self):
        for stream in self._streams:
            stream.flush()


class PipelineTimer:
    def __init__(self, timing_path, project):
        self.timing_path = timing_path
        self.project = project
        self._ensure_header()

    def _ensure_header(self):
        os.makedirs(os.path.dirname(self.timing_path), exist_ok=True)
        if not os.path.exists(self.timing_path) or os.path.getsize(self.timing_path) == 0:
            with open(self.timing_path, "w") as handle:
                handle.write("timestamp\tproject\tstage\tstatus\tduration_sec\tduration_human\tdetails\n")

    def record(self, stage, status, duration, details=""):
        safe_details = str(details or "").replace("\t", " ").replace("\n", " ")
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(self.timing_path, "a") as handle:
            handle.write(
                f"{ts}\t{self.project}\t{stage}\t{status}\t{duration:.3f}\t"
                f"{format_duration(duration)}\t{safe_details}\n"
            )

    @contextmanager
    def section(self, stage, details=""):
        start = time.perf_counter()
        print(f"[TIMING] START {stage}")
        try:
            yield
        except Exception:
            duration = time.perf_counter() - start
            self.record(stage, "failed", duration, details)
            print(f"[TIMING] FAIL  {stage}: {format_duration(duration)}")
            raise
        else:
            duration = time.perf_counter() - start
            self.record(stage, "ok", duration, details)
            print(f"[TIMING] DONE  {stage}: {format_duration(duration)}")


@dataclass
class PipelineTools:
    samtools: str
    pigz: str
    seqkit: str


@dataclass
class PipelineRuntime:
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
        toolkit_dir = config.get("toolkit_directory", ".")
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
        tmp_root = (config.get("performance_opts", {}) or {}).get("tmp_root")
        if not tmp_root:
            return tmp_merge_dir(out_dir)
        project = str(config.get("project", "sample"))
        safe_project = "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in project)
        out_hash = hashlib.sha1(os.path.abspath(out_dir).encode("utf-8")).hexdigest()[:10]
        return os.path.join(os.path.abspath(tmp_root), f"mfsflow_{safe_project}_{out_hash}", "tmp_merge")

    def install_src_path(self):
        toolkit_src_dir = os.path.join(self.toolkit_dir, "src")
        if os.path.isdir(toolkit_src_dir) and toolkit_src_dir not in sys.path:
            sys.path.insert(0, toolkit_src_dir)

    def resolve_script(self, script_name):
        direct = os.path.join(self.toolkit_dir, script_name)
        if os.path.exists(direct):
            return direct
        in_src = os.path.join(self.toolkit_dir, "src", script_name)
        if os.path.exists(in_src):
            return in_src
        raise FileNotFoundError(f"Script not found: {script_name}. Tried: {direct}, {in_src}")


def remove_path(path):
    if not path:
        return
    if os.path.exists(path):
        os.remove(path)
    bai = path + ".bai"
    if os.path.exists(bai):
        os.remove(bai)


def run_stage_cmd(cmd, stage_name, run_log, exec_env, timer, log_path, shell=False):
    if isinstance(cmd, list) and not shell:
        cmd_str = " ".join(map(str, cmd))
    else:
        cmd_str = str(cmd)
    start = time.perf_counter()
    print(f"[TIMING] START {stage_name}")
    res = subprocess.run(cmd, stdout=run_log, stderr=subprocess.STDOUT, shell=shell, env=exec_env)
    duration = time.perf_counter() - start
    if res.returncode != 0:
        timer.record(stage_name, "failed", duration, cmd_str)
        print(f"[TIMING] FAIL  {stage_name}: {format_duration(duration)}")
        run_log.flush()
        try:
            with open(log_path, "r") as lr:
                print(f"\n[ERROR] {stage_name} failed (rc={res.returncode}). Last 30 lines of log ({log_path}):\n", file=sys.stderr)
                print("".join(lr.readlines()[-30:]), file=sys.stderr)
        except Exception:
            pass
        raise RuntimeError(f"{stage_name} failed with exit code {res.returncode}.")
    timer.record(stage_name, "ok", duration, cmd_str)
    print(f"[TIMING] DONE  {stage_name}: {format_duration(duration)}")
