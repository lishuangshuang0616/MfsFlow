import os
import sys


MIN_PYTHON = (3, 8)


def require_supported_python():
    if sys.version_info < MIN_PYTHON:
        required = ".".join(map(str, MIN_PYTHON))
        current = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
        raise RuntimeError(f"Python {required}+ is required. Current interpreter: Python {current}")


def validate_input_files(config):
    files_to_check = [config["reference"]["STAR_index"], config["reference"]["GTF_file"]]
    for path in files_to_check:
        if not path or not os.path.exists(path):
            raise FileNotFoundError(f"Reference file not found: {path}")

    for fastq_group in ("file1", "file2"):
        fastq_names = config["sequence_files"][fastq_group]["name"]
        if not fastq_names:
            continue
        for fastq_path in fastq_names.split(","):
            fastq_path = fastq_path.strip()
            if fastq_path and not os.path.exists(fastq_path):
                raise FileNotFoundError(f"Fastq file not found: {fastq_path}")
