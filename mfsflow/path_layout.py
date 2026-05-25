import os


def config_dir(out_dir):
    return os.path.join(out_dir, "config")


def logs_dir(out_dir):
    return os.path.join(out_dir, "logs")


def barcode_dir(out_dir):
    return os.path.join(out_dir, "barcodes")


def expression_dir(out_dir):
    return os.path.join(out_dir, "expression")


def stats_dir(out_dir):
    return os.path.join(out_dir, "stats")


def intermediate_dir(out_dir):
    return os.path.join(out_dir, "intermediate")


def tmp_merge_dir(out_dir):
    return os.path.join(intermediate_dir(out_dir), "tmp_merge")


def outputs_dir(out_dir):
    return os.path.join(os.path.dirname(out_dir), "outs")


def ensure_layout(out_dir):
    for path in (
        out_dir,
        config_dir(out_dir),
        logs_dir(out_dir),
        barcode_dir(out_dir),
        expression_dir(out_dir),
        stats_dir(out_dir),
        intermediate_dir(out_dir),
        tmp_merge_dir(out_dir),
        outputs_dir(out_dir),
    ):
        os.makedirs(path, exist_ok=True)
