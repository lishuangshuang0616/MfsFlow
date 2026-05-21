import copy
import os

import yaml


def write_run_config(config):
    run_config = copy.deepcopy(config)
    toolkit_dir = config.get("toolkit_directory") or "."
    run_config["toolkit_directory"] = toolkit_dir
    software_dir = os.path.join(toolkit_dir, "software")

    def resolve_tool(name):
        candidate = os.path.join(software_dir, name)
        if os.path.exists(candidate) and os.access(candidate, os.X_OK):
            return os.path.abspath(candidate)
        return name

    run_config["samtools_exec"] = resolve_tool("samtools")
    run_config["pigz_exec"] = resolve_tool("pigz")
    run_config["seqkit_exec"] = resolve_tool("seqkit")
    run_config["STAR_exec"] = resolve_tool("STAR")
    run_config["featureCounts_exec"] = resolve_tool("featureCounts")

    class ForceStr:
        def __init__(self, value):
            self.value = value

    def force_str_representer(dumper, data):
        return dumper.represent_scalar("tag:yaml.org,2002:str", data.value, style='"')

    class RunConfigDumper(yaml.SafeDumper):
        pass

    def bool_representer(dumper, value):
        return dumper.represent_scalar("tag:yaml.org,2002:bool", "yes" if value else "no")

    def none_representer(dumper, _value):
        return dumper.represent_scalar("tag:yaml.org,2002:null", "~")

    yaml.add_representer(ForceStr, force_str_representer, Dumper=RunConfigDumper)
    yaml.add_representer(bool, bool_representer, Dumper=RunConfigDumper)
    yaml.add_representer(type(None), none_representer, Dumper=RunConfigDumper)

    ds = run_config["counting_opts"].get("downsampling", "0")
    if isinstance(ds, list):
        ds = ",".join(map(str, ds))
    else:
        ds = str(ds)
    run_config["counting_opts"]["downsampling"] = ForceStr(ds)

    final_yaml_path = os.path.join(config["out_dir"], "config", "run_config.yaml")
    os.makedirs(os.path.dirname(final_yaml_path), exist_ok=True)
    with open(final_yaml_path, "w") as f:
        yaml.dump(run_config, f, Dumper=RunConfigDumper, default_flow_style=False, sort_keys=False)
    return final_yaml_path
