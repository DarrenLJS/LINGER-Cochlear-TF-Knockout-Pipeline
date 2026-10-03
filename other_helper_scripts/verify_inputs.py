#!/usr/bin/env python3
"""
verify_inputs.py — fail-loud check that every raw input the workflow reads
exists on disk, before Snakemake is started.

What it checks
--------------
It loads config/config_eddie.yaml and walks the seven sections the workflow
actually reads:

    samples, extra_bulk_rna_inputs, extra_chromatin_prior_inputs,
    chromatin_prior_extra, model_construction_refs, held_out_inputs,
    negative_control_input

(tuning_inputs and generalization_test_panel are NOT read by any rule, so they
are deliberately not checked.)

For every entry it resolves {multiome_root} / {chrom_root} / {scratch} and checks:

  path            exists; if it is a directory it must hold at least one data
                  file (a folder with only filelist.txt, *_RAW.tar and
                  .extracted markers means the RAW.tar was never extracted)
  h5_path, atac_fragments, atac_peaks_bed
                  the file exists
  combined_prefix, rna_prefix, atac_prefix
                  at least one file starts with that prefix
  glob_pattern    matches at least one file under the entry's path
  *_file          (inside an entry that has a path) exists under that path

Nothing is downloaded, written or modified.

Exit status: 0 = every input found, 1 = something is missing or empty.

Usage (from the repo root; needs only PyYAML):
    python other_helper_scripts/verify_inputs.py --config config/config_eddie.yaml
"""
import argparse
import glob
import os
import sys

import yaml

SECTIONS = [
    "samples",
    "extra_bulk_rna_inputs",
    "extra_chromatin_prior_inputs",
    "chromatin_prior_extra",
    "model_construction_refs",
    "held_out_inputs",
    "negative_control_input",
]
FILE_KEYS = {"h5_path", "atac_fragments", "atac_peaks_bed"}
PREFIX_KEYS = {"combined_prefix", "rna_prefix", "atac_prefix"}
ID_KEYS = ("sample_id", "prior_id", "ref_id")
# files that are bookkeeping, not data
NON_DATA = ("filelist.txt",)
NON_DATA_SUFFIX = ("_RAW.tar", ".extracted")


def build_resolver(cfg):
    subs = {
        "{multiome_root}": cfg["multiome_root"],
        "{chrom_root}": cfg["chrom_root"],
        "{scratch}": cfg["scratch"],
    }

    def resolve(s):
        for k, v in subs.items():
            s = s.replace(k, v)
        return os.path.normpath(s)

    return resolve


def data_files(d):
    out = []
    for name in os.listdir(d):
        if name in NON_DATA or name.endswith(NON_DATA_SUFFIX):
            continue
        out.append(name)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--config", default="config/config_eddie.yaml")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    for k in ("multiome_root", "chrom_root", "scratch"):
        if k not in cfg:
            sys.exit(f"ERROR: '{k}' missing from {args.config}")
    resolve = build_resolver(cfg)

    checks = {}  # (kind, resolved target, extra) -> set of labels

    def add(kind, target, label, extra=""):
        checks.setdefault((kind, target, extra), set()).add(label)

    def walk(node, base, label):
        if isinstance(node, list):
            for x in node:
                walk(x, base, label)
            return
        if not isinstance(node, dict):
            return
        for idk in ID_KEYS:
            if isinstance(node.get(idk), str):
                label = node[idk]
                break
        here = base
        if isinstance(node.get("path"), str):
            p = resolve(node["path"])
            if "{" in p:
                print(f"NOTE  unresolved placeholder, skipped: {node['path']}")
            else:
                here = p
                add("path", p, label)
        for k, v in node.items():
            if isinstance(v, str):
                if "{" in v and k in (FILE_KEYS | PREFIX_KEYS):
                    r = resolve(v)
                    if "{" in r:
                        print(f"NOTE  unresolved placeholder, skipped: {v}")
                        continue
                if k in FILE_KEYS:
                    add("file", resolve(v), label)
                elif k in PREFIX_KEYS:
                    add("prefix", resolve(v), label)
                elif k == "glob_pattern" and here:
                    add("glob", os.path.join(here, v), label)
                elif k.endswith("_file") and here:
                    f = v if os.path.isabs(v) else os.path.join(here, v)
                    add("file", f, label)
            elif isinstance(v, (dict, list)):
                walk(v, here, label)

    for sec in SECTIONS:
        if sec not in cfg:
            sys.exit(f"ERROR: section '{sec}' missing from {args.config}")
        walk(cfg[sec], None, sec)

    print(f"config          : {args.config}")
    print(f"multiome_root   : {cfg['multiome_root']}")
    print(f"chrom_root      : {cfg['chrom_root']}")
    print(f"checks          : {len(checks)} distinct targets\n")

    problems = []
    for (kind, target, _), labels in sorted(checks.items()):
        who = ", ".join(sorted(labels))
        if kind == "path":
            if os.path.isdir(target):
                if not data_files(target):
                    problems.append(("EMPTY  (only filelist/RAW.tar: not extracted?)", target, who))
            elif not os.path.exists(target):
                problems.append(("MISSING dir/file", target, who))
        elif kind == "file":
            if not os.path.isfile(target):
                problems.append(("MISSING file", target, who))
        elif kind == "prefix":
            if not glob.glob(glob.escape(target) + "*"):
                problems.append(("MISSING prefix (no file starts with it)", target, who))
        elif kind == "glob":
            if not glob.glob(target):
                problems.append(("NO MATCH for glob_pattern", target, who))

    if problems:
        print(f"FAILED: {len(problems)} input problem(s)\n")
        for what, target, who in problems:
            print(f"  {what}\n    {target}\n    used by: {who}")
        sys.exit(1)
    print("OK: every raw input the workflow reads is present.")


if __name__ == "__main__":
    main()
