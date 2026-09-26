#!/usr/bin/env python3
"""
geo_lookup_full.py — Stage A metadata triage for the condition-mixing audit.

Extends the original geo_lookup.py (which hand-listed 22 extra_bulk_rna_inputs
accessions) to EVERY unique GSE accession referenced anywhere in
config_eddie.yaml (69 as of 2026-09-26), auto-discovered from the config
itself rather than hand-typed, and tags each one with which config
section(s)/sample_id(s) reference it. That mapping is what turns a bare GEO
title into something actionable: "does this dataset's own description
suggest multiple conditions, AND does the pipeline currently keep those
conditions separate or blend them into one entry."

This does NOT touch any downloaded data and needs no conda env beyond
`requests`/`pyyaml` — it queries NCBI's E-utilities (esearch+esummary against
the "gds" database) for each accession's real title/summary/n_samples.

Needs outbound internet access to eutils.ncbi.nlm.nih.gov (HTTPS, port
443). Run anywhere with internet and the repo checked out — a laptop or
Eddie's login node both work, no downloaded datasets required.

Usage:
    pip install requests pyyaml   # if not already available
    python geo_lookup_full.py --config /path/to/config_eddie.yaml
    # writes geo_titles_full.json next to wherever you run this from
"""
import argparse
import json
import re
import time

import requests
import yaml

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

# Cheap first-pass keyword flag on the GEO title+summary text. This is NOT
# the final verdict — it's just "worth a closer structural look" (Stage B
# is the actual verdict, since a risky-sounding title can already be
# handled correctly in config, e.g. GSE120462 was split into two entries
# specifically to avoid this).
CONDITION_KEYWORDS = re.compile(
    r"\b("
    r"wild.?type|wildtype|\bwt\b|knockout|\bko\b|knock.?out|"
    r"heterozygous|homozygous|\bhet\b|\bhomo\b|mutant|\bmut\b|"
    r"control|\bctrl\b|\bctl\b|treated|untreated|deficient|"
    r"conditional|floxed|cre|null|-/-|"
    r"aging|aged|young|old\b|postnatal day|\bp\d+\b|\bday \d+\b|month"
    r")\b",
    re.IGNORECASE,
)


def find_all_accessions(config):
    """Walk the whole parsed YAML looking for GSE-accession-shaped strings,
    so this stays correct automatically if config_eddie.yaml grows more
    entries later rather than drifting from a hand-maintained list."""
    found = set()

    def walk(obj):
        if isinstance(obj, str):
            for m in re.finditer(r"GSE\d+", obj):
                found.add(m.group(0))
        elif isinstance(obj, dict):
            for v in obj.values():
                walk(v)
        elif isinstance(obj, list):
            for v in obj:
                walk(v)

    walk(config)
    return found


def build_accession_index(config):
    """accession -> list of {section, sample_id} referencing it, so the
    output JSON tells you not just "what is GSE163798" but "which config
    entries load it, under which category"."""
    index = {}

    def add(acc, section, sample_id):
        index.setdefault(acc, []).append({"section": section, "sample_id": sample_id})

    def scan_list(section_name, entries):
        for e in entries or []:
            if not isinstance(e, dict):
                continue
            sample_id = e.get("sample_id") or e.get("ref_id")
            for v in e.values():
                if isinstance(v, str):
                    for m in re.finditer(r"GSE\d+", v):
                        add(m.group(0), section_name, sample_id)
            # multi_file_merge entries nest a "files" list with their own patterns
            for f in e.get("files", []) or []:
                for v in f.values() if isinstance(f, dict) else []:
                    if isinstance(v, str):
                        for m in re.finditer(r"GSE\d+", v):
                            add(m.group(0), section_name, sample_id)

    scan_list("samples", config.get("samples"))
    scan_list("extra_bulk_rna_inputs", config.get("extra_bulk_rna_inputs"))
    scan_list("model_construction_refs", config.get("model_construction_refs"))
    scan_list("tuning_inputs", config.get("tuning_inputs"))
    scan_list("held_out_inputs", config.get("held_out_inputs"))
    noc = config.get("negative_control_input")
    if noc:
        scan_list("negative_control_input", [noc])

    return index


def lookup(acc):
    r = requests.get(ESEARCH, params={
        "db": "gds", "term": f"{acc}[Accession]", "retmode": "json",
    }, timeout=15)
    r.raise_for_status()
    ids = r.json()["esearchresult"]["idlist"]
    if not ids:
        return {"accession": acc, "title": None, "summary": None, "error": "not found"}

    uid = ids[0]
    r = requests.get(ESUMMARY, params={
        "db": "gds", "id": uid, "retmode": "json",
    }, timeout=15)
    r.raise_for_status()
    doc = r.json()["result"][uid]
    return {
        "accession": acc,
        "title": doc.get("title"),
        "summary": doc.get("summary"),
        "taxon": doc.get("taxon"),
        "n_samples": doc.get("n_samples"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config_eddie.yaml",
                     help="Path to config_eddie.yaml (default: config/config_eddie.yaml, "
                          "i.e. run this from the repo root)")
    ap.add_argument("--out", default="geo_titles_full.json")
    args = ap.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    accessions = sorted(find_all_accessions(config))
    index = build_accession_index(config)
    print(f"Found {len(accessions)} unique GSE accessions in {args.config}")

    results = []
    for acc in accessions:
        try:
            info = lookup(acc)
        except Exception as e:
            info = {"accession": acc, "title": None, "summary": None, "error": str(e)}
        info["used_by"] = index.get(acc, [])
        text = f"{info.get('title') or ''} {info.get('summary') or ''}"
        matches = sorted(set(m.group(0).lower() for m in CONDITION_KEYWORDS.finditer(text)))
        info["condition_keyword_hits"] = matches
        info["flag_multi_condition_suspected"] = len(matches) > 0
        results.append(info)
        flag = "FLAG" if info["flag_multi_condition_suspected"] else "    "
        print(f"[{flag}] {acc} (used by {len(info['used_by'])} entr{'y' if len(info['used_by'])==1 else 'ies'}): "
              f"{info.get('title') or info.get('error')}")
        time.sleep(0.4)  # stay under NCBI's unauthenticated rate limit

    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)

    n_flagged = sum(r.get("flag_multi_condition_suspected") for r in results)
    print(f"\n{n_flagged}/{len(results)} accessions flagged on keyword heuristic alone "
          f"(this is a first pass, not a verdict — some are already split correctly in "
          f"config, e.g. GSE120462; Stage B settles it against the real data).")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
