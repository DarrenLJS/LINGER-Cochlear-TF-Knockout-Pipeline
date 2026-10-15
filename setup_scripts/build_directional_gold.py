#!/usr/bin/env python3
"""
build_directional_gold.py — build the FROZEN directional gold-standard table (Track 4 / Stage 4, ADDED 2026-10-15).

Reads the raw tables fetched by setup_scripts/download_gold_standard.sh (provenance files only) plus the
hand-curated cochlear table, and writes resources/directional_gold.tsv, which the pipeline reads and never
regenerates at run time. Run it ONCE on the login node and commit the output:

    python3 setup_scripts/build_directional_gold.py \
        --gold-dir /exports/eddie/scratch/s2906787/cochlear_datasets/gold_standard \
        --tf-index  $SCRATCH/module7_tf_activity/tf_activity_summary.tsv \
        --curated   resources/cochlear_curated_gold.tsv \
        --out       resources/directional_gold.tsv

SEMANTICS. `expected_ko_sign` is the expected sign of CAPS for a KNOCKOUT of the TF. CAPS > 0 means the knockout moves
the network toward the aging signature, so:
  Tier A  cochlear curated  (PROXY: hearing-loss / hair-cell-maintenance phenotypes are not "aging-ward"; every row
          needs a PMID or DOI + an evidence tier, rows without either are DROPPED and reported)
  Tier B  GenAge model organisms, Mus musculus only
            longevity influence = Pro-Longevity   -> the gene promotes longevity, its loss shortens life -> KO CAPS > 0 (+1)
            longevity influence = Anti-Longevity  -> its loss extends life                               -> KO CAPS < 0 (-1)
            anything else (Necessary for fitness, Unclear, Unannotated) -> dropped
  Tier C  CellAge (human cell-senescence screens; symbols matched case-insensitively to the mouse TF index)
            Senescence Effect = Inhibits  -> the gene suppresses senescence, KO induces it -> +1
            Senescence Effect = Induces   -> KO suppresses senescence                       -> -1
            Unclear -> dropped
  Tier D  GenAge human: membership only (expected_ko_sign = 0). Used for enrichment ("is this TF aging-annotated at
          all"), never for sign accuracy.
A TF whose rows within one tier disagree on direction is EXCLUDED from that tier (listed in the provenance table).
Only TFs present in the --tf-index (the 586-TF network index) are kept. A tier with fewer than --min-tier-tfs (10)
TFs is flagged `underpowered` in the provenance table; the metrics code reports it rather than hiding it.

Raw files are located by recursive search under --gold-dir (case-insensitive: *genage_human*.csv,
*genage_models*.csv, *cellage*.tsv|csv) and their md5 recorded. Column names are matched case-insensitively with
whitespace stripped, and a missing required column stops the build.
"""
import argparse
import glob
import hashlib
import os
import sys

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--gold-dir", required=True)
p.add_argument("--tf-index", required=True, help="tf_activity_summary.tsv (columns = TFs, plus `role`) or a 1-column TF list")
p.add_argument("--curated", default=None, help="resources/cochlear_curated_gold.tsv (Tier A)")
p.add_argument("--out", required=True)
p.add_argument("--min-tier-tfs", type=int, default=10)
a = p.parse_args()


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for ch in iter(lambda: fh.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()


def find(patterns):
    hits = []
    for pat in patterns:
        hits += [f for f in glob.glob(os.path.join(a.gold_dir, "**", pat), recursive=True) if os.path.isfile(f)]
    hits = sorted(set(hits))
    if not hits:
        sys.exit(f"no file matching {patterns} under {a.gold_dir} — run download_gold_standard.sh first")
    if len(hits) > 1:
        print(f"[WARN] {len(hits)} files match {patterns}; using {hits[0]}: {hits}")
    return hits[0]


def read(path):
    df = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    df.columns = [c.strip().lower() for c in df.columns]
    return df


def need(df, cols, path):
    miss = [c for c in cols if c not in df.columns]
    if miss:
        sys.exit(f"{path}: missing column(s) {miss}; found {list(df.columns)}")


# ---- TF index -----------------------------------------------------------------------------------------------------
ti = pd.read_csv(a.tf_index, sep="\t", nrows=0) if a.tf_index.endswith((".tsv", ".txt")) else None
if ti is not None and "role" in ti.columns:
    tfs = [c for c in ti.columns if c != "role" and not c.startswith("Unnamed")]
else:
    tfs = [x.strip() for x in open(a.tf_index) if x.strip()]
tf_by_upper = {}
for t in tfs:
    tf_by_upper.setdefault(t.upper(), t)
print(f"TF index: {len(tfs)} TFs")

rows, prov, excluded = [], [], []
inputs = {}


def finish_tier(tier, label, recs, n_raw, files):
    df = pd.DataFrame(recs, columns=["tf", "expected_ko_sign", "evidence", "ref"]) if recs else \
        pd.DataFrame(columns=["tf", "expected_ko_sign", "evidence", "ref"])
    n_conf = 0
    if len(df) and tier != "D":
        g = df.groupby("tf")["expected_ko_sign"].nunique()
        conf = g[g > 1].index.tolist()
        n_conf = len(conf)
        for t in conf:
            excluded.append({"tier": tier, "tf": t, "reason": "conflicting directions within tier"})
        df = df[~df["tf"].isin(conf)]
    out = (df.groupby("tf", as_index=False)
             .agg(expected_ko_sign=("expected_ko_sign", "first"),
                  evidence=("evidence", lambda s: ";".join(sorted(set(map(str, s))))),
                  ref=("ref", lambda s: ";".join(sorted(set(x for x in map(str, s) if x))))))
    out.insert(1, "tier", tier)
    out.insert(2, "source", label)
    n_tf = len(out)
    prov.append({"tier": tier, "source": label, "files": ";".join(files), "n_raw_rows": n_raw, "n_rows_in_tf_index": len(recs),
                 "n_tfs_conflict_excluded": n_conf, "n_tfs": n_tf,
                 "n_pos": int((out["expected_ko_sign"] == 1).sum()), "n_neg": int((out["expected_ko_sign"] == -1).sum()),
                 "underpowered": bool(n_tf < a.min_tier_tfs)})
    rows.append(out)
    print(f"Tier {tier} {label}: {n_raw} raw rows -> {len(recs)} in TF index -> {n_tf} TFs "
          f"({int((out['expected_ko_sign'] == 1).sum())} +, {int((out['expected_ko_sign'] == -1).sum())} -)"
          f"{'  UNDERPOWERED' if n_tf < a.min_tier_tfs else ''}")


# ---- Tier A: cochlear curated ---------------------------------------------------------------------------------------
recs, n_raw = [], 0
if a.curated and os.path.exists(a.curated):
    c = pd.read_csv(a.curated, sep="\t", dtype=str, keep_default_na=False, comment="#")
    c.columns = [x.strip().lower() for x in c.columns]
    need(c, ["tf", "expected_ko_sign", "evidence_tier", "pmid", "doi"], a.curated)
    n_raw = len(c)
    inputs[a.curated] = md5(a.curated)
    for _, r in c.iterrows():
        ref = r["pmid"].strip() or r["doi"].strip()
        if not ref or not r["evidence_tier"].strip():
            excluded.append({"tier": "A", "tf": r["tf"], "reason": "no PMID/DOI or no evidence tier (unsourced rows are omitted)"})
            continue
        if r["expected_ko_sign"].strip() not in ("1", "-1", "+1"):
            excluded.append({"tier": "A", "tf": r["tf"], "reason": f"expected_ko_sign {r['expected_ko_sign']!r} not +1/-1"})
            continue
        t = tf_by_upper.get(r["tf"].strip().upper())
        if t is None:
            excluded.append({"tier": "A", "tf": r["tf"], "reason": "not in the TF index"})
            continue
        recs.append((t, int(r["expected_ko_sign"]), f"{r['evidence_tier'].strip()}", ref))
finish_tier("A", "cochlear_curated(PROXY)", recs, n_raw, [a.curated] if a.curated else [])

# ---- Tier B: GenAge models (mouse) -----------------------------------------------------------------------------------
pm = find(["*genage_models*.csv", "*models_genes*.csv", "*genage_models*.tsv"])
inputs[pm] = md5(pm)
m = read(pm)
need(m, ["symbol", "organism", "longevity influence"], pm)
mm = m[m["organism"].str.contains("Mus musculus", case=False)]
recs = []
for _, r in mm.iterrows():
    infl = r["longevity influence"].strip().lower()
    sign = 1 if infl == "pro-longevity" else (-1 if infl == "anti-longevity" else 0)
    t = tf_by_upper.get(r["symbol"].strip().upper())
    if sign != 0 and t is not None:
        recs.append((t, sign, "GenAge-models:" + r["longevity influence"].strip(), "genage-id:" + r.get("genage id", "")))
finish_tier("B", "GenAge_models_mouse", recs, len(mm), [pm])

# ---- Tier C: CellAge ---------------------------------------------------------------------------------------------------
pc = find(["*cellage*.tsv", "*cellage*.csv", "*cellAge*.tsv", "*cellAge*.csv"])
inputs[pc] = md5(pc)
cdf = read(pc)
need(cdf, ["gene symbol", "senescence effect"], pc)
recs = []
for _, r in cdf.iterrows():
    eff = r["senescence effect"].strip().lower()
    sign = 1 if eff.startswith("inhibit") else (-1 if eff.startswith("induc") else 0)
    t = tf_by_upper.get(r["gene symbol"].strip().upper())
    if sign != 0 and t is not None:
        recs.append((t, sign, "CellAge:" + r["senescence effect"].strip(), r.get("reference", "")))
finish_tier("C", "CellAge_human_symbols_case_insensitive", recs, len(cdf), [pc])

# ---- Tier D: GenAge human, membership only ---------------------------------------------------------------------------
ph = find(["*genage_human*.csv", "*human_genes*.csv", "*genage_human*.tsv"])
inputs[ph] = md5(ph)
h = read(ph)
need(h, ["symbol"], ph)
recs = []
for _, r in h.iterrows():
    t = tf_by_upper.get(r["symbol"].strip().upper())
    if t is not None:
        recs.append((t, 0, "GenAge-human:member", "genage-id:" + r.get("genage id", "")))
finish_tier("D", "GenAge_human_membership_only", recs, len(h), [ph])

gold = pd.concat(rows, ignore_index=True)
os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
gold = gold.sort_values(["tier", "tf"]).reset_index(drop=True)
gold.to_csv(a.out, sep="\t", index=False)
base = os.path.splitext(a.out)[0]
pd.DataFrame(prov).to_csv(base + "_provenance.tsv", sep="\t", index=False)
pd.DataFrame(excluded, columns=["tier", "tf", "reason"]).to_csv(base + "_excluded.tsv", sep="\t", index=False)
pd.DataFrame([{"file": k, "md5": v} for k, v in inputs.items()]).to_csv(base + "_inputs_md5.tsv", sep="\t", index=False)
print(f"\nwrote {a.out} ({len(gold)} rows), {base}_provenance.tsv, {base}_excluded.tsv, {base}_inputs_md5.tsv")
print(pd.DataFrame(prov)[["tier", "n_tfs", "n_pos", "n_neg", "underpowered"]].to_string(index=False))
