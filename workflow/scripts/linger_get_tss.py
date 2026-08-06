"""
linger_get_tss.py — Module 4 step 2.

Runs LingerGRN.LINGER_tr.get_TSS() + RE_TG_dis(), the two calls the official
docs/scNN.md tutorial runs before training. Confirmed from source: BOTH
functions read/write relative `./data/...` paths regardless of the outdir
argument, so this script chdir's into snakemake.params.workdir (the same
directory linger_prep_pseudobulk.py wrote data/Peaks.txt etc into) before
calling them.

genome='mm10' is directly supported here (get_TSS reads
GRNdir + 'TSS_mm10.txt', shipped in LINGER's own reference bundle) — this is
the documented mouse-mode path, not something we're improvising.

CLI, not `script:` — 2026-07-25 fix, see linger_prep_pseudobulk.py docstring
for why (conda: directory-path misclassification under real SGE execution).
"""
import argparse
import os
import LingerGRN.LINGER_tr as LINGER_tr

p = argparse.ArgumentParser()
p.add_argument("--workdir", required=True)
p.add_argument("--grn-dir", required=True)
p.add_argument("--genome", required=True)
p.add_argument("--tss-distance-bp", type=int, required=True)
p.add_argument("--output-done", required=True)
args = p.parse_args()

WORKDIR = args.workdir
# LINGER_tr.get_TSS() builds its path via plain string concatenation
# (GRNdir + 'TSS_' + genome + '.txt'), NOT os.path.join, so GRNdir must
# already end in a separator or the filename gets glued onto the last
# path component (e.g. ".../provide_dataTSS_mm10.txt"). Normalize here
# regardless of what the caller passed in --grn-dir.
GRN_DIR = args.grn_dir.rstrip("/") + "/"
GENOME = args.genome
TSS_DIS = args.tss_distance_bp

assert os.path.isdir(GRN_DIR), (
    f"linger.grn_dir ({GRN_DIR}) does not exist. This must be LINGER's own "
    "downloadable reference bundle ('provide_data/' per docs/scNN.md) — "
    "see SETUP.md 'Module 4 prerequisite' section. This pipeline does not "
    "build this bundle itself."
)
tss_file = os.path.join(GRN_DIR, f"TSS_{GENOME}.txt")
assert os.path.isfile(tss_file), (
    f"{tss_file} not found in GRNdir — confirm the reference bundle actually "
    f"includes {GENOME} support (it does per LINGER's public docs, but a "
    "partial/incomplete download would fail exactly here)."
)

os.chdir(WORKDIR)
os.makedirs("data", exist_ok=True)

print(f"get_TSS(GRNdir={GRN_DIR}, genome={GENOME}, TSS_dis={TSS_DIS})")
LINGER_tr.get_TSS(GRN_DIR, GENOME, TSS_DIS)

print("RE_TG_dis(outdir=./)")
LINGER_tr.RE_TG_dis("./")  # writes data/RE_gene_distance.txt

with open(args.output_done, "w") as f:
    f.write("get_TSS + RE_TG_dis complete\n")
    f.write(f"workdir: {WORKDIR}\ngenome: {GENOME}\ntss_distance_bp: {TSS_DIS}\n")

print(f"\nDone. data/TSS_extend_1M.txt, data/RE_gene_distance.txt written under {WORKDIR}/data/")
