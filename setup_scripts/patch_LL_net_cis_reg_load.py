"""
patch_LL_net_cis_reg_load.py — fixes a second, distinct real bug in
installed LingerGRN's cell_type_specific_cis_reg(), the else branch
(single-celltype + method=='scNN' — exactly what grn_celltype_specific.py
uses, confirmed via the real traceback on Eddie).

Bug: unlike its sibling branch (celltype=='all' & method=='scNN'), which
correctly calls `distance,cisGRN,REs,TGs=load_RE_TG_scNN(outdir)` before
using those variables, this else branch never calls load_RE_TG_scNN() at
all — it uses distance/cisGRN/REs/TGs immediately, causing
UnboundLocalError. Not an ordering bug this time (unlike the earlier
cell_type_specific_TF_RE_binding fix) — the line is simply absent from
this branch. Fixed by adding it, matching the sibling branch exactly.

Run this ONCE on Eddie (separate from, and safe to run alongside,
patch_LL_net_RE_ordering.py — different function, different block):
    /exports/csce/eddie/biology/groups/bioinfmsc/anaconda/envs/s2906787/LINGER/bin/python \
        patch_LL_net_cis_reg_load.py

Idempotent: safe to re-run, does nothing if already patched. Makes a
.bak backup before writing (separate backup file from the earlier patch).
"""
import LingerGRN.LL_net as _mod
import os

target = _mod.__file__
print(f"Patching: {target}")

with open(target) as f:
    content = f.read()

BROKEN = """    else: 
        label0=celltype
        label0=str(label0)
        temp=adata_ATAC.X[np.array(label)==label0,:].mean(axis=0).T
        RE=pd.DataFrame(temp,index=adata_ATAC.var['gene_ids'].values,columns=['values'])
        temp=adata_RNA.X[np.array(label)==label0,:].mean(axis=0).T
        TG=pd.DataFrame(temp,index=adata_RNA.var['gene_ids'].values,columns=['values'])
        del temp
        result=cell_type_specific_cis_reg_scNN(distance,cisGRN,RE,TG,REs,TGs)
        result.to_csv(outdir+'cell_type_specific_cis_regulatory_'+label0+'.txt',sep='\\t',header=None,index=None)"""

FIXED = """    else: 
        # --- FIXED 2026-09-07: this branch never called load_RE_TG_scNN(),
        # unlike its sibling (celltype=='all' & method=='scNN') branch just
        # above, which does — causing UnboundLocalError on distance/cisGRN/
        # REs/TGs for every single-celltype scNN call (confirmed: this is
        # exactly grn_celltype_specific.py's call pattern). Added to match
        # the working sibling branch. ---
        distance,cisGRN,REs,TGs=load_RE_TG_scNN(outdir)
        label0=celltype
        label0=str(label0)
        temp=adata_ATAC.X[np.array(label)==label0,:].mean(axis=0).T
        RE=pd.DataFrame(temp,index=adata_ATAC.var['gene_ids'].values,columns=['values'])
        temp=adata_RNA.X[np.array(label)==label0,:].mean(axis=0).T
        TG=pd.DataFrame(temp,index=adata_RNA.var['gene_ids'].values,columns=['values'])
        del temp
        result=cell_type_specific_cis_reg_scNN(distance,cisGRN,RE,TG,REs,TGs)
        result.to_csv(outdir+'cell_type_specific_cis_regulatory_'+label0+'.txt',sep='\\t',header=None,index=None)"""

if "# --- FIXED 2026-09-07: this branch never called load_RE_TG_scNN" in content:
    print("Already patched — nothing to do.")
elif BROKEN not in content:
    print("ERROR: expected broken block not found verbatim in installed file.")
    print("Do NOT guess further — paste back the surrounding lines so the")
    print("patch can be adjusted:")
    print(f"    grep -n 'cell_type_specific_cis_reg_scNN(distance' {target}")
    raise SystemExit(1)
else:
    backup = target + ".bak_2026-09-07_cis_reg"
    if not os.path.exists(backup):
        with open(backup, "w") as f:
            f.write(content)
        print(f"Backup written: {backup}")
    content = content.replace(BROKEN, FIXED)
    with open(target, "w") as f:
        f.write(content)
    print("Patched successfully.")
