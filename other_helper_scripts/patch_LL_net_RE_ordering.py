"""
patch_LL_net_RE_ordering.py — fixes a real ordering bug in installed
LingerGRN's cell_type_specific_TF_RE_binding(), else branch (single-celltype
call path — exactly what grn_celltype_specific.py uses).

Bug, confirmed against the real installed file's own traceback:
    REoverlap=list(set(TFbinding.index)&set(RE.index))
references `RE` before it's ever assigned in this branch — the loop branch
(unused by our pipeline) computes RE correctly beforehand; the else branch
has the RE-computation block placed several lines too late. Fixed by moving
that block above the REoverlap line, matching the loop branch's working
order.

Run this ONCE on Eddie:
    /exports/csce/eddie/biology/groups/bioinfmsc/anaconda/envs/s2906787/LINGER/bin/python \
        patch_LL_net_RE_ordering.py

Idempotent: safe to re-run, does nothing if already patched. Makes a
.bak backup before writing.
"""
import LingerGRN.LL_net as _mod
import os

target = _mod.__file__
print(f"Patching: {target}")

with open(target) as f:
    content = f.read()

BROKEN = """        TFbinding = TFbinding[TFoverlap]
        REoverlap=list(set(TFbinding.index)&set(RE.index))
        TFbinding=TFbinding.loc[REoverlap]
        TFbinding1=np.zeros((mat.shape[0],len(TFoverlap)))
        REidx=pd.DataFrame(range(mat.shape[0]),index=mat.index)
        TFbinding1[REidx.loc[TFbinding.index][0].values,:]=TFbinding.values
        TFbinding1 = pd.DataFrame(TFbinding1,index=mat.index,columns=TFoverlap)
        print('Generate cell type specitic TF binding potential for cell type '+ str(label0)+'...')
        temp=adata_ATAC.X[np.array(label)==label0,:].mean(axis=0).T
        RE=pd.DataFrame(temp,index=adata_ATAC.var['gene_ids'].values,columns=['values'])
        temp=adata_RNA.X[np.array(label)==label0,:].mean(axis=0).T
        TG=pd.DataFrame(temp,index=adata_RNA.var['gene_ids'].values,columns=['values'])
        RE=RE.loc[REs]
        result=cell_type_specific_TF_RE_binding_score_scNN(mat,TFbinding,RE,TG,TFoverlap)"""

FIXED = """        TFbinding = TFbinding[TFoverlap]
        print('Generate cell type specitic TF binding potential for cell type '+ str(label0)+'...')
        temp=adata_ATAC.X[np.array(label)==label0,:].mean(axis=0).T
        RE=pd.DataFrame(temp,index=adata_ATAC.var['gene_ids'].values,columns=['values'])
        temp=adata_RNA.X[np.array(label)==label0,:].mean(axis=0).T
        TG=pd.DataFrame(temp,index=adata_RNA.var['gene_ids'].values,columns=['values'])
        RE=RE.loc[REs]
        # --- FIXED 2026-09-07: RE-computation block (above) moved before
        # REoverlap, which references RE.index. Original source had this
        # backwards, causing UnboundLocalError on every single-celltype
        # call (confirmed: the else branch here is exactly what
        # grn_celltype_specific.py uses, one celltype at a time). ---
        REoverlap=list(set(TFbinding.index)&set(RE.index))
        TFbinding=TFbinding.loc[REoverlap]
        TFbinding1=np.zeros((mat.shape[0],len(TFoverlap)))
        REidx=pd.DataFrame(range(mat.shape[0]),index=mat.index)
        TFbinding1[REidx.loc[TFbinding.index][0].values,:]=TFbinding.values
        TFbinding1 = pd.DataFrame(TFbinding1,index=mat.index,columns=TFoverlap)
        result=cell_type_specific_TF_RE_binding_score_scNN(mat,TFbinding,RE,TG,TFoverlap)"""

if "# --- FIXED 2026-09-07" in content:
    print("Already patched — nothing to do.")
elif BROKEN not in content:
    print("ERROR: expected broken block not found verbatim in installed file.")
    print("The installed version's exact whitespace/text may differ from what")
    print("this patch expects. Do NOT guess further — paste back the ~15 lines")
    print("surrounding line 565 of the actual file so the patch can be adjusted:")
    print(f"    sed -n '555,570p' {target}")
    raise SystemExit(1)
else:
    backup = target + ".bak_2026-09-07"
    if not os.path.exists(backup):
        with open(backup, "w") as f:
            f.write(content)
        print(f"Backup written: {backup}")
    content = content.replace(BROKEN, FIXED)
    with open(target, "w") as f:
        f.write(content)
    print("Patched successfully.")
