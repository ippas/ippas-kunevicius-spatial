#!/usr/bin/env python3
"""
Functional (over-representation) analysis of differentially expressed genes
Maternal FMT model (ASD vs neurotypical donors) – Visium spatial transcriptomics, offspring brain.

What the script does
  1. downloads the MSigDB v2026.1 gene sets (msigdbdf R package data, release v26.1.0) if not present,
  2. runs over-representation analysis (ORA) for 8 databases, separately for all / up / down genes,
  3. writes all result tables to an Excel file and the overview + GO:CC figures (PNG/PDF),
  4. writes the software versions and run parameters to run_info.txt.

Usage
  pip install -r requirements.txt
  python run_functional_analysis.py                       # uses input_DE_genes.tsv
  python run_functional_analysis.py my_genes.tsv out_dir  # own input / output folder

Input: tab-separated file with columns  gene  log2FC   (other columns are ignored;
       a gene listed in several clusters is counted once; direction = sign of log2FC;
       a gene that is up in one cluster and down in another is placed in both lists).
"""
import sys, os, hashlib, platform, urllib.request, datetime
import numpy as np, pandas as pd, pyreadr, scipy
from scipy.stats import hypergeom
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# ------------------------------------------------------------------ parameters
INPUT   = sys.argv[1] if len(sys.argv) > 1 else "input_DE_genes.tsv"
OUTDIR  = sys.argv[2] if len(sys.argv) > 2 else "results"
MIN_SIZE    = 10      # gene sets with fewer genes are not tested
MAX_SIZE    = None    # no upper size limit
MIN_OVERLAP = 2       # terms with < 2 overlapping genes get p = 1
FDR_CUT     = 0.10    # significance threshold (Benjamini–Hochberg)
MSIGDB_URL  = "https://raw.githubusercontent.com/igordot/msigdbdf/v26.1.0/R/sysdata.rda"
MSIGDB_FILE = "msigdb_v2026.1_sysdata.rda"
MSIGDB_GIT_BLOB = "3ee38e4360cf88b23994ea2b7d3838e68f3e7341"   # git blob hash of the file used for the presentation

# database name -> (species of the gene sets, selection rule on the gene-set table)
DATABASES = [
    ("GO Biological Process",    "mm", lambda d: d.gs_subcollection == "GO:BP"),
    ("GO Molecular Function",    "mm", lambda d: d.gs_subcollection == "GO:MF"),
    ("GO Cellular Component",    "mm", lambda d: d.gs_subcollection == "GO:CC"),
    ("Reactome",                 "mm", lambda d: d.gs_subcollection == "CP:REACTOME"),
    ("KEGG",                     "hs", lambda d: d.gs_subcollection.isin(["CP:KEGG_LEGACY", "CP:KEGG_MEDICUS"])),
    ("WikiPathways",             "mm", lambda d: d.gs_subcollection == "CP:WIKIPATHWAYS"),
    ("MSigDB Hallmarks",         "mm", lambda d: d.gs_collection == "MH"),
    ("Human Phenotype Ontology", "hs", lambda d: d.gs_subcollection == "HPO"),
]

# ------------------------------------------------------------------ helpers
def git_blob_hash(path):
    data = open(path, "rb").read()
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()

def get_msigdb():
    if not os.path.exists(MSIGDB_FILE):
        print(f"Downloading MSigDB v2026.1 gene sets (~13 MB) from\n  {MSIGDB_URL}")
        urllib.request.urlretrieve(MSIGDB_URL, MSIGDB_FILE)
    h = git_blob_hash(MSIGDB_FILE)
    if h != MSIGDB_GIT_BLOB:
        print(f"WARNING: gene-set file differs from the one used for the presentation ({h})")
    return pyreadr.read_r(MSIGDB_FILE)

def bh(p):
    """Benjamini–Hochberg adjusted p-values."""
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out

def pretty(term):
    for pre in ("GOBP_", "GOMF_", "GOCC_", "REACTOME_", "WP_", "KEGG_MEDICUS_", "KEGG_", "HALLMARK_", "HP_"):
        if term.startswith(pre): term = term[len(pre):]; break
    return term.replace("_", " ").capitalize()

# ------------------------------------------------------------------ data
os.makedirs(OUTDIR, exist_ok=True)
de = pd.read_csv(INPUT, sep="\t")
lists = {"all":  set(de.gene),
         "up":   set(de.loc[de.log2FC > 0, "gene"]),
         "down": set(de.loc[de.log2FC < 0, "gene"])}
direction = de.groupby("gene").log2FC.apply(lambda x: "up" if (x > 0).all() else ("down" if (x < 0).all() else "mixed"))
print(f"Input: {len(de)} rows, {len(lists['all'])} unique genes ({len(lists['up'])} up, {len(lists['down'])} down)")
msig = get_msigdb()

# ------------------------------------------------------------------ ORA
summary, full = [], {}
for name, sp, rule in DATABASES:
    det, mem = msig[f"gene_set_details_{sp}"], msig[f"gene_set_members_{sp}"]
    d = det[rule(det)]
    m = mem[mem.gs_id.isin(d.gs_id)]
    universe = set(m.db_gene_symbol)                           # background = all genes of this database
    sets = m.groupby("gs_id").db_gene_symbol.apply(set)
    sets = sets[sets.apply(len) >= MIN_SIZE]
    if MAX_SIZE: sets = sets[sets.apply(len) <= MAX_SIZE]
    info = d.set_index("gs_id")
    rec = {"database": name, "gene_sets": "mouse" if sp == "mm" else "human (symbols matched in upper case)",
           "background_genes": len(universe), "terms_tested": len(sets)}
    tabs = []
    for lab, genes in lists.items():
        sym = {g.upper(): g for g in genes} if sp == "hs" else {g: g for g in genes}
        q = set(sym) & universe; N, n = len(universe), len(q)
        rows = []
        for gid, gs in sets.items():
            hit = sorted(gs & q); k = len(hit)
            p = hypergeom.sf(k - 1, N, len(gs), n) if k >= MIN_OVERLAP else 1.0
            orig = [sym[h] for h in hit]
            rows.append({"list": lab, "term": info.at[gid, "gs_name"], "term_pretty": pretty(info.at[gid, "gs_name"]),
                         "source_id": info.at[gid, "gs_exact_source"], "overlap": k, "term_size": len(gs),
                         "input_genes_in_database": n, "background": N,
                         "fold_enrichment": (k / n) / (len(gs) / N) if n and k else 0.0, "p_value": p,
                         "genes": ",".join(orig),
                         "genes_up": ",".join(g for g in orig if direction.get(g) == "up"),
                         "genes_down": ",".join(g for g in orig if direction.get(g) == "down")})
        t = pd.DataFrame(rows); t["FDR_BH"] = bh(t.p_value); t = t.sort_values("p_value")
        top = t.iloc[0]
        rec.update({f"{lab}_input_genes_in_database": n, f"{lab}_n_significant": int((t.FDR_BH < FDR_CUT).sum()),
                    f"{lab}_top_term": top.term_pretty, f"{lab}_top_overlap": int(top.overlap), f"{lab}_top_FDR": top.FDR_BH})
        tabs.append(t)
    full[name] = pd.concat(tabs)
    summary.append(rec)
    print(f"  {name:26s} terms tested {len(sets):5d} | significant (all genes): {rec['all_n_significant']}")
summary = pd.DataFrame(summary)

# ------------------------------------------------------------------ tables
xlsx = os.path.join(OUTDIR, "functional_analysis_results.xlsx")
cols_doc = pd.DataFrame([
    ("list", "all = all changed genes; up / down = genes up- / down-regulated in ASD-donor offspring"),
    ("term / term_pretty", "gene-set name (MSigDB) / readable name"),
    ("source_id", "original identifier (GO ID, Reactome ID, WikiPathways ID, HPO ID ...)"),
    ("overlap", "number of input genes in the term"),
    ("term_size", "number of genes in the term (within the database background)"),
    ("input_genes_in_database", "number of input genes present in the database background"),
    ("background", "number of genes annotated in the database (universe of the test)"),
    ("fold_enrichment", "(overlap / input genes) / (term size / background)"),
    ("p_value", "one-sided hypergeometric test P(X >= overlap); 1 if overlap < 2"),
    ("FDR_BH", "Benjamini–Hochberg adjusted p-value over all tested terms of the database, per list"),
    ("genes / genes_up / genes_down", "overlapping genes; split by direction of change"),
], columns=["column", "meaning"])
with pd.ExcelWriter(xlsx) as w:
    summary.to_excel(w, sheet_name="summary", index=False)
    for name, t in full.items():
        t.to_excel(w, sheet_name=name[:31], index=False)
    de.to_excel(w, sheet_name="input_genes", index=False)
    cols_doc.to_excel(w, sheet_name="column_description", index=False)
summary.to_csv(os.path.join(OUTDIR, "summary.tsv"), sep="\t", index=False)
print("Tables written to", xlsx)

# ------------------------------------------------------------------ figures
# (1) overview: -log10(FDR) of the top term per database (all genes)
fig, ax = plt.subplots(figsize=(6.5, 3.6))
v = -np.log10(summary.all_top_FDR.clip(lower=1e-10))
cols = ["#2F5D8C" if k > 0 else "#B8C2CC" for k in summary.all_n_significant]
ax.barh(range(len(summary))[::-1], v, color=cols)
ax.set_yticks(range(len(summary))[::-1]); ax.set_yticklabels([f"{a}  ({b} sig.)" for a, b in zip(summary.database, summary.all_n_significant)], fontsize=9)
ax.axvline(-np.log10(FDR_CUT), color="#B2182B", ls="--", lw=1); ax.set_xlabel("−log10(FDR) of top term")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); [fig.savefig(os.path.join(OUTDIR, "overview_databases" + e), dpi=300) for e in (".png", ".pdf")]; plt.close(fig)
# (2) GO Cellular Component: significant terms (all genes)
cc = full["GO Cellular Component"]; cc = cc[(cc.list == "all") & (cc.FDR_BH < FDR_CUT)].sort_values("p_value")
if len(cc):
    fig, ax = plt.subplots(figsize=(6.5, 0.35 * len(cc) + 1.2))
    ax.barh(range(len(cc))[::-1], -np.log10(cc.FDR_BH), color="#5B7FA6")
    ax.set_yticks(range(len(cc))[::-1]); ax.set_yticklabels([f"{t} ({k})" for t, k in zip(cc.term_pretty, cc.overlap)], fontsize=9)
    ax.axvline(-np.log10(FDR_CUT), color="#B2182B", ls="--", lw=1); ax.set_xlabel("−log10(FDR)")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.tight_layout(); [fig.savefig(os.path.join(OUTDIR, "GO_CC_significant_terms" + e), dpi=300) for e in (".png", ".pdf")]; plt.close(fig)
    uniq = sorted(set(",".join(cc.genes).split(",")))
    print(f"GO:CC – {len(cc)} significant terms, {len(uniq)} unique genes: {', '.join(uniq)}")

# ------------------------------------------------------------------ run info
with open(os.path.join(OUTDIR, "run_info.txt"), "w") as f:
    f.write(f"date: {datetime.datetime.now():%Y-%m-%d %H:%M}\ninput: {INPUT}\n")
    f.write(f"gene sets: MSigDB v2026.1 ({MSIGDB_URL}), git blob {git_blob_hash(MSIGDB_FILE)}\n")
    f.write(f"min set size: {MIN_SIZE}, max set size: {MAX_SIZE}, min overlap: {MIN_OVERLAP}, FDR cut-off: {FDR_CUT}\n")
    f.write(f"python {platform.python_version()}, pandas {pd.__version__}, numpy {np.__version__}, scipy {scipy.__version__}, pyreadr {pyreadr.__version__}, matplotlib {matplotlib.__version__}\n")
print("Done. Results in", OUTDIR)
