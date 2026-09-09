"""Visualizations for the AMR genome->resistance predictor (saves PNGs to figures/)."""
import json, os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

os.makedirs("figures", exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.titlweight" if False else "axes.titleweight": "bold",
                     "figure.dpi": 130, "savefig.bbox": "tight"})
INK, RES, SUS, MISS = "#1b2a4a", "#c0392b", "#2e86c1", "#e9ecef"

mat = pd.read_csv("data/tb_matrix.csv", index_col=0)
geno = pd.read_csv("data/tb_genotypes.csv", index_col=0)
metrics = pd.read_csv("data/model_metrics.csv", index_col=0)
DRUGS = list(metrics.index)
mat = mat[DRUGS]

# ---------- 1. Resistance heatmap (isolates x drugs), clustered ----------
sub = mat.dropna(thresh=8).copy()
sub = sub.sample(min(90, len(sub)), random_state=3)
fill = sub.fillna(-1).values
# order isolates by total resistance for a clean gradient
order = np.argsort(np.nansum(sub.values, axis=1))[::-1]
fill = fill[order]
cmap = ListedColormap([MISS, SUS, RES])
fig, ax = plt.subplots(figsize=(8.6, 9))
ax.imshow(np.where(fill < 0, 0, fill + 1), cmap=cmap, aspect="auto",
          vmin=0, vmax=2, interpolation="nearest")
ax.set_xticks(range(len(DRUGS)))
ax.set_xticklabels([d[:4].upper() for d in DRUGS], rotation=0)
ax.set_yticks([]); ax.set_ylabel(f"{len(sub)} M. tuberculosis isolates (BV-BRC)")
ax.set_title("In-silico resistance profile\nReal phenotypes, "
             f"{len(DRUGS)} drugs", loc="left")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=RES, label="Resistant"), Patch(color=SUS, label="Susceptible"),
                   Patch(color=MISS, label="Not tested")],
          loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=3, frameon=False)
plt.savefig("figures/01_resistance_heatmap.png"); plt.close()

# ---------- 2. Co-resistance matrix (Jaccard among Resistant calls) ----------
M = mat.copy()
co = np.zeros((len(DRUGS), len(DRUGS)))
for i, a in enumerate(DRUGS):
    for j, b in enumerate(DRUGS):
        if a == b:
            co[i, j] = 1.0; continue
        both = M[[a, b]].dropna()
        if len(both):
            ra = (both[a] == 1)
            co[i, j] = ((both[a] == 1) & (both[b] == 1)).sum() / max(int(ra.sum()), 1)
fig, ax = plt.subplots(figsize=(8.5, 7.2))
im = ax.imshow(co, cmap="rocket_r" if False else "magma_r", vmin=0, vmax=1)
ax.set_xticks(range(len(DRUGS))); ax.set_xticklabels([d[:4].upper() for d in DRUGS], rotation=90)
ax.set_yticks(range(len(DRUGS))); ax.set_yticklabels([d.capitalize() for d in DRUGS])
for i in range(len(DRUGS)):
    for j in range(len(DRUGS)):
        ax.text(j, i, f"{co[i,j]:.2f}", ha="center", va="center",
                color="white" if co[i, j] > 0.5 else INK, fontsize=7)
ax.set_title("Co-resistance map  —  P(resistant to column | resistant to row)", loc="left")
fig.colorbar(im, fraction=0.046, pad=0.04)
plt.savefig("figures/02_coresistance_map.png"); plt.close()

# ---------- 3. Model performance: AUROC + very-major-error ----------
fig, ax = plt.subplots(figsize=(9.5, 5.6))
order = metrics.sort_values("AUROC").index
y = np.arange(len(order))
ax.barh(y, metrics.loc[order, "AUROC"], color=INK, height=0.6)
for i, d in enumerate(order):
    ax.text(metrics.loc[d, "AUROC"] + .004, i,
            f"{metrics.loc[d,'AUROC']:.2f}  (VME {metrics.loc[d,'very_major_error']:.0%})",
            va="center", fontsize=8)
ax.set_yticks(y); ax.set_yticklabels([d.capitalize() for d in order])
ax.set_xlim(0.8, 1.0); ax.set_xlabel("AUROC (held-out 30% test set)")
ax.axvline(0.90, color=RES, ls="--", lw=1)
ax.set_title("Genome→resistance prediction performance per drug\n"
             f"Macro-AUROC {metrics.AUROC.mean():.3f} · mean specificity {metrics.specificity.mean():.0%}",
             loc="left")
plt.savefig("figures/03_model_performance.png"); plt.close()

# ---------- 4. Transmission / phylogenetic tree (genotype distance dendrogram) ----------
# Use a sample of isolates; distance = Hamming over genotype; color tips by MDR status.
samp = mat.dropna(subset=["isoniazid", "rifampin"])
samp = samp[~samp.index.duplicated()].sample(min(60, len(samp)), random_state=11)
G = geno[~geno.index.duplicated()].reindex(samp.index)
samp = samp.loc[G.index]            # keep perfectly aligned
D = pdist(G.values, metric="hamming")
Z = linkage(D, method="average")
mdr = ((samp["isoniazid"] == 1) & (samp["rifampin"] == 1)).values
fig, ax = plt.subplots(figsize=(8.8, 9.5))
dend = dendrogram(Z, orientation="left", no_labels=True, color_threshold=0.0,
                  above_threshold_color=INK, ax=ax)
# annotate leaves with MDR color dots
leaf_order = dend["leaves"]
ylocs = np.arange(5, 10 * len(leaf_order) + 5, 10)
for yk, idx in zip(ylocs, leaf_order):
    ax.plot(-0.002, yk, "o", color=RES if mdr[idx] else SUS, ms=6,
            clip_on=False)
ax.set_title("Genotype-distance phylogeny (60 isolates)\n"
             "red = MDR (INH+RIF resistant) · blue = non-MDR", loc="left")
ax.set_xlabel("Hamming distance over resistance-marker genotype")
ax.spines[["top", "right"]].set_visible(False); ax.set_yticks([])
plt.savefig("figures/04_transmission_tree.png"); plt.close()

print("Saved figures:")
for f in sorted(os.listdir("figures")):
    print("  figures/" + f)
