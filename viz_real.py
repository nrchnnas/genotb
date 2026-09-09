"""Real-data figures (CRyPTIC). Overwrites figures/ so the proposal picks them up."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

os.makedirs("figures", exist_ok=True)
plt.rcParams.update({"font.size":10,"axes.titleweight":"bold","figure.dpi":130,"savefig.bbox":"tight"})
INK,RES,SUS,MISS,MID = "#1b2a4a","#c0392b","#2e86c1","#e9ecef","#e08a1e"

geno = pd.read_csv("cryptic/genotypes_real.csv").set_index("ENA_RUN")
phe  = pd.read_csv("cryptic/pheno_subset.csv").set_index("ENA_RUN")
met  = pd.read_csv("cryptic/model_metrics_real.csv", index_col=0)
impm = pd.read_csv("cryptic/importance_matrix_real.csv", index_col=0)
DRUGS = list(met.index)
NAME3 = {"isoniazid":"INH","rifampin":"RIF","ethambutol":"EMB","ethionamide":"ETH",
         "levofloxacin":"LEV","moxifloxacin":"MXF","amikacin":"AMI","kanamycin":"KAN"}
phe = phe.rename(columns={v:k for k,v in NAME3.items()})  # 3-letter -> long
phmat = phe[[d for d in DRUGS if d in phe.columns]].applymap(
    lambda x: 1 if x=="R" else (0 if x=="S" else np.nan))

# ---------- 1. resistance heatmap (real phenotypes) ----------
sub = phmat.dropna(thresh=6).sample(min(90,len(phmat)), random_state=3)
order = np.argsort(np.nansum(sub.values,axis=1))[::-1]
fill = sub.fillna(-1).values[order]
cmap = ListedColormap([MISS,SUS,RES])
fig,ax=plt.subplots(figsize=(7.6,9))
ax.imshow(np.where(fill<0,0,fill+1),cmap=cmap,aspect="auto",vmin=0,vmax=2,interpolation="nearest")
ax.set_xticks(range(len(sub.columns))); ax.set_xticklabels([NAME3[d] for d in sub.columns])
ax.set_yticks([]); ax.set_ylabel(f"{len(sub)} CRyPTIC M. tuberculosis isolates")
ax.set_title("Real resistance profile (CRyPTIC)\nphenotypes from 96-well MIC plates",loc="left")
ax.legend(handles=[Patch(color=RES,label="Resistant"),Patch(color=SUS,label="Susceptible"),
    Patch(color=MISS,label="Not tested")],loc="upper center",bbox_to_anchor=(0.5,-0.04),ncol=3,frameon=False)
plt.savefig("figures/01_resistance_heatmap.png"); plt.close()

# ---------- 2. co-resistance ----------
co=np.zeros((len(DRUGS),len(DRUGS)))
for i,a in enumerate(DRUGS):
    for j,b in enumerate(DRUGS):
        if a==b: co[i,j]=1; continue
        bb=phmat[[a,b]].dropna()
        if len(bb): co[i,j]=((bb[a]==1)&(bb[b]==1)).sum()/max(int((bb[a]==1).sum()),1)
fig,ax=plt.subplots(figsize=(7.4,6.4))
im=ax.imshow(co,cmap="magma_r",vmin=0,vmax=1)
ax.set_xticks(range(len(DRUGS))); ax.set_xticklabels([NAME3[d] for d in DRUGS],rotation=90)
ax.set_yticks(range(len(DRUGS))); ax.set_yticklabels([d.capitalize() for d in DRUGS])
for i in range(len(DRUGS)):
    for j in range(len(DRUGS)):
        ax.text(j,i,f"{co[i,j]:.2f}",ha="center",va="center",
                color="white" if co[i,j]>0.5 else INK,fontsize=7)
ax.set_title("Co-resistance (real CRyPTIC)  —  P(resistant col | resistant row)",loc="left")
fig.colorbar(im,fraction=0.046,pad=0.04); plt.savefig("figures/02_coresistance_map.png"); plt.close()

# ---------- 3. model performance ----------
fig,ax=plt.subplots(figsize=(9.5,5.2))
order=met.sort_values("AUROC").index; y=np.arange(len(order))
ax.barh(y,met.loc[order,"AUROC"],color=INK,height=0.6)
for i,d in enumerate(order):
    ax.text(met.loc[d,"AUROC"]+.004,i,f"{met.loc[d,'AUROC']:.2f}  (VME {met.loc[d,'very_major_error']:.0%})",
            va="center",fontsize=8)
ax.set_yticks(y); ax.set_yticklabels([d.capitalize() for d in order])
ax.set_xlim(0.8,1.0); ax.set_xlabel("AUROC (held-out 30% test set)"); ax.axvline(0.90,color=RES,ls="--",lw=1)
ax.set_title("REAL genotype→phenotype performance (CRyPTIC variant calls)\n"
    f"Macro-AUROC {met.AUROC.mean():.3f} · mean specificity {met.specificity.mean():.0%} · {len(geno)} isolates",loc="left")
plt.savefig("figures/03_model_performance.png"); plt.close()

# ---------- 4. transmission tree (real genotype distance) ----------
samp = phmat.dropna(subset=["isoniazid","rifampin"])
samp = samp[~samp.index.duplicated()].sample(min(60,len(samp)),random_state=11)
Gs = geno.reindex(samp.index)
D=pdist(Gs.values,metric="hamming"); Z=linkage(D,method="average")
mdr=((samp["isoniazid"]==1)&(samp["rifampin"]==1)).values
fig,ax=plt.subplots(figsize=(8.6,9.3))
dend=dendrogram(Z,orientation="left",no_labels=True,color_threshold=0.0,above_threshold_color=INK,ax=ax)
ylocs=np.arange(5,10*len(dend["leaves"])+5,10)
for yk,idx in zip(ylocs,dend["leaves"]):
    ax.plot(-0.002,yk,"o",color=RES if mdr[idx] else SUS,ms=6,clip_on=False)
ax.set_title("Real genotype-distance phylogeny (60 CRyPTIC isolates)\nred = MDR (INH+RIF) · blue = non-MDR",loc="left")
ax.set_xlabel("Hamming distance over real resistance-marker genotype")
ax.spines[["top","right"]].set_visible(False); ax.set_yticks([])
plt.savefig("figures/04_transmission_tree.png"); plt.close()

# ---------- 5. feature-importance heatmap (biology recovery) ----------
# show the key resistance features in a readable order
key=["hot_rpoB_RRDR","hot_katG_315","hot_inhA_promoter","gene_fabG1","gene_ethA",
     "hot_embB_306","gene_embB","hot_gyrA_QRDR","gene_gyrB","hot_rrs_1401","hot_eis_promoter","gene_pncA"]
key=[k for k in key if k in impm.index]
M=impm.loc[key,DRUGS]
fig,ax=plt.subplots(figsize=(8.8,7.2))
im=ax.imshow(M.values,cmap="YlGnBu",aspect="auto",vmin=0,vmax=float(M.values.max()))
ax.set_xticks(range(len(DRUGS))); ax.set_xticklabels([d.capitalize() for d in DRUGS],rotation=40,ha="right")
ax.set_yticks(range(len(key))); ax.set_yticklabels([k.replace("hot_","").replace("gene_","") for k in key])
for i in range(len(key)):
    for j in range(len(DRUGS)):
        v=M.values[i,j]
        if v>0.05: ax.text(j,i,f"{v:.2f}",ha="center",va="center",fontsize=7,
                           color="white" if v>M.values.max()*0.5 else INK)
ax.set_title("What the model learned: feature importance by drug\n"
             "(the model recovers the known resistance gene for each drug, from real data)",loc="left")
fig.colorbar(im,fraction=0.046,pad=0.04,label="Gini importance")
plt.savefig("figures/05_feature_importance.png"); plt.close()

print("Saved real figures:", sorted(os.listdir("figures")))
print("Macro-AUROC %.3f, isolates %d"%(met.AUROC.mean(), len(geno)))
