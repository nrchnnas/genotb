"""Final figures: full CRyPTIC cohort (12,245) x amino-acid features. Overwrites figures/."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

os.makedirs("figures", exist_ok=True)
plt.rcParams.update({"font.size":10,"axes.titleweight":"bold","figure.dpi":130,"savefig.bbox":"tight"})
INK,RES,SUS,MISS = "#1b2a4a","#c0392b","#2e86c1","#e9ecef"

geno = pd.read_csv("cryptic/aa_features_v2.csv").set_index("ENA_RUN")
phe  = pd.read_csv("cryptic/pheno_full.csv").set_index("ENA_RUN")
met  = pd.read_csv("cryptic/model_metrics_full.csv", index_col=0)
impm = pd.read_csv("cryptic/importance_matrix_full.csv", index_col=0)
DRUGS=list(met.index)
N3={"isoniazid":"INH","rifampin":"RIF","ethambutol":"EMB","ethionamide":"ETH","levofloxacin":"LEV",
    "moxifloxacin":"MXF","amikacin":"AMI","kanamycin":"KAN","rifabutin":"RFB"}
phe=phe.rename(columns={v:k for k,v in N3.items()})
phmat=phe[[d for d in DRUGS if d in phe.columns]].apply(
    lambda col: col.map(lambda x:1 if x=="R" else (0 if x=="S" else np.nan)))

# 1. resistance heatmap
sub=phmat.dropna(thresh=7).sample(min(90,len(phmat)),random_state=3)
order=np.argsort(np.nansum(sub.values,axis=1))[::-1]; fill=sub.fillna(-1).values[order]
fig,ax=plt.subplots(figsize=(7.8,9))
ax.imshow(np.where(fill<0,0,fill+1),cmap=ListedColormap([MISS,SUS,RES]),aspect="auto",vmin=0,vmax=2,interpolation="nearest")
ax.set_xticks(range(len(sub.columns))); ax.set_xticklabels([N3[d] for d in sub.columns])
ax.set_yticks([]); ax.set_ylabel(f"{len(sub)} of 12,245 CRyPTIC isolates")
ax.set_title("Susceptibility profile, representative isolates",loc="left")
ax.legend(handles=[Patch(color=RES,label="Resistant"),Patch(color=SUS,label="Susceptible"),Patch(color=MISS,label="Not tested")],
    loc="upper center",bbox_to_anchor=(0.5,-0.04),ncol=3,frameon=False)
plt.savefig("figures/01_resistance_heatmap.png"); plt.close()

# 2. co-resistance
co=np.zeros((len(DRUGS),len(DRUGS)))
for i,a in enumerate(DRUGS):
    for j,b in enumerate(DRUGS):
        if a==b: co[i,j]=1; continue
        bb=phmat[[a,b]].dropna()
        if len(bb): co[i,j]=((bb[a]==1)&(bb[b]==1)).sum()/max(int((bb[a]==1).sum()),1)
fig,ax=plt.subplots(figsize=(7.8,6.8))
im=ax.imshow(co,cmap="magma_r",vmin=0,vmax=1)
ax.set_xticks(range(len(DRUGS))); ax.set_xticklabels([N3[d] for d in DRUGS],rotation=90)
ax.set_yticks(range(len(DRUGS))); ax.set_yticklabels([d.capitalize() for d in DRUGS])
for i in range(len(DRUGS)):
    for j in range(len(DRUGS)):
        ax.text(j,i,f"{co[i,j]:.2f}",ha="center",va="center",color="white" if co[i,j]>0.5 else INK,fontsize=6.5)
ax.set_title("Co-resistance: P(resistant to column | resistant to row)",loc="left")
fig.colorbar(im,fraction=0.046,pad=0.04); plt.savefig("figures/02_coresistance_map.png"); plt.close()

# 3. performance
fig,ax=plt.subplots(figsize=(9.5,5.6))
order=met.sort_values("AUROC").index; y=np.arange(len(order))
ax.barh(y,met.loc[order,"AUROC"],color=INK,height=0.62)
for i,d in enumerate(order):
    ax.text(met.loc[d,"AUROC"]+.004,i,f"{met.loc[d,'AUROC']:.2f}  (VME {met.loc[d,'very_major_error']:.0%})",va="center",fontsize=8)
ax.set_yticks(y); ax.set_yticklabels([d.capitalize() for d in order])
ax.set_xlim(0.8,1.0); ax.set_xlabel("AUROC (held-out 30% test set)"); ax.axvline(0.90,color=RES,ls="--",lw=1)
ax.set_title(""
    f"Macro-AUROC {met.AUROC.mean():.3f} · mean sensitivity {met.sensitivity.mean():.0%} · 12,245 isolates · 9 drugs",loc="left")
plt.savefig("figures/03_model_performance.png"); plt.close()

# 4. transmission tree
samp=phmat.dropna(subset=["isoniazid","rifampin"]); samp=samp[~samp.index.duplicated()].sample(min(60,len(samp)),random_state=11)
Gs=geno.reindex(samp.index)
Z=linkage(pdist(Gs.values,metric="hamming"),method="average")
mdr=((samp["isoniazid"]==1)&(samp["rifampin"]==1)).values
fig,ax=plt.subplots(figsize=(8.6,9.3))
dend=dendrogram(Z,orientation="left",no_labels=True,color_threshold=0.0,above_threshold_color=INK,ax=ax)
for yk,idx in zip(np.arange(5,10*len(dend["leaves"])+5,10),dend["leaves"]):
    ax.plot(-0.002,yk,"o",color=RES if mdr[idx] else SUS,ms=6,clip_on=False)
ax.set_title("Genotype-distance phylogeny (60 isolates); red, multidrug-resistant",loc="left")
ax.set_xlabel("Hamming distance over amino-acid resistance genotype"); ax.spines[["top","right"]].set_visible(False); ax.set_yticks([])
plt.savefig("figures/04_transmission_tree.png"); plt.close()

# 5. feature importance (AA resolution)
key=["rpoB_RRDR","rpoB_S450L","rpoB_435_445","katG_S315","inhA_promoter","ethA_LoF","ethA_nonsyn",
     "embB_306","embB_406_497","gyrA_QRDR","gyrA_D94","gyrB_QRDR","rrs_1401","rrs_1484","eis_promoter","rpsL_43_88"]
key=[k for k in key if k in impm.index]
M=impm.loc[key,DRUGS]
fig,ax=plt.subplots(figsize=(9.4,7.6))
im=ax.imshow(M.values,cmap="YlGnBu",aspect="auto",vmin=0,vmax=float(M.values.max()))
ax.set_xticks(range(len(DRUGS))); ax.set_xticklabels([d.capitalize() for d in DRUGS],rotation=40,ha="right")
ax.set_yticks(range(len(key))); ax.set_yticklabels(key)
for i in range(len(key)):
    for j in range(len(DRUGS)):
        v=M.values[i,j]
        if v>0.05: ax.text(j,i,f"{v:.2f}",ha="center",va="center",fontsize=6.5,color="white" if v>M.values.max()*0.5 else INK)
ax.set_title("Amino-acid feature importance by drug",loc="left")
fig.colorbar(im,fraction=0.046,pad=0.04,label="Gini importance"); plt.savefig("figures/05_feature_importance.png"); plt.close()

print("done. macro-AUROC %.3f, isolates %d, drugs %d"%(met.AUROC.mean(),len(geno),len(DRUGS)))
