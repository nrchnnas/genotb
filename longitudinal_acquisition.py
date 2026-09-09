"""
Within-host longitudinal pre-resistance analysis (CRyPTIC serial isolates).

Among 229 patients with >1 sequenced isolate, we ask a genuinely temporal/within-host
question: given a SUSCEPTIBLE isolate, does its genomic background (excluding the
drug's own target gene) predict that the SAME PATIENT harbours/acquires a RESISTANT
isolate for that drug?  Unit = (patient, drug) with >=1 susceptible isolate;
label = 1 if the patient also has a resistant isolate (acquisition/discordance),
else 0 (stable-susceptible). Features come from the patient's susceptible isolate.

Evaluation is grouped by patient (GroupKFold) so no patient appears in both train
and test, with a patient-level bootstrap 95% CI on the AUROC.

Caveat (stated in the manuscript): without collection dates, within-host
discordance may reflect de novo acquisition, heteroresistance, or mixed infection;
this is a within-host emergence signal, not a dated phylogenetic forecast.
"""
import re, csv, collections, numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from scipy.stats import fisher_exact

G = pd.read_csv("cryptic/aa_features_v3.csv").set_index("ENA_RUN")
reuse = pd.read_csv("cryptic/reuse.csv")
FEATURES = list(G.columns)
DRUGS = ["INH","RIF","EMB","LEV","MXF","AMI","KAN","ETH","RFB"]
TARGET = {
 "INH":["katG_S315","katG_LoF","inhA_promoter","inhA_I194_S94"],
 "RIF":["rpoB_S450L","rpoB_450_any","rpoB_RRDR","rpoB_435_445"],
 "RFB":["rpoB_S450L","rpoB_450_any","rpoB_RRDR","rpoB_435_445"],
 "EMB":["embB_306","embB_406_497","embB_354_328"],
 "ETH":["ethA_LoF","ethA_nonsyn","inhA_promoter","inhA_I194_S94"],
 "LEV":["gyrA_QRDR","gyrA_D94","gyrA_A90","gyrB_QRDR"],
 "MXF":["gyrA_QRDR","gyrA_D94","gyrA_A90","gyrB_QRDR"],
 "AMI":["rrs_1401","rrs_1484","eis_promoter"],
 "KAN":["rrs_1401","rrs_1484","eis_promoter"],
}
def parse(u):
    m=re.search(r"(site\.\d+\.subj(?:ect)?\.[^.]+).*?iso\.?(\d+)", u)
    return (m.group(1),int(m.group(2))) if m else (None,None)

reuse["subj"],reuse["isoN"]=zip(*reuse["UNIQUEID"].map(parse))
multi=reuse[reuse["subj"].notna()].groupby("subj").filter(lambda d: d["ENA_RUN"].nunique()>1)

# build (patient, drug) units from susceptible isolates
rows=[]
for subj,grp in multi.groupby("subj"):
    grp=grp[grp["ENA_RUN"].isin(G.index)]
    if grp["ENA_RUN"].nunique()<2: continue
    for d in DRUGS:
        col=d+"_BINARY_PHENOTYPE"
        ph={r["ENA_RUN"]:r[col] for _,r in grp.iterrows()}
        S=[run for run,p in ph.items() if p=="S"]
        anyR=any(p=="R" for p in ph.values())
        anyS=len(S)>0
        if not anyS: continue                # need a susceptible isolate to predict from
        rep=sorted(S, key=lambda run:int(grp.loc[grp.ENA_RUN==run,"isoN"].iloc[0]))[0]
        feat=G.loc[rep]
        if isinstance(feat,pd.DataFrame): feat=feat.iloc[0]
        # background = all features except this drug's target genes
        bg={c:int(feat[c]) for c in FEATURES}      # full background
        for tcol in TARGET[d]: bg[tcol]=0           # zero this drug's own target genes
        rows.append(dict(subj=subj, drug=d, label=int(anyR), **bg))

df=pd.DataFrame(rows).fillna(0)
# add drug one-hot
df=pd.concat([df, pd.get_dummies(df["drug"],prefix="drug")],axis=1)
feat_cols=[c for c in df.columns if c not in ("subj","drug","label")]
X=df[feat_cols].values; y=df["label"].values; groups=df["subj"].values
print("units (patient,drug):",len(df)," acquisition-positive:",int(y.sum()),
      " patients:",df['subj'].nunique())

# grouped CV out-of-fold predictions
oof=np.zeros(len(y));
gkf=GroupKFold(n_splits=5)
for tr,te in gkf.split(X,y,groups):
    clf=GradientBoostingClassifier(n_estimators=150,max_depth=2,random_state=7).fit(X[tr],y[tr])
    oof[te]=clf.predict_proba(X[te])[:,1]
auc=roc_auc_score(y,oof)
# patient-level bootstrap CI
rng=np.random.default_rng(0); pts=df["subj"].unique(); aucs=[]
idx_by_pt={p:np.where(groups==p)[0] for p in pts}
for _ in range(2000):
    samp=rng.choice(pts,size=len(pts),replace=True)
    ii=np.concatenate([idx_by_pt[p] for p in samp])
    if len(np.unique(y[ii]))<2: continue
    aucs.append(roc_auc_score(y[ii],oof[ii]))
lo,hi=np.percentile(aucs,[2.5,97.5])
print(f"\nWithin-host acquisition forecast (background, target gene excluded):")
print(f"  grouped-CV AUROC = {auc:.3f}  (95% bootstrap CI {lo:.3f}-{hi:.3f})")

# descriptive: odds ratios for key background markers (acquisition vs stable)
print("\nBackground marker enrichment in susceptible isolates that later show resistance:")
for m in ["rpoC_comp","mutator","efflux_Tap","whiB7_var"]:
    if m not in df: continue
    a=((df[m]==1)&(df.label==1)).sum(); b=((df[m]==1)&(df.label==0)).sum()
    c=((df[m]==0)&(df.label==1)).sum(); e=((df[m]==0)&(df.label==0)).sum()
    orr,p=fisher_exact([[a,b],[c,e]])
    print(f"  {m:12s} OR={orr:4.2f}  p={p:.3f}  (acq:{a}/{a+c}  stable:{b}/{b+e})")

pd.DataFrame(dict(auc=[auc],lo=[lo],hi=[hi],n=[len(df)],pos=[int(y.sum())],
                  patients=[df['subj'].nunique()])).to_csv("cryptic/longitudinal_metrics.csv",index=False)
np.save("cryptic/_oof.npy",np.c_[y,oof])
print("\nsaved cryptic/longitudinal_metrics.csv")
