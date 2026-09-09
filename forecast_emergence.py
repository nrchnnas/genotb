"""
GenoTB-Forecast: pre-resistance / resistance-amplification signatures.

Honest framing: CRyPTIC is cross-sectional, so this is NOT a prospective time
forecast. It tests whether a strain's GENOMIC BACKGROUND — excluding the drug's
own target gene — anticipates resistance to that drug. A signal there means the
strain is evolutionarily 'primed', which is the surveillance-actionable proxy for
emergence risk. Two experiments:

  (A) Background-only anticipation: drop drug T's target-gene features, predict T.
      Compare AUROC to the full (target-gene) model.
  (B) Cryptic/pre-resistance: among isolates the catalogue calls SUSCEPTIBLE for T
      (no target mutation), can background still flag the truly-resistant ones?

Plus an evolvability score: do background features (compensatory rpoC, mutators,
efflux, first-line resistance breadth) rise with resistance breadth?
"""
import json, numpy as np, pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score

G = pd.read_csv("cryptic/aa_features_v3.csv").set_index("ENA_RUN")
Y = pd.read_csv("cryptic/pheno_full.csv").set_index("ENA_RUN")
df = G.join(Y, how="inner")
FEATURES = list(G.columns)

# drug -> its own direct target-gene features (to EXCLUDE for background model)
TARGETS = {
 "LEV": ["gyrA_QRDR","gyrA_D94","gyrA_A90","gyrB_QRDR"],
 "MXF": ["gyrA_QRDR","gyrA_D94","gyrA_A90","gyrB_QRDR"],
 "AMI": ["rrs_1401","rrs_1484","eis_promoter"],
 "KAN": ["rrs_1401","rrs_1484","eis_promoter"],
 "EMB": ["embB_306","embB_406_497","embB_354_328"],
 "ETH": ["ethA_LoF","ethA_nonsyn","inhA_promoter","inhA_I194_S94"],
}
NAME={"LEV":"levofloxacin","MXF":"moxifloxacin","AMI":"amikacin","KAN":"kanamycin",
      "EMB":"ethambutol","ETH":"ethionamide"}

def fit_auc(X, y):
    Xtr,Xte,ytr,yte = train_test_split(X,y,test_size=0.30,random_state=7,stratify=y)
    clf=GradientBoostingClassifier(n_estimators=200,max_depth=3,random_state=7).fit(Xtr,ytr)
    return roc_auc_score(yte, clf.predict_proba(Xte)[:,1])

rows=[]
for d3,dl in NAME.items():
    y=df[d3].map({"R":1,"S":0}); keep=y.notna()
    yy=y[keep].astype(int);
    if yy.sum()<50: continue
    Xfull=G.loc[keep[keep].index]
    bg_cols=[c for c in FEATURES if c not in TARGETS[d3]]
    auc_full=fit_auc(Xfull, yy)
    auc_bg  =fit_auc(Xfull[bg_cols], yy)
    # (B) cryptic resistance: isolates with NO target mutation
    no_target=(Xfull[TARGETS[d3]].sum(axis=1)==0)
    sub=Xfull[no_target]; ysub=yy[no_target]
    cryptic=""
    if ysub.nunique()==2 and ysub.sum()>=25:
        auc_cr=fit_auc(sub[bg_cols], ysub)
        cryptic=f"{auc_cr:.3f} (n_R={int(ysub.sum())})"
    rows.append(dict(drug=dl, n=int(keep.sum()), pctR=round(100*yy.mean(),1),
                     AUROC_full=round(auc_full,3), AUROC_background_only=round(auc_bg,3),
                     cryptic_resistance_AUROC=cryptic))

res=pd.DataFrame(rows)
res.to_csv("cryptic/forecast_metrics.csv", index=False)
pd.set_option("display.width",160)
print("=== (A,B) Anticipating resistance from background (no target gene) ===")
print(res.to_string(index=False))

# ---- evolvability score: background features vs resistance breadth ----
firstline=["INH","RIF","EMB"]
breadth=df[firstline].apply(lambda c:c.map({"R":1,"S":0})).sum(axis=1)
# classify
def cat(row):
    inh=row.get("INH")=="R"; rif=row.get("RIF")=="R"
    fq=(row.get("LEV")=="R") or (row.get("MXF")=="R")
    inj=(row.get("AMI")=="R") or (row.get("KAN")=="R")
    if inh and rif and (fq and inj): return "XDR"
    if inh and rif and (fq or inj): return "pre-XDR"
    if inh and rif: return "MDR"
    if inh or rif: return "mono/poly"
    return "pan-S"
df["_cat"]=df.apply(cat,axis=1)
bg=["rpoC_comp","mutator","efflux_Tap","whiB7_var"]
order=["pan-S","mono/poly","MDR","pre-XDR","XDR"]
eg=df.groupby("_cat")[bg].mean().reindex(order).round(3)
eg.to_csv("cryptic/evolvability_by_class.csv")
print("\n=== Evolvability background features by resistance class (mean prevalence) ===")
print(eg.to_string())
print("\nClass counts:", df["_cat"].value_counts().reindex(order).to_dict())
