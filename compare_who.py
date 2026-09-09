"""Compare codon-level features vs codon + per-allele WHO catalogue features."""
import numpy as np, pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, confusion_matrix

G  = pd.read_csv("cryptic/aa_features_v3.csv").set_index("ENA_RUN")
W  = pd.read_csv("cryptic/who_features.csv").set_index("ENA_RUN")
Y  = pd.read_csv("cryptic/pheno_full.csv").set_index("ENA_RUN")
GW = G.join(W, how="inner")
NAME={"INH":"isoniazid","RIF":"rifampin","EMB":"ethambutol","ETH":"ethionamide",
      "LEV":"levofloxacin","MXF":"moxifloxacin","AMI":"amikacin","KAN":"kanamycin","RFB":"rifabutin"}

def evaln(X,y):
    Xtr,Xte,ytr,yte=train_test_split(X,y,test_size=0.30,random_state=7,stratify=y)
    clf=GradientBoostingClassifier(n_estimators=250,max_depth=3,random_state=7).fit(Xtr,ytr)
    p=clf.predict_proba(Xte)[:,1]; pred=(p>=.5).astype(int)
    tn,fp,fn,tp=confusion_matrix(yte,pred,labels=[0,1]).ravel()
    return roc_auc_score(yte,p), (tp/(tp+fn) if tp+fn else np.nan), (fn/(tp+fn) if tp+fn else np.nan)

rows=[]
for d3,dl in NAME.items():
    y=Y.join(GW,how="inner")[d3].map({"R":1,"S":0}) if d3 in Y.columns else None
    yy=Y[d3].reindex(GW.index).map({"R":1,"S":0}); keep=yy.notna()
    yk=yy[keep].astype(int); idx=keep[keep].index
    if yk.nunique()<2 or yk.sum()<50: continue
    aB,sB,vB=evaln(G.loc[idx],yk)                       # codon-level only
    aW,sW,vW=evaln(GW.loc[idx],yk)                      # + per-allele WHO
    rows.append(dict(drug=dl, AUROC_codon=round(aB,3), AUROC_codon_plus_WHO=round(aW,3),
                     dAUROC=round(aW-aB,3), sens_codon=round(sB,3), sens_WHO=round(sW,3),
                     dSens=round(sW-sB,3)))
res=pd.DataFrame(rows); res.to_csv("cryptic/who_compare.csv",index=False)
pd.set_option("display.width",160); print(res.to_string(index=False))
print("\nMean AUROC: codon %.3f -> +WHO %.3f (Δ %+.3f)"%(
    res.AUROC_codon.mean(),res.AUROC_codon_plus_WHO.mean(),res.dAUROC.mean()))
print("Mean sensitivity: codon %.3f -> +WHO %.3f (Δ %+.3f)"%(
    res.sens_codon.mean(),res.sens_WHO.mean(),res.dSens.mean()))
