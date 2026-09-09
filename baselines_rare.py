"""Baselines the generative approach must beat, on the high-dim mutation matrix.
Discriminative (logistic regression) vs SMOTE oversampling, for rare-resistance
drugs. Reports held-out AUROC / sensitivity / specificity / very-major error."""
import numpy as np, pandas as pd, scipy.sparse as sp
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, confusion_matrix
from imblearn.over_sampling import SMOTE

M=sp.load_npz("cryptic/mut_matrix.npz")
rows=open("cryptic/mut_rows.txt").read().split("\n")
Y=pd.read_csv("cryptic/pheno_full.csv").set_index("ENA_RUN")
ridx={r:i for i,r in enumerate(rows)}
DRUGS=["AMI","KAN","MXF","ETH","EMB","LEV","INH","RIF"]

def metrics(y,p):
    pred=(p>=.5).astype(int)
    tn,fp,fn,tp=confusion_matrix(y,pred,labels=[0,1]).ravel()
    sens=tp/(tp+fn) if tp+fn else np.nan; spec=tn/(tn+fp) if tn+fp else np.nan
    return roc_auc_score(y,p), sens, spec, (fn/(tp+fn) if tp+fn else np.nan)

out=[]
for d in DRUGS:
    yy=Y[d].map({"R":1,"S":0}); yy=yy.dropna()
    idx=[ridx[r] for r in yy.index if r in ridx]
    keep=[r for r in yy.index if r in ridx]
    X=M[idx]; y=yy.loc[keep].astype(int).values
    if y.sum()<40: continue
    Xtr,Xte,ytr,yte=train_test_split(X,y,test_size=0.30,random_state=7,stratify=y)
    # (1) discriminative baseline
    lr=LogisticRegression(max_iter=2000,C=1.0).fit(Xtr,ytr)
    a,s,sp_,v=metrics(yte,lr.predict_proba(Xte)[:,1])
    # (2) SMOTE oversampling
    try:
        Xr,yr=SMOTE(random_state=7,k_neighbors=min(5,int(ytr.sum())-1)).fit_resample(Xtr,ytr)
        lr2=LogisticRegression(max_iter=2000,C=1.0).fit(Xr,yr)
        a2,s2,sp2,v2=metrics(yte,lr2.predict_proba(Xte)[:,1])
    except Exception as e:
        a2=s2=sp2=v2=np.nan
    out.append(dict(drug=d, n=len(y), pctR=round(100*y.mean(),1),
        base_AUROC=round(a,3), base_sens=round(s,3), base_VME=round(v,3),
        smote_AUROC=round(a2,3), smote_sens=round(s2,3), smote_VME=round(v2,3)))

res=pd.DataFrame(out); res.to_csv("cryptic/gen_baselines.csv",index=False)
pd.set_option("display.width",160); print(res.to_string(index=False))
