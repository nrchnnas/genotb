"""
Is the 'background anticipates resistance' result pre-resistance, or just
co-resistance with MDR status?

For each drug we remove its own target genes (as in the original analysis) and
compare:
  bg          : full background model (what was reported)
  mdr_only    : a SINGLE binary feature, 'is this strain MDR (INH-R and RIF-R)?'
  bg_in_mdr   : background model evaluated only within MDR strains
  bg_in_nonmdr: background model evaluated only within non-MDR strains
If mdr_only approaches bg, the reported signal is largely confounding.
Resumable; appends to cryptic/confound_check.csv.
"""
import os, time, numpy as np, pandas as pd, scipy.sparse as sp
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

OUT="cryptic/confound_check.csv"; BUDGET=140; START=time.time()
COLS=["drug","seed","bg","mdr_only","bg_in_mdr","bg_in_nonmdr","n_test_mdr","n_test_non"]
M=sp.load_npz("cryptic/mut_matrix.npz")
rows=open("cryptic/mut_rows.txt").read().split("\n")
cols=open("cryptic/mut_cols.txt").read().split("\n")
Y=pd.read_csv("cryptic/pheno_full.csv").set_index("ENA_RUN")
ridx={r:i for i,r in enumerate(rows)}
TARGET={"MXF":("gyrA_","gyrB_"),"LEV":("gyrA_","gyrB_"),
        "AMI":("rrs_","eis_"),"KAN":("rrs_","eis_"),
        "EMB":("embB_","embA_"),"ETH":("ethA_","inhA_","fabG1_")}
SEEDS=[0,1,2,3,4]

done=set()
if os.path.exists(OUT):
    d0=pd.read_csv(OUT); done={(r.drug,r.seed) for r in d0.itertuples()}
else:
    pd.DataFrame(columns=COLS).to_csv(OUT,index=False)

recs=[]; stop=False
for d,pref in TARGET.items():
    if stop: break
    yy=Y[d].map({"R":1,"S":0}).dropna()
    keep=[r for r in yy.index if r in ridx]; idx=np.array([ridx[r] for r in keep])
    y=yy.loc[keep].astype(int).values
    bgcols=[j for j,c in enumerate(cols) if not c.startswith(pref)]
    X=M[idx][:,bgcols].toarray().astype(np.float32)
    inh=(Y["INH"].reindex(keep).values=="R"); rif=(Y["RIF"].reindex(keep).values=="R")
    mdr=(inh&rif).astype(int)
    for seed in SEEDS:
        if (d,seed) in done: continue
        if time.time()-START>BUDGET: stop=True; break
        tr,te=train_test_split(np.arange(len(y)),test_size=.3,random_state=seed,stratify=y)
        mdl=LogisticRegression(max_iter=2000).fit(X[tr],y[tr])
        p_te=mdl.predict_proba(X[te])[:,1]
        r={"drug":d,"seed":seed,"bg":roc_auc_score(y[te],p_te),
           "mdr_only":roc_auc_score(y[te],mdr[te])}
        for nm,mask in [("bg_in_mdr",mdr==1),("bg_in_nonmdr",mdr==0)]:
            sub=np.array([i for i in te if mask[i]])
            ok=len(sub)>30 and len(set(y[sub]))>1
            r[nm]=roc_auc_score(y[sub],mdl.predict_proba(X[sub])[:,1]) if ok else np.nan
            r["n_test_mdr" if nm=="bg_in_mdr" else "n_test_non"]=len(sub)
        recs.append(r)
        pd.DataFrame(recs)[COLS].to_csv(OUT,mode="a",header=False,index=False); recs=[]
        print(f"{d} seed{seed} bg={r['bg']:.3f} mdr_only={r['mdr_only']:.3f}",flush=True)

t=pd.read_csv(OUT)
print("\nrows:",len(t),"DONE" if len(t)>=len(TARGET)*len(SEEDS) else "MORE")
