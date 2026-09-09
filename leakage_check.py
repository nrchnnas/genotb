"""
Pitfall 3: clonal / population-structure leakage.

M. tuberculosis is clonal: isolates from the same lineage or transmission cluster
are near-identical. A random train/test split therefore places near-duplicates on
both sides, and the model can succeed by recognising the clone rather than the
resistance mechanism. The honest evaluation blocks splits by population structure.

We approximate structure by projecting the mutation matrix (TruncatedSVD) and
clustering (KMeans), then compare:
   random   : ordinary stratified random split
   blocked  : cluster-blocked split (no cluster spans train and test)
Resumable; appends to cryptic/leakage_check.csv.
"""
import os, time, numpy as np, pandas as pd, scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans
from sklearn.model_selection import train_test_split, GroupShuffleSplit
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score

OUT="cryptic/leakage_check.csv"; BUDGET=140; START=time.time()
COLS=["drug","seed","split","AUROC","AP","n_test","n_pos_test"]
M=sp.load_npz("cryptic/mut_matrix.npz")
rows=open("cryptic/mut_rows.txt").read().split("\n")
Y=pd.read_csv("cryptic/pheno_full.csv").set_index("ENA_RUN")
ridx={r:i for i,r in enumerate(rows)}

# population-structure proxy: SVD -> KMeans clusters over ALL isolates
if os.path.exists("cryptic/clusters.npy"):
    clus=np.load("cryptic/clusters.npy")
else:
    Z=TruncatedSVD(n_components=50,random_state=0).fit_transform(M)
    clus=KMeans(n_clusters=60,random_state=0,n_init=4).fit_predict(Z)
    np.save("cryptic/clusters.npy",clus)
    print("clusters built; sizes:",np.bincount(clus)[:10],"...")

DRUGS=["INH","RIF","EMB","ETH","LEV","MXF","AMI","KAN"]
SEEDS=[0,1,2,3,4]
done=set()
if os.path.exists(OUT):
    d0=pd.read_csv(OUT); done={(r.drug,r.seed,r.split) for r in d0.itertuples()}
else:
    pd.DataFrame(columns=COLS).to_csv(OUT,index=False)

recs=[]; stop=False
for d in DRUGS:
    if stop: break
    yy=Y[d].map({"R":1,"S":0}).dropna()
    keep=[r for r in yy.index if r in ridx]; idx=np.array([ridx[r] for r in keep])
    X=M[idx].toarray().astype(np.float32); y=yy.loc[keep].astype(int).values
    g=clus[idx]
    for seed in SEEDS:
        if time.time()-START>BUDGET: stop=True; break
        if (d,seed,"random") in done and (d,seed,"blocked") in done: continue
        # random split
        tr,te=train_test_split(np.arange(len(y)),test_size=.3,random_state=seed,stratify=y)
        m=LogisticRegression(max_iter=2000).fit(X[tr],y[tr]); p=m.predict_proba(X[te])[:,1]
        recs.append(dict(drug=d,seed=seed,split="random",AUROC=roc_auc_score(y[te],p),
                         AP=average_precision_score(y[te],p),n_test=len(te),n_pos_test=int(y[te].sum())))
        # cluster-blocked split
        gss=GroupShuffleSplit(n_splits=1,test_size=.3,random_state=seed)
        tr2,te2=next(gss.split(X,y,groups=g))
        if len(set(y[te2]))>1 and len(set(y[tr2]))>1:
            m2=LogisticRegression(max_iter=2000).fit(X[tr2],y[tr2]); p2=m2.predict_proba(X[te2])[:,1]
            recs.append(dict(drug=d,seed=seed,split="blocked",AUROC=roc_auc_score(y[te2],p2),
                             AP=average_precision_score(y[te2],p2),n_test=len(te2),n_pos_test=int(y[te2].sum())))
        pd.DataFrame(recs)[COLS].to_csv(OUT,mode="a",header=False,index=False); recs=[]
        print(f"{d} seed{seed} done ({int(time.time()-START)}s)",flush=True)

t=pd.read_csv(OUT); print("\nrows:",len(t),"DONE" if len(t)>=len(DRUGS)*len(SEEDS)*2 else "MORE")
