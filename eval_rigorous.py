"""
Independent re-evaluation of the "does resampling / generative augmentation help
rare-variant resistance prediction?" question.

Fixes three defects in the earlier comparison:
  1. It scored sensitivity at a FIXED 0.5 threshold. Oversampling mechanically
     inflates minority probabilities, so that metric rewards miscalibration, not
     discrimination. We report threshold-free metrics (AUROC, average precision)
     AND a matched operating point (sensitivity at specificity = 0.95).
  2. It omitted the trivial control: class_weight='balanced', which achieves the
     same threshold shift for free. If SMOTE only matches this, SMOTE adds nothing.
  3. It used a single split/seed. Differences of ~0.02 on a few hundred positives
     are noise. We use repeated stratified splits and report mean +/- 95% CI.

Resumable: appends per (drug, seed, method) to cryptic/eval_rigorous.csv.
"""
import os, time, numpy as np, pandas as pd, scipy.sparse as sp
import torch, torch.nn as nn
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, confusion_matrix
from imblearn.over_sampling import SMOTE

OUT="cryptic/eval_rigorous.csv"; BUDGET=140; START=time.time()
M=sp.load_npz("cryptic/mut_matrix.npz")
rows=open("cryptic/mut_rows.txt").read().split("\n")
Y=pd.read_csv("cryptic/pheno_full.csv").set_index("ENA_RUN")
ridx={r:i for i,r in enumerate(rows)}
DRUGS=["AMI","KAN","MXF","ETH","LEV","EMB"]
SEEDS=[0,1,2,3,4]

class VAE(nn.Module):
    def __init__(s,d,h=256,z=32):
        super().__init__(); s.e1=nn.Linear(d,h); s.mu=nn.Linear(h,z); s.lv=nn.Linear(h,z)
        s.d1=nn.Linear(z,h); s.out=nn.Linear(h,d); s.z=z
    def enc(s,x): h=torch.relu(s.e1(x)); return s.mu(h),s.lv(h)
    def dec(s,z): return torch.sigmoid(s.out(torch.relu(s.d1(z))))
    def forward(s,x):
        mu,lv=s.enc(x); return s.dec(mu+torch.exp(0.5*lv)*torch.randn_like(lv)),mu,lv

def vae_sample(Xpos,n,seed,epochs=40):
    torch.manual_seed(seed)
    X=torch.tensor(Xpos,dtype=torch.float32)
    m=VAE(X.shape[1]); opt=torch.optim.Adam(m.parameters(),1e-3)
    for _ in range(epochs):
        perm=torch.randperm(len(X))
        for i in range(0,len(X),128):
            b=X[perm[i:i+128]]; rec,mu,lv=m(b)
            loss=(nn.functional.binary_cross_entropy(rec,b,reduction="sum")
                  -0.5*torch.sum(1+lv-mu.pow(2)-lv.exp()))/len(b)
            opt.zero_grad(); loss.backward(); opt.step()
    if n<=0: return np.zeros((0,Xpos.shape[1]),np.float32)
    with torch.no_grad(): p=m.dec(torch.randn(n,m.z)).numpy()
    return (np.random.rand(*p.shape)<p).astype(np.float32)

def sens_at_spec(y,p,target=0.95):
    fpr,tpr,_=roc_curve(y,p); ok=fpr<=(1-target)
    return float(tpr[ok].max()) if ok.any() else 0.0

def score(y,p):
    pred=(p>=.5).astype(int)
    tn,fp,fn,tp=confusion_matrix(y,pred,labels=[0,1]).ravel()
    return dict(AUROC=roc_auc_score(y,p), AP=average_precision_score(y,p),
                sens_at_spec95=sens_at_spec(y,p),
                sens_at_0p5=(tp/(tp+fn) if tp+fn else np.nan),
                spec_at_0p5=(tn/(tn+fp) if tn+fp else np.nan))

done=set()
if os.path.exists(OUT):
    d0=pd.read_csv(OUT)
    done={(r.drug,r.seed,r.method) for r in d0.itertuples()}
else:
    pd.DataFrame(columns=["drug","seed","method","n","pctR","AUROC","AP",
        "sens_at_spec95","sens_at_0p5","spec_at_0p5"]).to_csv(OUT,index=False)

recs=[]
stop=False
for d in DRUGS:
    if stop: break
    yy=Y[d].map({"R":1,"S":0}).dropna()
    keep=[r for r in yy.index if r in ridx]; idx=np.array([ridx[r] for r in keep])
    Xd=M[idx].toarray().astype(np.float32); yd=yy.loc[keep].astype(int).values
    for seed in SEEDS:
        if time.time()-START>BUDGET: stop=True; break
        need=[m for m in ["baseline","class_weight","smote","vae"] if (d,seed,m) not in done]
        if not need: continue
        Xtr,Xte,ytr,yte=train_test_split(Xd,yd,test_size=0.30,random_state=seed,stratify=yd)
        base=dict(drug=d,seed=seed,n=len(yd),pctR=round(100*yd.mean(),1))
        for meth in need:
            try:
                if meth=="baseline":
                    mdl=LogisticRegression(max_iter=2000).fit(Xtr,ytr); p=mdl.predict_proba(Xte)[:,1]
                elif meth=="class_weight":
                    mdl=LogisticRegression(max_iter=2000,class_weight="balanced").fit(Xtr,ytr)
                    p=mdl.predict_proba(Xte)[:,1]
                elif meth=="smote":
                    Xr,yr=SMOTE(random_state=seed,k_neighbors=5).fit_resample(Xtr,ytr)
                    mdl=LogisticRegression(max_iter=2000).fit(Xr,yr); p=mdl.predict_proba(Xte)[:,1]
                else:
                    syn=vae_sample(Xtr[ytr==1], int((ytr==0).sum()-(ytr==1).sum()), seed)
                    Xa=np.vstack([Xtr,syn]); ya=np.concatenate([ytr,np.ones(len(syn),int)])
                    mdl=LogisticRegression(max_iter=2000).fit(Xa,ya); p=mdl.predict_proba(Xte)[:,1]
                recs.append({**base,"method":meth,**score(yte,p)})
            except Exception as e:
                print("ERR",d,seed,meth,e,flush=True)
        # append incrementally so partial progress survives a timeout
        if recs:
            COLS=["drug","seed","method","n","pctR","AUROC","AP","sens_at_spec95","sens_at_0p5","spec_at_0p5"]
            pd.DataFrame(recs)[COLS].to_csv(OUT,mode="a",header=False,index=False)
            recs=[]
        print(f"{d} seed{seed} done ({int(time.time()-START)}s)",flush=True)
tot=pd.read_csv(OUT)
print(f"\nrows so far: {len(tot)} | combos done: {tot.groupby(['drug','method']).size().to_dict()}")
print("DONE" if len(tot)>=len(DRUGS)*len(SEEDS)*4 else "MORE")
