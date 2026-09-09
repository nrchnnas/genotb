"""
Second generative use: a VAE latent REPRESENTATION learned unsupervised over all
12,245 isolates. Does encoding genotypes into a 32-d latent space help rare-drug
prediction (raw vs latent vs raw+latent)? A fair test of representation learning.
"""
import numpy as np, pandas as pd, scipy.sparse as sp, torch, torch.nn as nn
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, confusion_matrix
torch.manual_seed(0); np.random.seed(0)

M=sp.load_npz("cryptic/mut_matrix.npz"); rows=open("cryptic/mut_rows.txt").read().split("\n")
Y=pd.read_csv("cryptic/pheno_full.csv").set_index("ENA_RUN")
ridx={r:i for i,r in enumerate(rows)}; D=M.shape[1]
Xall=torch.tensor(M.toarray(),dtype=torch.float32)

class VAE(nn.Module):
    def __init__(s,d,h=256,z=32):
        super().__init__(); s.e1=nn.Linear(d,h); s.mu=nn.Linear(h,z); s.lv=nn.Linear(h,z)
        s.d1=nn.Linear(z,h); s.out=nn.Linear(h,d)
    def enc(s,x): h=torch.relu(s.e1(x)); return s.mu(h),s.lv(h)
    def dec(s,z): return torch.sigmoid(s.out(torch.relu(s.d1(z))))
    def forward(s,x):
        mu,lv=s.enc(x); z=mu+torch.exp(0.5*lv)*torch.randn_like(lv); return s.dec(z),mu,lv

m=VAE(D); opt=torch.optim.Adam(m.parameters(),1e-3)
for ep in range(50):
    perm=torch.randperm(len(Xall))
    for i in range(0,len(Xall),256):
        b=Xall[perm[i:i+256]]; rec,mu,lv=m(b)
        bce=nn.functional.binary_cross_entropy(rec,b,reduction="sum")
        kl=-0.5*torch.sum(1+lv-mu.pow(2)-lv.exp())
        loss=(bce+kl)/len(b); opt.zero_grad(); loss.backward(); opt.step()
with torch.no_grad(): Zall=m.enc(Xall)[0].numpy()      # latent means
print("VAE trained; latent shape",Zall.shape)

def metrics(y,p):
    pred=(p>=.5).astype(int); tn,fp,fn,tp=confusion_matrix(y,pred,labels=[0,1]).ravel()
    return roc_auc_score(y,p),(tp/(tp+fn) if tp+fn else np.nan),(fn/(tp+fn) if tp+fn else np.nan)

Mden=M.toarray().astype(np.float32)
DRUGS=["AMI","KAN","MXF","ETH","EMB","LEV","INH","RIF"]
out=[]
for d in DRUGS:
    yy=Y[d].map({"R":1,"S":0}).dropna(); keep=[r for r in yy.index if r in ridx]
    idx=np.array([ridx[r] for r in keep]); y=yy.loc[keep].astype(int).values
    if y.sum()<40: continue
    Xr=Mden[idx]; Xz=Zall[idx]; Xrz=np.hstack([Xr,Xz])
    tr,te=train_test_split(np.arange(len(y)),test_size=0.30,random_state=7,stratify=y)
    r={}
    for name,XX in [("raw",Xr),("latent",Xz),("raw+latent",Xrz)]:
        lr=LogisticRegression(max_iter=2000,C=1.0).fit(XX[tr],y[tr])
        a,s,v=metrics(y[te],lr.predict_proba(XX[te])[:,1]); r[name]=(a,s,v)
    out.append(dict(drug=d,pctR=round(100*y.mean(),1),
        raw_sens=round(r["raw"][1],3), lat_sens=round(r["latent"][1],3), rl_sens=round(r["raw+latent"][1],3),
        raw_VME=round(r["raw"][2],3), lat_VME=round(r["latent"][2],3), rl_VME=round(r["raw+latent"][2],3),
        raw_AUROC=round(r["raw"][0],3), lat_AUROC=round(r["latent"][0],3)))
res=pd.DataFrame(out); res.to_csv("cryptic/vae_latent_results.csv",index=False)
pd.set_option("display.width",180)
print("\n=== sensitivity: raw vs latent vs raw+latent ===")
print(res[["drug","pctR","raw_sens","lat_sens","rl_sens"]].to_string(index=False))
print("\n=== AUROC: raw vs latent-only ===")
print(res[["drug","raw_AUROC","lat_AUROC"]].to_string(index=False))
