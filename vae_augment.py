"""
Generative oversampling with a Bernoulli VAE (the VAE analogue of SMOTE):
for each rare-resistance drug, train a VAE on the RESISTANT training genotypes,
sample synthetic resistant genotypes to balance the training set, train a
logistic classifier, and evaluate on the held-out REAL test set.

Pre-registered comparison: does VAE augmentation beat SMOTE (and the raw
baseline) on held-out sensitivity / very-major error for rare drugs?
"""
import numpy as np, pandas as pd, scipy.sparse as sp, torch, torch.nn as nn
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, confusion_matrix
from imblearn.over_sampling import SMOTE
torch.manual_seed(0); np.random.seed(0)

M=sp.load_npz("cryptic/mut_matrix.npz"); rows=open("cryptic/mut_rows.txt").read().split("\n")
Y=pd.read_csv("cryptic/pheno_full.csv").set_index("ENA_RUN")
ridx={r:i for i,r in enumerate(rows)}; D=M.shape[1]
DRUGS=["AMI","KAN","MXF","ETH","EMB","LEV"]

class VAE(nn.Module):
    def __init__(s,d,h=256,z=32):
        super().__init__()
        s.e1=nn.Linear(d,h); s.mu=nn.Linear(h,z); s.lv=nn.Linear(h,z)
        s.d1=nn.Linear(z,h); s.out=nn.Linear(h,d)
    def enc(s,x): h=torch.relu(s.e1(x)); return s.mu(h),s.lv(h)
    def dec(s,z): return torch.sigmoid(s.out(torch.relu(s.d1(z))))
    def forward(s,x):
        mu,lv=s.enc(x); std=torch.exp(0.5*lv); z=mu+std*torch.randn_like(std)
        return s.dec(z),mu,lv

def train_vae(Xpos,epochs=60,z=32):
    X=torch.tensor(Xpos,dtype=torch.float32)
    m=VAE(X.shape[1],256,z); opt=torch.optim.Adam(m.parameters(),1e-3)
    for ep in range(epochs):
        perm=torch.randperm(len(X))
        for i in range(0,len(X),128):
            b=X[perm[i:i+128]]
            rec,mu,lv=m(b)
            bce=nn.functional.binary_cross_entropy(rec,b,reduction="sum")
            kl=-0.5*torch.sum(1+lv-mu.pow(2)-lv.exp())
            loss=(bce+kl)/len(b); opt.zero_grad(); loss.backward(); opt.step()
    return m

def sample(m,n,z=32):
    with torch.no_grad():
        probs=m.dec(torch.randn(n,z)).numpy()
    return (np.random.rand(n,probs.shape[1])<probs).astype(np.float32)

def metrics(y,p):
    pred=(p>=.5).astype(int); tn,fp,fn,tp=confusion_matrix(y,pred,labels=[0,1]).ravel()
    return roc_auc_score(y,p),(tp/(tp+fn) if tp+fn else np.nan),(fn/(tp+fn) if tp+fn else np.nan)

out=[]
for d in DRUGS:
    yy=Y[d].map({"R":1,"S":0}).dropna()
    keep=[r for r in yy.index if r in ridx]; idx=[ridx[r] for r in keep]
    X=M[idx].toarray().astype(np.float32); y=yy.loc[keep].astype(int).values
    if y.sum()<40: continue
    Xtr,Xte,ytr,yte=train_test_split(X,y,test_size=0.30,random_state=7,stratify=y)
    # VAE trained on minority (resistant) training genotypes
    Xpos=Xtr[ytr==1]; nneg=int((ytr==0).sum()); npos=int((ytr==1).sum())
    m=train_vae(Xpos, epochs=60)
    syn=sample(m, max(nneg-npos,0))                      # balance to majority count
    Xaug=np.vstack([Xtr, syn]); yaug=np.concatenate([ytr, np.ones(len(syn),int)])
    lr=LogisticRegression(max_iter=2000,C=1.0).fit(Xaug,yaug)
    a,s,v=metrics(yte, lr.predict_proba(Xte)[:,1])
    out.append(dict(drug=d,n=len(y),pctR=round(100*y.mean(),1),
        vae_AUROC=round(a,3),vae_sens=round(s,3),vae_VME=round(v,3),
        n_synth=len(syn)))
    print("done",d,"vae_sens %.3f vae_VME %.3f"%(s,v),flush=True)

res=pd.DataFrame(out); res.to_csv("cryptic/vae_results.csv",index=False)
# merge with baselines for a side-by-side
base=pd.read_csv("cryptic/gen_baselines.csv")
comp=base.merge(res[["drug","vae_sens","vae_VME","vae_AUROC"]],on="drug",how="right")
comp.to_csv("cryptic/gen_comparison.csv",index=False)
pd.set_option("display.width",180)
print("\n=== sensitivity: baseline vs SMOTE vs VAE ===")
print(comp[["drug","pctR","base_sens","smote_sens","vae_sens"]].to_string(index=False))
print("\n=== very-major error: baseline vs SMOTE vs VAE (lower better) ===")
print(comp[["drug","base_VME","smote_VME","vae_VME"]].to_string(index=False))
