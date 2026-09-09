"""
MIC-value regression (the technical-novelty layer): predict the actual minimum
inhibitory concentration from the genome, not just resistant/susceptible.

MICs are on a doubling-dilution ladder, so we regress log2(MIC) and report the
metric clinical microbiology uses for quantitative AST: ESSENTIAL AGREEMENT (EA)
= predicted within +-1 doubling dilution of the true MIC. Most deployed tools and
student projects stop at binary R/S; predicting the number is harder and tells
you HOW resistant (relevant to dosing, e.g. high-dose isoniazid).
"""
import re, numpy as np, pandas as pd, csv
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import r2_score

G = pd.read_csv("cryptic/aa_features_v3.csv").set_index("ENA_RUN")
reuse = pd.read_csv("cryptic/reuse.csv").set_index("ENA_RUN")
DRUGS = {"INH":"isoniazid","RIF":"rifampin","EMB":"ethambutol","LEV":"levofloxacin",
         "MXF":"moxifloxacin","AMI":"amikacin","KAN":"kanamycin","ETH":"ethionamide","RFB":"rifabutin"}

def parse_mic(v):
    if not isinstance(v,str) or v in ("","NA"): return np.nan
    s = v.strip().replace("<=","").replace(">=","").replace("<","").replace(">","")
    try: x = float(s)
    except: return np.nan
    return np.log2(x) if x>0 else np.nan

rows=[]
for code,name in DRUGS.items():
    col=code+"_MIC"
    if col not in reuse.columns: continue
    y=reuse[col].map(parse_mic)
    d=G.join(y.rename("mic"),how="inner").dropna(subset=["mic"])
    if len(d)<300: continue
    X=d[G.columns]; yy=d["mic"]
    Xtr,Xte,ytr,yte=train_test_split(X,yy,test_size=0.30,random_state=7)
    reg=GradientBoostingRegressor(n_estimators=300,max_depth=3,learning_rate=0.1,random_state=7).fit(Xtr,ytr)
    pred=reg.predict(Xte)
    err=np.abs(pred-yte)
    ea=float((err<=1.0).mean())          # within +-1 doubling dilution
    ea2=float((err<=2.0).mean())
    mae=float(err.mean())
    r2=r2_score(yte,pred)
    rows.append(dict(drug=name,n=len(d),essential_agreement=round(ea,3),
                     within_2dil=round(ea2,3),MAE_dilutions=round(mae,2),R2=round(r2,3)))

res=pd.DataFrame(rows)
res.to_csv("cryptic/mic_metrics.csv",index=False)
pd.set_option("display.width",140)
print("=== MIC-value regression (genome -> MIC), held-out 30% ===")
print(res.to_string(index=False))
print("\nMean essential agreement (within +-1 dilution): %.1f%%" % (100*res.essential_agreement.mean()))
print("Mean R2: %.3f" % res.R2.mean())
