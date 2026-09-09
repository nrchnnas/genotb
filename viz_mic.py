"""MIC-regression figure: essential agreement per drug + predicted-vs-true ladder."""
import re, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
plt.rcParams.update({"font.size":10,"axes.titleweight":"bold","figure.dpi":130,"savefig.bbox":"tight"})
INK,RES,SUS,MID="#1b2a4a","#c0392b","#2e86c1","#e08a1e"

met=pd.read_csv("cryptic/mic_metrics.csv")
G=pd.read_csv("cryptic/aa_features_v3.csv").set_index("ENA_RUN")
reuse=pd.read_csv("cryptic/reuse.csv").set_index("ENA_RUN")
def parse(v):
    if not isinstance(v,str) or v in ("","NA"): return np.nan
    s=v.replace("<=","").replace(">=","").replace("<","").replace(">","")
    try:x=float(s)
    except:return np.nan
    return np.log2(x) if x>0 else np.nan

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,5.4))
# Panel A
d=met.sort_values("essential_agreement"); y=np.arange(len(d))
ax1.barh(y,d["within_2dil"],color="#cfe0f0",label="within ±2 dilutions")
ax1.barh(y,d["essential_agreement"],color=INK,label="essential agreement (±1 dilution)")
for i,v in enumerate(d["essential_agreement"]): ax1.text(v+0.01,i,f"{v:.0%}",va="center",fontsize=8)
ax1.set_yticks(y); ax1.set_yticklabels([s.capitalize() for s in d["drug"]])
ax1.set_xlim(0,1.05); ax1.set_xlabel("Agreement with measured MIC (held-out)")
ax1.legend(fontsize=8,loc="lower right")
ax1.set_title("A · MIC prediction accuracy (held-out)",loc="left")

# Panel B: predicted vs true MIC ladder for isoniazid
y=reuse["INH_MIC"].map(parse); dd=G.join(y.rename("mic"),how="inner").dropna(subset=["mic"])
X=dd[G.columns]; yy=dd["mic"]
Xtr,Xte,ytr,yte=train_test_split(X,yy,test_size=0.30,random_state=7)
reg=GradientBoostingRegressor(n_estimators=300,max_depth=3,random_state=7).fit(Xtr,ytr)
pred=reg.predict(Xte)
hb=ax2.hexbin(yte,pred,gridsize=22,cmap="Blues",mincnt=1)
lim=[min(yte.min(),pred.min())-0.5,max(yte.max(),pred.max())+0.5]
ax2.plot(lim,lim,color=RES,lw=1.5,ls="--",label="perfect")
ax2.plot(lim,[v+1 for v in lim],color=MID,lw=1,ls=":")
ax2.plot(lim,[v-1 for v in lim],color=MID,lw=1,ls=":",label="±1 dilution")
ax2.set_xlabel("Measured MIC (log2)"); ax2.set_ylabel("Predicted MIC (log2)")
ax2.legend(fontsize=8,loc="upper left")
ax2.set_title("B · Isoniazid: predicted vs measured MIC",loc="left")
fig.colorbar(hb,ax=ax2,fraction=0.046,pad=0.04,label="isolates")
plt.savefig("figures/07_mic.png"); plt.close()
print("saved figures/07_mic.png; mean EA %.1f%%"%(100*met.essential_agreement.mean()))
