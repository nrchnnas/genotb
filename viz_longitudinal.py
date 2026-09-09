"""Within-host longitudinal figure: ROC + comparison to cross-sectional."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
plt.rcParams.update({"font.size":10,"axes.titleweight":"bold","figure.dpi":130,"savefig.bbox":"tight"})
INK,RES,MID="#1b2a4a","#c0392b","#e08a1e"

m=pd.read_csv("cryptic/longitudinal_metrics.csv").iloc[0]
yo=np.load("cryptic/_oof.npy"); y=yo[:,0]; oof=yo[:,1]
fpr,tpr,_=roc_curve(y,oof)
fc=pd.read_csv("cryptic/forecast_metrics.csv")

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5.2))
# A: ROC
ax1.plot(fpr,tpr,color=INK,lw=2.2,
         label=f"within-host model (AUROC {m.auc:.2f}, 95% CI {m.lo:.2f}–{m.hi:.2f})")
ax1.plot([0,1],[0,1],color="#999",ls="--",lw=1,label="chance")
ax1.set_xlabel("False-positive rate"); ax1.set_ylabel("True-positive rate")
ax1.legend(fontsize=8,loc="lower right")
ax1.set_title(f"A · Predicting within-host resistance acquisition\n{int(m.pos)} events across {int(m.patients)} patients with serial isolates",loc="left")

# B: cross-sectional vs within-host (honest comparison)
labels=["Cross-sectional\nbackground\n(co-occurring)","Within-host\nprospective\n(acquisition)"]
vals=[fc.AUROC_background_only.mean(), m.auc]
errs=[[0,m.auc-m.lo],[0,m.hi-m.auc]]
x=np.arange(2)
bars=ax2.bar(x,vals,color=[MID,INK],width=0.55)
ax2.errorbar(1,m.auc,yerr=[[m.auc-m.lo],[m.hi-m.auc]],color="k",capsize=5,fmt="none")
for i,v in enumerate(vals): ax2.text(i,v+0.01,f"{v:.2f}",ha="center",fontsize=10)
ax2.axhline(0.5,color="#999",ls=":",lw=1)
ax2.set_xticks(x); ax2.set_xticklabels(labels,fontsize=9)
ax2.set_ylim(0.4,1.0); ax2.set_ylabel("AUROC")
ax2.set_title("B · Anticipation is genuine but modest\nprospective < cross-sectional, as expected",loc="left")
plt.savefig("figures/08_longitudinal.png"); plt.close()
print(f"saved figures/08_longitudinal.png  AUROC {m.auc:.3f} ({m.lo:.3f}-{m.hi:.3f})")
