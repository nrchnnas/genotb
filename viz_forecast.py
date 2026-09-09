"""Figures for the emergence/pre-resistance forecaster."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":10,"axes.titleweight":"bold","figure.dpi":130,"savefig.bbox":"tight"})
INK,RES,SUS,MID,ACC="#1b2a4a","#c0392b","#2e86c1","#e08a1e","#6aa84f"

res=pd.read_csv("cryptic/forecast_metrics.csv")
eg=pd.read_csv("cryptic/evolvability_by_class.csv",index_col=0)

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,5.6))

# Panel A: full vs background-only AUROC
d=res.copy(); x=np.arange(len(d)); w=0.36
ax1.bar(x-w/2, d["AUROC_full"], w, color=INK, label="Full model (sees target gene)")
ax1.bar(x+w/2, d["AUROC_background_only"], w, color=MID, label="Background only (target gene removed)")
# cryptic markers
cr=d["cryptic_resistance_AUROC"].str.extract(r"([\d.]+)").astype(float)[0]
ax1.plot(x+w/2, cr, "o", color=RES, ms=7, label="Cryptic-resistance AUROC (catalogue-susceptible only)")
for i,v in enumerate(d["AUROC_background_only"]): ax1.text(i+w/2, v+0.01, f"{v:.2f}", ha="center", fontsize=8)
ax1.set_xticks(x); ax1.set_xticklabels([s[:4].upper() for s in d["drug"]])
ax1.set_ylim(0.5,1.0); ax1.axhline(0.5,color="#999",ls=":",lw=1)
ax1.set_ylabel("AUROC (held-out)"); ax1.legend(fontsize=8,loc="lower left")
ax1.set_title("A · Resistance prediction from background (target gene removed)",loc="left")

# Panel B: evolvability features rise with escalation
order=["pan-S","mono/poly","MDR","pre-XDR","XDR"]
eg=eg.reindex(order)
ax2.plot(order, eg["rpoC_comp"], "-o", color=RES, lw=2.4, label="rpoC compensatory mutation")
ax2.plot(order, eg["efflux_Tap"], "-o", color=SUS, lw=2.4, label="Tap efflux variant")
ax2.plot(order, eg["whiB7_var"], "-o", color=ACC, lw=2.0, label="whiB7 regulator variant")
ax2.set_ylabel("Prevalence in strains"); ax2.set_ylim(0,0.72)
ax2.legend(fontsize=8,loc="upper left")
ax2.set_title("B · Background-variant prevalence by resistance category",loc="left")
ax2.tick_params(axis="x",rotation=20)
plt.savefig("figures/06_forecast.png"); plt.close()
print("saved figures/06_forecast.png")
print("background-only AUROC range:", round(res.AUROC_background_only.min(),2),"-",round(res.AUROC_background_only.max(),2))
