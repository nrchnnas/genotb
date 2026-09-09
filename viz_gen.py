"""Figure: does generative modeling help rare-variant resistance? (honest negative)"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":10,"axes.titleweight":"bold","figure.dpi":130,"savefig.bbox":"tight"})
INK,SM,VA,VL="#1b2a4a","#1f8a4c","#c0392b","#e08a1e"

b=pd.read_csv("cryptic/gen_baselines.csv")
v=pd.read_csv("cryptic/vae_results.csv")[["drug","vae_sens","vae_VME"]]
l=pd.read_csv("cryptic/vae_latent_results.csv")[["drug","rl_sens","rl_VME"]]
df=b.merge(v,on="drug").merge(l,on="drug")
df=df[df.pctR<25].reset_index(drop=True)   # rare/low-data drugs
x=np.arange(len(df)); w=0.2

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,5.4))
# sensitivity (higher better)
ax1.bar(x-1.5*w, df.base_sens, w, color=INK,  label="Discriminative baseline")
ax1.bar(x-0.5*w, df.smote_sens,w, color=SM,   label="SMOTE oversampling")
ax1.bar(x+0.5*w, df.vae_sens,  w, color=VA,   label="VAE augmentation")
ax1.bar(x+1.5*w, df.rl_sens,   w, color=VL,   label="VAE latent (+raw)")
ax1.set_xticks(x); ax1.set_xticklabels(df.drug); ax1.set_ylim(0,1)
ax1.set_ylabel("Held-out sensitivity"); ax1.legend(fontsize=8,loc="lower right")
ax1.set_title("A · Sensitivity on rare-resistance drugs\nSMOTE wins; generative models do not help",loc="left")
# VME (lower better)
ax2.bar(x-1.5*w, df.base_VME, w, color=INK)
ax2.bar(x-0.5*w, df.smote_VME,w, color=SM)
ax2.bar(x+0.5*w, df.vae_VME,  w, color=VA)
ax2.bar(x+1.5*w, df.rl_VME,   w, color=VL)
ax2.set_xticks(x); ax2.set_xticklabels(df.drug)
ax2.set_ylabel("Very-major error (lower better)")
ax2.set_title("B · Very-major error\nSMOTE lowest across the board",loc="left")
plt.savefig("figures/09_generative_comparison.png"); plt.close()
print("saved figures/09_generative_comparison.png")
print("\nmean sensitivity — baseline %.3f | SMOTE %.3f | VAE-aug %.3f | VAE-latent %.3f"%(
    df.base_sens.mean(), df.smote_sens.mean(), df.vae_sens.mean(), df.rl_sens.mean()))
