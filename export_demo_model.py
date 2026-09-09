import json, math, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
G=pd.read_csv("cryptic/genotypes_real.csv").set_index("ENA_RUN")
Y=pd.read_csv("cryptic/pheno_subset.csv").set_index("ENA_RUN")
df=G.join(Y,how="inner")
NAME={"INH":"isoniazid","RIF":"rifampin","EMB":"ethambutol","ETH":"ethionamide",
      "LEV":"levofloxacin","MXF":"moxifloxacin","AMI":"amikacin","KAN":"kanamycin"}
FEATS=[("rpoB RRDR (e.g. S450L)","hot_rpoB_RRDR"),("katG codon 315","hot_katG_315"),
       ("inhA promoter","hot_inhA_promoter"),("fabG1 (inhA promoter)","gene_fabG1"),
       ("embB M306","hot_embB_306"),("embB (other)","gene_embB"),
       ("gyrA QRDR (90/94)","hot_gyrA_QRDR"),("gyrB","gene_gyrB"),
       ("rrs a1401g","hot_rrs_1401"),("eis promoter","hot_eis_promoter"),
       ("ethA","gene_ethA"),("rpsL","gene_rpsL")]
cols=[c for _,c in FEATS]
models={}; prev={}
for d3,dl in NAME.items():
    y=df[d3].map({"R":1,"S":0}); keep=y.notna()
    yi=y[keep].astype(int); Xi=df.loc[keep,cols].values
    if yi.nunique()<2: continue
    lr=LogisticRegression(max_iter=1000,C=1.0).fit(Xi,yi)
    models[dl]={"coef":[round(float(c),3) for c in lr.coef_[0]],
                "intercept":round(float(lr.intercept_[0]),3)}
    prev[dl]=round(float(yi.mean()),3)
json.dump({"features":[f for f,_ in FEATS],"drugs":models,"prevalence":prev,
  "note":"Real logistic-regression coefficients trained on CRyPTIC genotype(VCF)->phenotype pairs."},
  open("demo_model.json","w"),indent=1)
def p(d,on):
    m=models[d]; z=m["intercept"]+sum(c*(1 if i in on else 0) for i,c in enumerate(m["coef"]))
    return 1/(1+math.exp(-z))
print("RIF none %.2f  rpoB %.2f"%(p("rifampin",set()),p("rifampin",{0})))
print("INH none %.2f  katG %.2f"%(p("isoniazid",set()),p("isoniazid",{1})))
print("LEV none %.2f  gyrA %.2f"%(p("levofloxacin",set()),p("levofloxacin",{6})))
print("AMI none %.2f  rrs  %.2f"%(p("amikacin",set()),p("amikacin",{8})))
print("OK")
