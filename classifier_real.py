"""
REAL genome->phenotype classifier for M. tuberculosis (CRyPTIC).

Both halves are now real and independent:
  - Features: presence of resistance-gene variants parsed from each isolate's
    own VCF of real variant calls (cryptic/genotypes_real.csv).
  - Labels:   real CRyPTIC binary phenotypes (cryptic/pheno_subset.csv).

No simulation, no conditioning of genotype on phenotype. The model is trained on
70% of isolates and evaluated on a held-out 30% it never saw, so performance is a
genuine generalization estimate and the recovered feature importances are real.
"""
import json, numpy as np, pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, confusion_matrix

G = pd.read_csv("cryptic/genotypes_real.csv").set_index("ENA_RUN")
Y = pd.read_csv("cryptic/pheno_subset.csv").set_index("ENA_RUN")
df = G.join(Y, how="inner")
FEATURES = list(G.columns)
DRUGNAME = {"INH":"isoniazid","RIF":"rifampin","EMB":"ethambutol","ETH":"ethionamide",
            "LEV":"levofloxacin","MXF":"moxifloxacin","AMI":"amikacin","KAN":"kanamycin"}
DRUGS = list(DRUGNAME)

results, importances, full_imp = {}, {}, {}
for d in DRUGS:
    y = df[d].map({"R":1,"S":0})
    keep = y.notna()
    Xi, yi = G.loc[keep.index[keep]], y[keep].astype(int)
    if yi.nunique() < 2 or yi.sum() < 30:
        continue
    Xtr, Xte, ytr, yte = train_test_split(Xi, yi, test_size=0.30,
                                          random_state=7, stratify=yi)
    clf = GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                     learning_rate=0.1, random_state=7)
    clf.fit(Xtr, ytr)
    prob = clf.predict_proba(Xte)[:, 1]
    pred = (prob >= 0.5).astype(int)
    auc = roc_auc_score(yte, prob)
    tn, fp, fn, tp = confusion_matrix(yte, pred, labels=[0, 1]).ravel()
    sens = tp/(tp+fn) if (tp+fn) else float("nan")
    spec = tn/(tn+fp) if (tn+fp) else float("nan")
    vme = fn/(tp+fn) if (tp+fn) else float("nan")
    me = fp/(tn+fp) if (tn+fp) else float("nan")
    results[DRUGNAME[d]] = dict(n=int(keep.sum()), n_test=int(len(yte)),
        pct_R=round(100*yi.mean(),1), AUROC=round(auc,3),
        sensitivity=round(sens,3), specificity=round(spec,3),
        very_major_error=round(vme,3), major_error=round(me,3))
    imp = pd.Series(clf.feature_importances_, index=FEATURES).sort_values(ascending=False)
    importances[DRUGNAME[d]] = imp[imp > 0].head(5).round(3).to_dict()
    full_imp[DRUGNAME[d]] = pd.Series(clf.feature_importances_, index=FEATURES)

res = pd.DataFrame(results).T[["n","n_test","pct_R","AUROC","sensitivity",
                               "specificity","very_major_error","major_error"]]
res.to_csv("cryptic/model_metrics_real.csv")
json.dump(importances, open("cryptic/feature_importance_real.json","w"), indent=2)
pd.DataFrame(full_imp).to_csv("cryptic/importance_matrix_real.csv")

pd.set_option("display.width",140)
print("=== REAL genotype->phenotype performance (CRyPTIC, held-out 30%) ===")
print(f"Isolates with parsed VCF + phenotype: {len(df)}")
print(res.to_string())
print("\nMacro-AUROC %.3f | mean sensitivity %.3f | mean specificity %.3f"
      % (res.AUROC.mean(), res.sensitivity.mean(), res.specificity.mean()))
print("\n=== Top real features recovered per drug ===")
for d, imp in importances.items():
    print(f"  {d:13s} -> {list(imp)[:3]}")
