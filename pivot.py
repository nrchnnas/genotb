import json, pandas as pd, collections
rows=[json.loads(l) for l in open("data/tb_amr.jsonl")]
df=pd.DataFrame(rows)
# normalize drug names
df['antibiotic']=df['antibiotic'].str.lower().str.strip()
# map rifampin/rifampicin
df['antibiotic']=df['antibiotic'].replace({'rifampicin':'rifampin'})
df['pheno']=df['resistant_phenotype'].map({'Resistant':1,'Susceptible':0})
df=df.dropna(subset=['pheno'])
print("records:",len(df),"unique genomes:",df['genome_id'].nunique())
# pick drugs of interest with enough data
top=df['antibiotic'].value_counts()
print("\nTop drugs by record count:\n", top.head(14))
# collapse duplicate (genome,drug): majority vote
g=df.groupby(['genome_id','antibiotic'])['pheno'].mean().round().reset_index()
mat=g.pivot(index='genome_id',columns='antibiotic',values='pheno')
DRUGS=['isoniazid','rifampin','ethambutol','pyrazinamide','streptomycin',
       'ofloxacin','moxifloxacin','levofloxacin','ciprofloxacin','amikacin','kanamycin','capreomycin']
DRUGS=[d for d in DRUGS if d in mat.columns]
mat=mat[DRUGS]
print("\nMatrix shape:",mat.shape)
print("\nPer-drug: n tested, %resistant:")
for d in DRUGS:
    s=mat[d].dropna(); 
    if len(s): print(f"  {d:14s} n={len(s):6d}  %R={100*s.mean():4.1f}")
# MDR = resistant to isoniazid AND rifampin
both=mat[['isoniazid','rifampin']].dropna()
mdr=((both['isoniazid']==1)&(both['rifampin']==1)).mean()
print(f"\nIsolates with both INH+RIF results: {len(both)}; MDR (R to both): {100*mdr:.1f}%")
mat.to_csv("data/tb_matrix.csv")
print("saved data/tb_matrix.csv")
