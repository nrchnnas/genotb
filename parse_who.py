"""Per-isolate WHO Group-1/2 catalogue hit per drug (true per-allele scoring).
Resumable; writes cryptic/who_features.csv. Run until DONE."""
import os, glob, time
from aa_engine import Engine
from who_catalogue import load

OUT="cryptic/who_features.csv"; BUDGET=38; START=time.time()
assoc, lof_genes, nuc = load()
DRUGS=["INH","RIF","RFB","EMB","ETH","LEV","MXF","AMI","KAN"]
# nucleotide/promoter rules expressed in engine token names
NUC_RULE={"AMI":{"rrs_1401","rrs_1484"},"KAN":{"rrs_1401","rrs_1484","eis_promoter"},
          "INH":{"fabG1_promoter"},"ETH":{"fabG1_promoter"}}
E=Engine("cryptic/H37Rv.fasta")

def who_hits(toks):
    out={}
    for d in DRUGS:
        hit = any(t in assoc[d] for t in toks)
        if not hit:
            hit = any(f"{g}_LoF" in toks for g in lof_genes[d])
        if not hit and d in NUC_RULE:
            hit = bool(NUC_RULE[d] & toks)
        out["who_"+d]=int(hit)
    return out

done=set()
if os.path.exists(OUT):
    for i,l in enumerate(open(OUT)):
        if i: done.add(l.split(",",1)[0])
else:
    open(OUT,"w").write("ENA_RUN,"+",".join("who_"+d for d in DRUGS)+"\n")

vcfs=sorted(glob.glob("cryptic/vcf/*.vcf.gz"))
todo=[p for p in vcfs if os.path.basename(p)[:-7] not in done]
n=0
with open(OUT,"a") as f:
    for p in todo:
        if time.time()-START>BUDGET: break
        run=os.path.basename(p)[:-7]
        try: h=who_hits(E.tokens(p))
        except Exception: continue
        f.write(run+","+",".join(str(h["who_"+d]) for d in DRUGS)+"\n"); n+=1
print(f"processed {n}; total {len(done)+n}/{len(vcfs)}")
print("DONE" if len(done)+n>=len(vcfs) else "MORE")
