"""Pass 1: emit each isolate's full set of specific amino-acid mutation tokens
to a JSONL (resumable). Used to build the high-dimensional mutation matrix for
the generative (VAE) rare-variant experiment."""
import os, glob, time, json
from aa_engine import Engine

OUT="cryptic/mut_tokens.jsonl"; BUDGET=40; START=time.time()
E=Engine("cryptic/H37Rv.fasta")

def spec_tokens(toks):
    # keep specific amino-acid changes (gene_<AA><pos><AA>) and LoF flags
    out=[]
    for t in toks:
        if any(ch.isdigit() for ch in t) and t.count("_")==1:
            out.append(t)
        elif t.endswith("_LoF"):
            out.append(t)
    return out

done=set()
if os.path.exists(OUT):
    for l in open(OUT):
        try: done.add(l.split("\t",1)[0])
        except: pass

vcfs=sorted(glob.glob("cryptic/vcf/*.vcf.gz"))
todo=[p for p in vcfs if os.path.basename(p)[:-7] not in done]
n=0
with open(OUT,"a") as f:
    for p in todo:
        if time.time()-START>BUDGET: break
        run=os.path.basename(p)[:-7]
        try: t=spec_tokens(E.tokens(p))
        except Exception: continue
        f.write(run+"\t"+",".join(sorted(set(t)))+"\n"); n+=1
tot=len(done)+n
print(f"processed {n}; total {tot}/{len(vcfs)}")
print("DONE" if tot>=len(vcfs) else "MORE")
