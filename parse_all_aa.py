"""Parse ALL CRyPTIC VCFs through the AA engine -> codon-restricted feature matrix.
Resumable: appends to cryptic/aa_features_v2.csv, skips runs already done. Self-limits
to ~38s per invocation so it fits the sandbox; call repeatedly until 'DONE'."""
import os, glob, csv, time
from aa_engine import Engine
from catalogue import feature_vector, FEATURES

OUT = "cryptic/aa_features_v2.csv"
START = time.time(); BUDGET = 38

E = Engine("cryptic/H37Rv.fasta")
done = set()
if os.path.exists(OUT):
    with open(OUT) as f:
        for i, line in enumerate(f):
            if i == 0: continue
            done.add(line.split(",", 1)[0])
else:
    with open(OUT, "w") as f:
        f.write("ENA_RUN," + ",".join(FEATURES) + "\n")

vcfs = sorted(glob.glob("cryptic/vcf/*.vcf.gz"))
todo = [p for p in vcfs if os.path.basename(p)[:-7] not in done]
n = 0
with open(OUT, "a") as f:
    for p in todo:
        if time.time() - START > BUDGET:
            break
        run = os.path.basename(p)[:-7]
        try:
            toks = E.tokens(p)
            fv = feature_vector(toks)
        except Exception:
            continue
        f.write(run + "," + ",".join(str(fv[k]) for k in FEATURES) + "\n")
        n += 1
print(f"processed {n} this run; total done {len(done)+n} / {len(vcfs)}")
print("DONE" if len(done) + n >= len(vcfs) else "MORE")
