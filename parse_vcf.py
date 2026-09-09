"""Parse real CRyPTIC VCFs -> binary genotype feature matrix (real variant calls)."""
import os, gzip, glob, csv, bisect
from genes_h37rv import GENES, HOTSPOTS

WINDOWS = {**{f"gene_{g}": rng for g, rng in GENES.items()},
           **{f"hot_{h}": rng for h, rng in HOTSPOTS.items()}}
FEATURES = list(WINDOWS)

def variants(path):
    """yield POS (int) for PASS, non-reference calls"""
    op = gzip.open(path, "rt")
    for line in op:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        pos, flt, fmt, smp = int(f[1]), f[6], f[8], f[9]
        if flt != "PASS":
            continue
        gt = smp.split(":")[0]
        # non-reference genotype (e.g. 1/1, 1/2); skip 0/0 and ./.
        alleles = gt.replace("|", "/").split("/")
        if all(a in ("0", ".") for a in alleles):
            continue
        yield pos

vcfs = sorted(glob.glob("cryptic/vcf/*.vcf.gz"))
print("VCFs to parse:", len(vcfs))

# precompute sorted window edges per feature for fast membership
wins = [(name, s, e) for name, (s, e) in WINDOWS.items()]

rows = {}
for i, path in enumerate(vcfs):
    run = os.path.basename(path)[:-7]  # strip .vcf.gz
    feat = {name: 0 for name in FEATURES}
    try:
        for pos in variants(path):
            for name, s, e in wins:
                if s <= pos <= e:
                    feat[name] = 1
    except Exception as ex:
        continue
    rows[run] = feat
    if (i + 1) % 500 == 0:
        print("  parsed", i + 1, flush=True)

# write matrix
with open("cryptic/genotypes_real.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["ENA_RUN"] + FEATURES)
    for run, feat in rows.items():
        w.writerow([run] + [feat[c] for c in FEATURES])
print("samples parsed:", len(rows), "-> cryptic/genotypes_real.csv")

# quick sanity: prevalence of key features
import statistics
tot = len(rows)
for key in ["hot_rpoB_RRDR", "hot_katG_315", "hot_gyrA_QRDR", "hot_rrs_1401", "gene_pncA"]:
    p = sum(r[key] for r in rows.values()) / tot
    print(f"  {key:18s} present in {100*p:4.1f}% of isolates")
