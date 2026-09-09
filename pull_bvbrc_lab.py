"""Re-pull BV-BRC AMR records for M. tuberculosis restricted to evidence =
'Laboratory Method'. The original pull did not filter on evidence, and 92% of
BV-BRC R/S records for taxon 1773 are 'Computational Method' (XGBoost-predicted
MICs), not measurements. Resumable JSONL."""
import urllib.request, json, os

API="https://www.bv-brc.org/api/genome_amr/"
OUT="data/tb_amr_lab.jsonl"
FIELDS="genome_id,genome_name,antibiotic,resistant_phenotype,evidence,laboratory_typing_method"
UA={"Accept":"application/json","User-Agent":"Mozilla/5.0 (research; reproducibility check)"}

def page(start,n=25000):
    q=('eq(taxon_id,1773)&in(resistant_phenotype,(Resistant,Susceptible))'
       '&eq(evidence,%22Laboratory%20Method%22)'
       f'&select({FIELDS})&sort(+id)&limit({n},{start})')
    r=urllib.request.Request(API+"?"+q,headers=UA)
    return json.load(urllib.request.urlopen(r,timeout=45))

start=sum(1 for _ in open(OUT)) if os.path.exists(OUT) else 0
for _ in range(6):
    try: d=page(start)
    except Exception as e: print("err",e); break
    if not d: break
    with open(OUT,"a") as f:
        for rec in d: f.write(json.dumps(rec)+"\n")
    start+=len(d); print("total",start,flush=True)
    if len(d)<25000: break
print("DONE" if start>=100000 else "MORE")
