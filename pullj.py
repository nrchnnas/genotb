import urllib.request, json, os
API="https://www.bv-brc.org/api/genome_amr/"
F="data/tb_amr.jsonl"
fields="genome_id,genome_name,antibiotic,resistant_phenotype,measurement_value,measurement_unit"
def page(start,n=25000):
    q=("eq(taxon_id,1773)&in(resistant_phenotype,(Resistant,Susceptible))"
       f"&select({fields})&sort(+id)&limit({n},{start})")
    req=urllib.request.Request(API+"?"+q,headers={"Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=38) as r: return json.load(r)
start=sum(1 for _ in open(F)) if os.path.exists(F) else 0
for i in range(3):
    try: d=page(start)
    except Exception as e: print("err",e); break
    if not d: break
    with open(F,"a") as fh:
        for rec in d: fh.write(json.dumps(rec)+"\n")
    start+=len(d); print("total",start,flush=True)
    if len(d)<25000: break
