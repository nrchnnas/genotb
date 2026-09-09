import urllib.request, json, os
API="https://www.bv-brc.org/api/genome_amr/"
F="data/tb_amr_raw.json"
fields="genome_id,genome_name,antibiotic,resistant_phenotype,measurement_value,measurement_unit"
def page(start,n=25000):
    q=("eq(taxon_id,1773)&in(resistant_phenotype,(Resistant,Susceptible))"
       f"&select({fields})&sort(+id)&limit({n},{start})")
    req=urllib.request.Request(API+"?"+q,headers={"Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=38) as r: return json.load(r)
rows=json.load(open(F)) if os.path.exists(F) else []
start=len(rows)
for i in range(3):
    try: d=page(start)
    except Exception as e: print("err",e); break
    if not d: break
    rows.extend(d); start+=len(d)
    json.dump(rows,open(F,"w"))      # save EACH page
    print("saved total",len(rows),flush=True)
    if len(d)<25000: break
