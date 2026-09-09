import urllib.request, json
API="https://www.bv-brc.org/api/genome_amr/"
fields="genome_id,genome_name,antibiotic,resistant_phenotype,measurement_value,measurement_unit"
def page(start,n=25000):
    q=("eq(taxon_id,1773)&in(resistant_phenotype,(Resistant,Susceptible))"
       f"&select({fields})&sort(+id)&limit({n},{start})")
    req=urllib.request.Request(API+"?"+q,headers={"Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=40) as r: return json.load(r)
rows=json.load(open("data/tb_amr_raw.json"))
start=len(rows)
for i in range(4):
    try: d=page(start)
    except Exception as e: print("err",e); break
    if not d: break
    rows.extend(d); start+=len(d)
    print("now",len(rows),flush=True)
    if len(d)<25000: break
json.dump(rows, open("data/tb_amr_raw.json","w"))
print("saved",len(rows))
