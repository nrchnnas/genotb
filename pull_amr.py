import urllib.request, urllib.parse, json, sys, time
API="https://www.bv-brc.org/api/genome_amr/"
DRUGS=["isoniazid","rifampin","ethambutol","pyrazinamide","streptomycin",
       "ofloxacin","moxifloxacin","amikacin","kanamycin","capreomycin","ciprofloxacin","levofloxacin"]
fields="genome_id,genome_name,antibiotic,resistant_phenotype,measurement_value,measurement_unit"
def page(start,n=25000):
    q=("eq(taxon_id,1773)&in(resistant_phenotype,(Resistant,Susceptible))"
       f"&select({fields})&sort(+id)&limit({n},{start})")
    url=API+"?"+q
    req=urllib.request.Request(url,headers={"Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=40) as r:
        return json.load(r)
rows=[]
start=0
for i in range(9):
    try:
        d=page(start)
    except Exception as e:
        print("err",e); break
    if not d: break
    rows.extend(d); start+=len(d)
    print(f"page {i}: +{len(d)} total {len(rows)}", flush=True)
    if len(d)<25000: break
print("TOTAL rows", len(rows))
json.dump(rows, open("data/tb_amr_raw.json","w"))
