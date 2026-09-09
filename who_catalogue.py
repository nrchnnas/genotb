"""
Parse the WHO catalogue v2 (machine-readable, via the TB-Profiler tbdb database)
into resistance-associated allele sets per drug, expressed in the same token
format emitted by aa_engine. Enables true per-allele scoring rather than
gene/codon-level features.
"""
import csv, re, collections

AA3to1 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
 "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P",
 "Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Ter":"*","Stop":"*"}

# map catalogue drug names -> our 3-letter codes
DRUGMAP = {"isoniazid":"INH","rifampicin":"RIF","rifabutin":"RFB","ethambutol":"EMB",
 "ethionamide":"ETH","levofloxacin":"LEV","moxifloxacin":"MXF","amikacin":"AMI",
 "kanamycin":"KAN"}
LOF_TYPES = {"frameshift_variant","stop_gained","feature_ablation","start_lost",
             "transcript_ablation"}

def _protein_token(gene, mut):
    """p.Ser450Leu -> gene_S450L ; p.Met1? / p.*450 -> LoF marker"""
    m = re.match(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*|\?)", mut)
    if not m: return None
    ref = AA3to1.get(m.group(1)); pos = m.group(2); alt = m.group(3)
    if ref is None: return None
    if alt in ("?","*","Stop","Ter") or alt=="*":
        return f"{gene}_LoF"
    a1 = AA3to1.get(alt)
    if a1 is None: return None
    if a1=="*": return f"{gene}_LoF"
    return f"{gene}_{ref}{pos}{a1}"

def load(path="cryptic/tbdb_mutations.csv"):
    """returns: assoc[drug3] = set(tokens) ; lof_genes[drug3] = set(genes)"""
    assoc = collections.defaultdict(set)
    lof_genes = collections.defaultdict(set)
    nuc = collections.defaultdict(set)   # nucleotide tokens (rrs/promoter) kept raw
    for r in csv.DictReader(open(path)):
        conf = (r.get("confidence") or "")
        # resistance-associated = WHO Group 1/2 ("Assoc w R", incl interim)
        if not conf.startswith("Assoc w R"): continue
        drug = DRUGMAP.get((r.get("drug") or "").lower())
        if not drug: continue
        gene = r["Gene"]; mut = r["Mutation"]; typ = r.get("type","")
        if typ in LOF_TYPES or mut in LOF_TYPES:
            lof_genes[drug].add(gene); continue
        if mut.startswith("p."):
            tok = _protein_token(gene, mut)
            if tok:
                (lof_genes[drug].add(gene) if tok.endswith("_LoF") else assoc[drug].add(tok))
        elif mut.startswith("c.") or mut.startswith("n."):
            nuc[drug].add(f"{gene}:{mut}")
    return assoc, lof_genes, nuc

if __name__ == "__main__":
    a,l,n = load()
    for d in ["RIF","INH","LEV","AMI","EMB"]:
        ex = [t for t in a[d] if "_" in t][:6]
        print(f"{d}: {len(a[d])} protein alleles, LoF genes {sorted(l[d])[:6]} | e.g. {ex}")
    # sanity: canonical mutations present?
    for tok,drug in [("rpoB_S450L","RIF"),("katG_S315T","INH"),("gyrA_D94G","LEV"),("embB_M306V","EMB")]:
        print(f"  {tok} in {drug} assoc set: {tok in a[drug]}")
