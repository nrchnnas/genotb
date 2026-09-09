"""
Turn the amino-acid token set from aa_engine into clean, codon-restricted,
WHO-catalogue-aligned ML features. Restricting to known resistance codons is
what removes the neutral phylogenetic SNPs (e.g. gyrA S95T, present in ~all
isolates) that contaminate a naive 'any nonsynonymous' feature.
"""
import re
_TOK = re.compile(r"^([A-Za-z0-9]+)_([A-Z])(\d+)([A-Z\*])?$")

def _codons(tokens):
    """yield (gene, codon_int, alt_or_None) for amino-acid tokens"""
    for t in tokens:
        m = _TOK.match(t)
        if m:
            yield m.group(1), int(m.group(3)), m.group(4)

# resistance codon windows (WHO catalogue regions)
def feature_vector(tokens):
    cod = list(_codons(tokens))
    def inrange(gene, lo, hi):
        return int(any(g == gene and lo <= c <= hi for g, c, _ in cod))
    def has_codon(gene, codon):
        return int(any(g == gene and c == codon for g, c, _ in cod))
    def specific(gene, codon, alt):
        return int(any(g == gene and c == codon and a == alt for g, c, a in cod))
    flag = lambda name: int(name in tokens)

    f = {}
    # --- rifampin ---
    f["rpoB_S450L"]     = specific("rpoB", 450, "L")
    f["rpoB_450_any"]   = has_codon("rpoB", 450)
    f["rpoB_RRDR"]      = inrange("rpoB", 426, 452)          # whole RRDR
    f["rpoB_435_445"]   = int(has_codon("rpoB",435) or has_codon("rpoB",445))
    # --- isoniazid / ethionamide ---
    f["katG_S315"]      = has_codon("katG", 315)
    f["katG_LoF"]       = flag("katG_LoF")
    f["inhA_promoter"]  = flag("fabG1_promoter")
    f["inhA_I194_S94"]  = int(has_codon("inhA",194) or has_codon("inhA",94))
    f["ethA_LoF"]       = flag("ethA_LoF")
    f["ethA_nonsyn"]    = flag("ethA_nonsyn")
    # --- ethambutol ---
    f["embB_306"]       = has_codon("embB", 306)
    f["embB_406_497"]   = int(has_codon("embB",406) or has_codon("embB",497))
    f["embB_354_328"]   = int(has_codon("embB",354) or has_codon("embB",328))
    # --- fluoroquinolones ---
    f["gyrA_QRDR"]      = inrange("gyrA", 88, 94)
    f["gyrA_D94"]       = has_codon("gyrA", 94)
    f["gyrA_A90"]       = has_codon("gyrA", 90)
    f["gyrB_QRDR"]      = inrange("gyrB", 446, 502)
    # --- aminoglycosides / injectables ---
    f["rrs_1401"]       = flag("rrs_1401")
    f["rrs_1484"]       = flag("rrs_1484")
    f["eis_promoter"]   = flag("eis_promoter")
    # --- streptomycin ---
    f["rpsL_43_88"]     = int(has_codon("rpsL",43) or has_codon("rpsL",88))
    f["gid_LoF"]        = flag("gid_LoF")
    # --- bedaquiline / clofazimine / linezolid ---
    f["Rv0678_LoF"]     = int(flag("Rv0678_LoF") or any(g=="Rv0678" for g,_,_ in cod))
    f["atpE_nonsyn"]    = flag("atpE_nonsyn")
    f["rplC_154"]       = has_codon("rplC", 154)
    # --- evolvability / background (NOT direct drug targets) ---
    f["rpoC_comp"]      = int(any(g=="rpoC" for g,_,_ in cod))   # compensatory; stabilises RIF-R
    f["whiB7_var"]      = int(any(g=="whiB7" for g,_,_ in cod) or flag("whiB7_LoF"))
    f["efflux_Tap"]     = int(any(g=="Rv1258c" for g,_,_ in cod) or flag("Rv1258c_LoF"))
    f["mutator"]        = int(any(g in ("mutT1","nucS") for g,_,_ in cod)
                              or flag("mutT1_LoF") or flag("nucS_LoF"))
    return f

# stable feature order
FEATURES = list(feature_vector(set()).keys())
