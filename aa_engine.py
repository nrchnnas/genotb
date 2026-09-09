"""
Amino-acid-resolution variant engine for M. tuberculosis (H37Rv / AL123456.3).

For each VCF variant it computes the actual codon / amino-acid change (or the
nucleotide change for the rRNA genes), then emits interpretable, WHO-catalogue-
aligned features. SNVs are translated exactly; in-frame/frameshift indels and
stop-gains in loss-of-function genes (pncA, katG, ethA, gidB, Rv0678) are flagged
as LoF, which is how resistance actually arises at those loci.
"""
import gzip

# ---- H37Rv CDS / rRNA coordinates (1-based inclusive) and strand ----
# kind: 'cds' = translate to amino acid; 'rrna' = nucleotide-level
GENES = {
    "gyrB":   (5240, 7267, "+", "cds"),
    "gyrA":   (7302, 9818, "+", "cds"),
    "rpoB":   (759807, 763325, "+", "cds"),
    "rpsL":   (781560, 781934, "+", "cds"),
    "Rv0678": (778990, 779487, "+", "cds"),
    "rplC":   (800809, 801462, "+", "cds"),
    "atpE":   (1461044, 1461290, "+", "cds"),
    "rrs":    (1471846, 1473382, "+", "rrna"),
    "rrl":    (1473658, 1476795, "+", "rrna"),
    "fabG1":  (1673440, 1674183, "+", "cds"),
    "inhA":   (1674202, 1675011, "+", "cds"),
    "pncA":   (2288681, 2289241, "-", "cds"),
    "katG":   (2153889, 2156111, "-", "cds"),
    "ahpC":   (2726088, 2726780, "+", "cds"),
    "eis":    (2714124, 2715332, "-", "cds"),
    "ethA":   (4326004, 4327473, "-", "cds"),
    "embA":   (4243233, 4246517, "+", "cds"),
    "embB":   (4246514, 4249810, "+", "cds"),
    "embC":   (4239863, 4243147, "+", "cds"),
    "gid":    (4407528, 4408202, "-", "cds"),
    # ---- evolvability / "pre-resistance" background genes ----
    "rpoC":   (763370, 767320, "+", "cds"),   # compensatory mutations stabilise RIF resistance
    "whiB7":  (3568402, 3568779, "-", "cds"), # intrinsic-resistance / efflux regulator
    "Rv1258c":(1406492, 1407985, "-", "cds"), # Tap efflux pump
    "mutT1":  (3334214, 3334741, "+", "cds"), # DNA-repair / mutator
    "nucS":   (1602876, 1603577, "+", "cds"), # mismatch repair (mutator when lost)
}
# promoter windows (upstream non-coding) -> nucleotide-level features
PROMOTERS = {
    "fabG1_promoter": (1673330, 1673440),   # inhA/mabA promoter (-15, -8 etc.)
    "eis_promoter":   (2715332, 2715432),    # upstream of eis (minus strand)
    "ahpC_promoter":  (2726088, 2726200),
}

COMP = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
def rc(s): return "".join(COMP.get(b,"N") for b in reversed(s))
CODON = {}
_b="TCAG"
_aa="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i,a in enumerate(_b):
    for j,b2 in enumerate(_b):
        for k,c in enumerate(_b):
            CODON[a+b2+c]=_aa[i*16+j*4+k]

class Engine:
    def __init__(self, fasta):
        self.seq = "".join(l.strip() for l in open(fasta) if not l.startswith(">"))
        # gene lookup by coordinate (sorted)
        self.items = sorted(((s,e,name,strand,kind) for name,(s,e,strand,kind) in GENES.items()))

    def base(self, pos):  # 1-based
        return self.seq[pos-1]

    def gene_at(self, pos):
        for s,e,name,strand,kind in self.items:
            if s <= pos <= e:
                return (name,s,e,strand,kind)
        return None

    def aa_change(self, name, s, e, strand, pos, ref, alt):
        """return (aa_ref, aa_num, aa_alt) for a SNV, or None"""
        if strand == "+":
            idx = pos - s
            cstart = s + (idx//3)*3
            codon = self.seq[cstart-1:cstart+2]
            off = idx % 3
            mut = list(codon); mut[off] = alt; mut="".join(mut)
            aa_num = idx//3 + 1
        else:
            d = e - pos
            cstart = e - (d//3)*3 - 2      # lowest genomic coord of codon
            gcodon = self.seq[cstart-1:cstart+2]
            codon = rc(gcodon)
            off = d % 3                      # offset on coding strand
            altc = COMP.get(alt,"N")
            mut = list(codon); mut[off] = altc; mut="".join(mut)
            aa_num = d//3 + 1
        if len(codon)!=3 or len(mut)!=3: return None
        return CODON.get(codon,"X"), aa_num, CODON.get(mut,"X")

    # CRyPTIC FILTER values we accept: high-confidence + minor-allele (heteroresistance)
    ACCEPT_FILTER = {"PASS", "MIN_FRS", "MIN_GCP", "MIN_DP", "."}

    def _snv(self, pos, refb, altb, feats):
        """process a single-base substitution at pos and add tokens"""
        if altb not in "ACGT" or refb == altb:
            return
        for pname,(ps,pe) in PROMOTERS.items():
            if ps<=pos<=pe: feats.add(pname)
        g=self.gene_at(pos)
        if not g: return
        name,s,e,strand,kind=g
        if kind=="rrna":
            rel=pos-s+1
            feats.add(f"{name}_n{rel}")
            if name=="rrs" and abs(pos-1473246)<=2: feats.add("rrs_1401")
            if name=="rrs" and abs(pos-1473329)<=2: feats.add("rrs_1484")
            return
        ch=self.aa_change(name,s,e,strand,pos,refb,altb)
        if not ch: return
        aref,anum,aalt=ch
        if aref==aalt: return                       # synonymous
        feats.add(f"{name}_{aref}{anum}{aalt}")
        feats.add(f"{name}_nonsyn")
        feats.add(f"{name}_{aref}{anum}")           # codon-level
        if aalt=="*": feats.add(f"{name}_LoF"); feats.add(f"{name}_stop")

    def tokens(self, path):
        """parse a VCF.gz -> set of feature tokens for this isolate"""
        feats=set()
        try: op=gzip.open(path,"rt")
        except Exception: return feats
        for line in op:
            if line[0]=="#": continue
            f=line.rstrip("\n").split("\t")
            if len(f)<10 or f[6] not in self.ACCEPT_FILTER: continue
            pos=int(f[1]); ref=f[3]; alt=f[4].split(",")[0]
            gt=f[9].split(":")[0].replace("|","/").split("/")
            if all(a in ("0",".") for a in gt): continue
            if len(ref)==len(alt):
                # SNV or multi-nucleotide substitution: decompose per position
                for i in range(len(ref)):
                    if ref[i]!=alt[i]:
                        self._snv(pos+i, ref[i], alt[i], feats)
            else:
                # indel -> loss of function in CDS, position-level for promoters/rRNA
                for pname,(ps,pe) in PROMOTERS.items():
                    if ps<=pos<=pe: feats.add(pname)
                g=self.gene_at(pos)
                if not g: continue
                name,s,e,strand,kind=g
                if kind=="rrna":
                    feats.add(f"{name}_n{pos-s+1}")
                else:
                    feats.add(f"{name}_LoF"); feats.add(f"{name}_indel")
        return feats
