"""H37Rv (NC_000962.3, 1-based) coordinate windows for M. tuberculosis
resistance determinants. Used to turn raw VCF variant positions into real,
interpretable genotype features. Whole-gene windows + a few high-resolution
hotspot windows for the canonical codons.
Coordinates from the H37Rv annotation / WHO catalogue loci."""

# gene -> (start, end)  : "any variant in this gene" feature
GENES = {
    "gyrA":   (7302, 9818),
    "gyrB":   (5240, 7267),
    "rpoB":   (759807, 763325),
    "rpsL":   (781560, 781934),
    "rrs":    (1471846, 1473382),
    "fabG1":  (1673440, 1674183),   # inhA promoter region (mabA)
    "inhA":   (1674202, 1675011),
    "rplC":   (800809, 801462),     # linezolid
    "rrl":    (1473658, 1476795),   # 23S rRNA, linezolid
    "atpE":   (1461044, 1461290),   # bedaquiline
    "Rv0678": (778990, 779487),     # mmpR, bedaquiline/clofazimine
    "pncA":   (2288681, 2289241),
    "katG":   (2153889, 2156111),
    "ahpC":   (2726088, 2726780),
    "eis":    (2714124, 2715332),
    "ethA":   (4326004, 4327476),
    "embA":   (4243223, 4246517),
    "embB":   (4246514, 4249810),
    "gid":    (4407528, 4408202),
    "tlyA":   (1917940, 1918746),
}

# high-resolution hotspot windows -> (start, end) : specific codon/region feature
HOTSPOTS = {
    "rpoB_RRDR":      (760872, 761160),   # codons 426-452 (incl. S450)
    "katG_315":       (2155140, 2155200),
    "gyrA_QRDR":      (7565, 7583),        # codons 90-94 ONLY (excludes neutral S95T at 7585)
    "rrs_1401":       (1473236, 1473256),
    "inhA_promoter":  (1674200, 1674260),
    "embB_306":       (4247410, 4247450),
    "eis_promoter":   (2715332, 2715632),
    "rpsL_43_88":     (781670, 781820),
}
