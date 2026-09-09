const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
        ImageRun, Header, Footer, AlignmentType, LevelFormat, HeadingLevel,
        TableOfContents, BorderStyle, WidthType, ShadingType, PageNumber,
        PageBreak } = require("docx");

const INK = "1b2a4a", ACCENT = "2e3a52", RES = "c0392b";
const CW = 9360;

// ---- load metrics ----
function loadCsv(path){
  const L = fs.readFileSync(path,"utf8").trim().split("\n");
  const head = L[0].split(",");
  return L.slice(1).map(r=>{const c=r.split(","); const o={}; head.forEach((h,i)=>o[h||"drug"]=c[i]); return o;});
}
const clf = loadCsv("cryptic/model_metrics_full.csv");        // key "" -> drug
const clfRows = clf.map(r=>({drug:r.drug, n:+r.n, pctR:+r.pct_R, auc:+r.AUROC,
  sens:+r.sensitivity, spec:+r.specificity, vme:+r.very_major_error, me:+r.major_error}));
const mic = loadCsv("cryptic/mic_metrics.csv");
const macroAUC=(clfRows.reduce((s,d)=>s+d.auc,0)/clfRows.length).toFixed(3);
const meanSpec=(clfRows.reduce((s,d)=>s+d.spec,0)/clfRows.length);
const meanSens=(clfRows.reduce((s,d)=>s+d.sens,0)/clfRows.length);
const micEA=(mic.reduce((s,d)=>s+ +d.essential_agreement,0)/mic.length);

// ---- helpers ----
const H = (t) => new Paragraph({ heading: HeadingLevel.HEADING_1, children:[new TextRun(t)] });
const H2 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_2, children:[new TextRun(t)] });
const P = (runs) => new Paragraph({ spacing:{after:120, line:300}, alignment:AlignmentType.JUSTIFIED,
  children:(Array.isArray(runs)?runs:[new TextRun(runs)]) });
const t  = (s)=> new TextRun(s);
const it = (s)=> new TextRun({text:s, italics:true});         // gene names / species
const bl = (s)=> new TextRun({text:s, bold:true});            // structured-abstract labels only
function img(path,w,h){ return new Paragraph({ alignment:AlignmentType.CENTER, spacing:{before:120,after:60},
  children:[ new ImageRun({ type:"png", data:fs.readFileSync(path), transformation:{width:w,height:h},
    altText:{title:"figure",description:"figure",name:"figure"} }) ]}); }
function cap(n, rest){ return new Paragraph({ spacing:{after:200}, alignment:AlignmentType.JUSTIFIED,
  children:[ new TextRun({text:"Figure "+n+". ", bold:true, size:18}),
             new TextRun({text:rest, size:18}) ]}); }

function table(headRow, dataRows, widths){
  const cell=(txt,o={})=>new TableCell({ width:{size:o.w,type:WidthType.DXA},
    shading:o.fill?{fill:o.fill,type:ShadingType.CLEAR}:undefined,
    margins:{top:50,bottom:50,left:100,right:100},
    borders:{top:{style:BorderStyle.SINGLE,size:1,color:"BBBBBB"},bottom:{style:BorderStyle.SINGLE,size:1,color:"BBBBBB"},
             left:{style:BorderStyle.NONE},right:{style:BorderStyle.NONE}},
    children:[new Paragraph({alignment:o.align||AlignmentType.LEFT,
      children:[new TextRun({text:txt,bold:!!o.bold,italics:!!o.italic,color:o.color,size:18})]})]});
  const headCells = headRow.map((h,i)=>cell(h,{w:widths[i],bold:true,
    align:i?AlignmentType.CENTER:AlignmentType.LEFT,fill:"e9edf3"}));
  const rows=[ new TableRow({tableHeader:true, children:headCells}) ];
  dataRows.forEach(r=>{
    const cells=r.map((c,i)=>{
      const obj = (c!==null && typeof c==="object");
      return cell(obj?c.v:String(c), {w:widths[i],
        align:i?AlignmentType.CENTER:AlignmentType.LEFT,
        italic: obj && c.italic, color: obj && c.color ? c.color : undefined});
    });
    rows.push(new TableRow({children:cells}));
  });
  return new Table({width:{size:CW,type:WidthType.DXA},columnWidths:widths,rows});
}
function tableCap(n,rest){ return new Paragraph({ spacing:{before:60,after:160},
  children:[ new TextRun({text:"Table "+n+". ", bold:true, size:18}), new TextRun({text:rest,size:18}) ]}); }

const cap1=s=>s.charAt(0).toUpperCase()+s.slice(1);
const children=[];

// ===================== TITLE =====================
children.push(
  new Paragraph({ spacing:{before:1200, after:160}, alignment:AlignmentType.CENTER, children:[
    new TextRun({text:"Genome-based prediction of drug susceptibility, minimum inhibitory concentration, and pre-resistance risk in ",bold:true,size:34,color:INK}),
    new TextRun({text:"Mycobacterium tuberculosis",bold:true,italics:true,size:34,color:INK}),
    new TextRun({text:": a reproducible analysis of 12,245 clinical isolates",bold:true,size:34,color:INK})]}),
  new Paragraph({ alignment:AlignmentType.CENTER, spacing:{after:360},
    children:[new TextRun({text:"Open data (CRyPTIC, BV-BRC) · all code and analyses reproducible",italics:true,size:20,color:"55607a"})]}),
);

// ===================== ABSTRACT =====================
children.push(new Paragraph({heading:HeadingLevel.HEADING_1,children:[new TextRun("Abstract")]}));
children.push(P([
  bl("Background. "),
  t("Phenotypic drug-susceptibility testing (DST) of "), it("Mycobacterium tuberculosis"),
  t(" requires weeks of culture, delaying effective therapy; in Thailand the reported median interval from specimen to DST result is 114 days. Whole-genome sequencing can infer resistance within hours, and catalogue-based tools are already deployed, but most operate on binary susceptibility and do not quantify the risk that a currently-treatable isolate will acquire resistance. "),
  bl("Methods. "),
  t("Using the CRyPTIC compendium (n = 12,245 isolates with matched genome and 96-well broth-microdilution phenotypes), we translated each isolate’s variant calls to amino-acid substitutions, restricted features to World Health Organization (WHO) catalogue resistance loci, and trained, with 30% held-out evaluation, (i) per-drug susceptibility classifiers, (ii) regression models for minimum inhibitory concentration (MIC), and (iii) a pre-resistance analysis in which a drug’s own target gene was withheld. "),
  bl("Results. "),
  t("Susceptibility prediction across nine drugs reached a mean area under the receiver-operating curve (AUROC) of "),
  t(macroAUC), t(" (specificity "), t((meanSpec*100).toFixed(0)+"%"),
  t("). MIC was predicted to within one doubling dilution (essential agreement) in "),
  t((micEA*100).toFixed(0)+"%"),
  t(" of held-out isolates. With the target gene removed, genomic background predicted second-line resistance at AUROC 0.83–0.94, but a specificity analysis attributes most of this to multidrug-resistant status rather than to pre-resistance: a single MDR indicator alone reached AUROC 0.813 versus 0.894 for the full background model, and discrimination fell to 0.725 within non-MDR isolates. Background features nonetheless retained discrimination among catalogue-susceptible isolates (moxifloxacin AUROC 0.71). Compensatory "), it("rpoC"),
  t(" mutations increased monotonically with resistance breadth (0.32 in pan-susceptible to 0.64 in pre-extensively-resistant isolates). In a within-host analysis of 226 patients with serial isolates, a susceptible isolate’s background predicted subsequent resistance with an AUROC of 0.66 (95% CI 0.60–0.73), lower than the cross-sectional estimate, as expected for a prospective task. "),
  bl("Conclusions. "),
  t("An interpretable, fully open pipeline reproduces deployed susceptibility prediction, extends it to quantitative MIC, and characterises a pre-resistance signal on the complete cohort; a within-host analysis confirms this signal is genuine but modest. Fully dated, prospective forecasting by longitudinal survival analysis on time-stamped isolates remains the principal next step."),
]));
children.push(new Paragraph({children:[new PageBreak()]}));

// ===================== 1. INTRODUCTION =====================
children.push(H("1. Introduction"));
children.push(P([
  t("Drug-resistant tuberculosis remains a leading cause of death from antimicrobial resistance. Effective treatment depends on knowing the susceptibility profile of the infecting isolate, but phenotypic DST for "),
  it("M. tuberculosis"), t(" is slow: the organism grows slowly, and solid-media culture with phenotypic testing routinely takes two to eight weeks. In high-burden settings the delay is consequential. In Thailand, a WHO high-burden country with an estimated tuberculosis incidence of 157 per 100,000 (approximately 113,000 cases in 2023), the median interval from specimen collection to physician review of a DST result has been measured at 114 days, and modelling indicates that more than 10,000 drug-resistant cases per year are managed as drug-susceptible. The same diagnostic latency contributes to mortality in melioidosis ("),
  it("Burkholderia pseudomallei"), t("), an intrinsically resistant pathogen causing an estimated 2,800 deaths annually in Thailand with case-fatality near 35%."),
]));
children.push(P([
  t("Whole-genome sequencing (WGS) can shorten this interval by inferring resistance directly from genotype. The relationship between specific mutations and resistance in "),
  it("M. tuberculosis"), t(" is comparatively well characterised, and curated catalogues and tools—including the WHO mutation catalogue, TB-Profiler, and Mykrobe—already translate sequence to a binary susceptibility call with high accuracy for first-line agents. The CRyPTIC consortium has further assembled matched genotypic and quantitative phenotypic data for thousands of isolates. Two limitations motivate the present work. First, most deployed pipelines report only resistant/susceptible, discarding the quantitative MIC that informs dosing. Second, catalogue-based prediction is, by construction, retrospective: it detects resistance that has already arisen and is silent on which currently-susceptible isolates are predisposed to acquire it. Recent work has shown that such pre-resistance signatures exist and can be recovered from genomic background."),
]));
children.push(P([
  t("We therefore conducted a reproducible analysis of the full CRyPTIC cohort with three aims: (i) to benchmark interpretable susceptibility prediction across first- and second-line drugs; (ii) to extend prediction from binary calls to quantitative MIC; and (iii) to characterise a pre-resistance signal by testing whether genomic background, excluding a drug’s own target gene, anticipates resistance to that drug. "),
  it("M. tuberculosis"), t(" was selected over other priority pathogens on the basis of data availability and mechanistic interpretability: querying the BV-BRC resource returned 2,137,022 laboratory susceptibility records for "),
  it("M. tuberculosis"), t(" against only 22 for "), it("B. pseudomallei"),
  t(", which is retained here as a motivating clinical context rather than a modelling target."),
]));

// ===================== 2. METHODS =====================
children.push(H("2. Methods"));
children.push(H2("2.1 Data"));
children.push(P([
  t("Genotype–phenotype pairs were obtained from the CRyPTIC compendium (release June 2022; reuse table), which provides, per isolate, a variant-call file (VCF) relative to the H37Rv reference (GenBank AL123456.3) and binary phenotypes with quantitative MICs measured on a 96-well broth-microdilution plate. Of 12,287 isolates with an available VCF, 12,245 (99.7%) were retrieved and analysed. Population-level resistance prevalences were additionally summarised from 225,000 laboratory records (44,163 isolates) in the BV-BRC "),
  it("genome_amr"), t(" resource (taxon 1773)."),
]));
children.push(H2("2.2 Variant translation and features"));
children.push(P([
  t("Each PASS or minor-allele variant was mapped to its gene and, for coding genes, translated to an amino-acid substitution using the H37Rv gene models, with correct handling of minus-strand genes, multi-nucleotide substitutions, and frameshift or nonsense changes (encoded as loss of function). Features were then restricted to WHO-catalogue resistance loci and codons—for example, the "),
  it("rpoB"), t(" rifampicin-resistance-determining region, "), it("katG"), t(" codon 315, the "),
  it("gyrA"), t(" quinolone-resistance-determining region, and "), it("rrs"),
  t(" position 1401—so that neutral phylogenetic variants (such as "), it("gyrA"),
  t(" Ser95Thr, present in most isolates) were excluded by construction. Background genes relevant to evolvability were additionally encoded: the compensatory RNA-polymerase subunit "),
  it("rpoC"), t(", the efflux regulator "), it("whiB7"), t(", the "), it("Rv1258c"),
  t(" efflux pump, and the mismatch-repair genes "), it("mutT1"), t(" and "), it("nucS"), t("."),
]));
children.push(H2("2.3 Models and evaluation"));
children.push(P([
  t("For binary susceptibility, gradient-boosted decision-tree classifiers were trained per drug on 70% of isolates and evaluated on the remaining 30%, stratified by phenotype. We report AUROC, sensitivity, specificity, very-major error (a resistant isolate called susceptible) and major error. Models were not informed of which locus corresponds to which drug, so the concordance of learned feature importances with established biology constitutes an internal validation. MIC was modelled by gradient-boosted regression of log2(MIC); performance is reported as essential agreement (prediction within one doubling dilution of the measured value), agreement within two dilutions, mean absolute error in dilutions, and the coefficient of determination. For the pre-resistance analysis, the feature set for each drug was reduced by removing that drug’s own target-gene features, and a classifier was trained on the residual background; a secondary analysis was restricted to isolates carrying no target-gene mutation. To test whether finer resolution improves prediction, we additionally scored every isolate against the full WHO catalogue version 2 (49,331 graded mutations, obtained in machine-readable form), matching each isolate's amino-acid changes to catalogue alleles graded “associated with resistance” and adding the resulting per-drug indicators to the codon-level feature set. All code, feature definitions, and trained-model parameters are provided for reproduction."),
]));
children.push(H2("2.4 Within-host longitudinal analysis"));
children.push(P([
  t("To move beyond cross-sectional inference, we used the 229 patients (226 after restriction to isolates with parsed genotypes) in CRyPTIC that contributed more than one sequenced isolate. For each patient and drug with at least one susceptible isolate, we defined a within-host acquisition outcome (1 if any isolate from that patient was resistant to the drug, 0 if all were susceptible) and predicted it from the genomic background of the patient’s susceptible isolate, with that drug’s own target-gene features set to zero. Units were pooled across drugs (n = 1,645; 138 acquisition events) with drug identity included as a covariate. Discrimination was estimated by patient-grouped five-fold cross-validation (no patient in both training and test folds) and a patient-level bootstrap 95% confidence interval (2,000 resamples); marker-level associations were tested by Fisher’s exact test. Because CRyPTIC does not provide collection dates, within-host susceptible/resistant discordance may reflect de novo acquisition, heteroresistance, or mixed infection; we therefore interpret this as a within-host emergence signal rather than a dated evolutionary trajectory."),
]));

// ===================== 3. RESULTS =====================
children.push(new Paragraph({pageBreakBefore:true,heading:HeadingLevel.HEADING_1,children:[new TextRun("3. Results")]}));
children.push(H2("3.1 Binary susceptibility prediction"));
children.push(P([
  t("Across nine drugs, susceptibility prediction achieved a mean AUROC of "),
  t(macroAUC), t(" with mean sensitivity "), t((meanSens*100).toFixed(0)+"%"),
  t(" and specificity "), t((meanSpec*100).toFixed(0)+"%"),
  t(" on held-out isolates (Table 2, Figure 1). Performance was highest for first-line agents (rifabutin AUROC 0.97, rifampicin 0.97, isoniazid 0.97) and lowest for kanamycin (0.89), consistent with the more diffuse genetic basis of aminoglycoside resistance. Very-major-error rates were correspondingly higher for kanamycin (30%) and ethionamide (25%) and are reported without adjustment. Adding per-allele indicators from the full WHO catalogue to the codon-level features changed mean AUROC by only +0.002 (0.948 to 0.950) and mean sensitivity by +0.003, with the largest gains for rifabutin and moxifloxacin; the codon-level representation therefore already captures the predictive content of the common catalogued alleles in this cohort."),
]));
children.push(table(
  ["Drug","Isolates","R (%)","AUROC","Sens.","Spec.","VME"],
  clfRows.map(d=>[{v:cap1(d.drug),italic:true}, d.n.toLocaleString(), d.pctR.toFixed(0),
    {v:d.auc.toFixed(2)}, (d.sens*100).toFixed(0)+"%", (d.spec*100).toFixed(0)+"%",
    {v:(d.vme*100).toFixed(0)+"%", color:d.vme>0.25?RES:undefined}]),
  [1900,1100,900,1100,1180,1180,1100]));
children.push(tableCap("2","Held-out (30%) susceptibility-prediction performance per drug. Sens., sensitivity; Spec., specificity; VME, very-major error (resistant isolate predicted susceptible)."));
children.push(img("figures/03_model_performance.png",560,308));
children.push(cap("1","Per-drug AUROC on the held-out test set with very-major-error rates; the dashed line marks AUROC = 0.90."));
children.push(P([
  t("Learned feature importances recovered the established molecular basis of resistance without supervision (Figure 2): rifampicin and rifabutin to the "),
  it("rpoB"), t(" resistance-determining region, isoniazid to "), it("katG"), t(" codon 315, ethambutol to "),
  it("embB"), t(" codon 306, the fluoroquinolones to the "), it("gyrA"),
  t(" quinolone-resistance-determining region, and the aminoglycosides to "), it("rrs"),
  t(". Off-diagonal contributions (for example "), it("rpoB"),
  t(" weighting for ethambutol) reflect co-resistance in multidrug-resistant strains rather than misattribution."),
]));
children.push(img("figures/05_feature_importance.png",470,382));
children.push(cap("2","Feature importance by drug. The diagonal indicates recovery of the canonical resistance locus for each drug from genotype alone."));

children.push(H2("3.2 Quantitative MIC prediction"));
children.push(P([
  t("Extending prediction from binary calls to the measured concentration, regression of log2(MIC) achieved a mean essential agreement of "),
  t((micEA*100).toFixed(0)+"%"),
  t(" (within one doubling dilution) and 93% within two dilutions across nine drugs (Table 3, Figure 3). Agreement was highest for amikacin (91%) and rifabutin (88%) and lowest for moxifloxacin (57%). Quantitative prediction is clinically relevant because the MIC, rather than a binary call, informs dose selection (for example, the use of high-dose isoniazid in the presence of low-level resistance)."),
]));
children.push(table(
  ["Drug","Essential agreement","Within ±2 dil.","R²"],
  mic.map(m=>[{v:cap1(m.drug),italic:true}, (+m.essential_agreement*100).toFixed(0)+"%",
    (+m.within_2dil*100).toFixed(0)+"%", (+m.R2).toFixed(2)]),
  [2600,2400,2400,1960]));
children.push(tableCap("3","Held-out MIC-regression performance. Essential agreement, prediction within one doubling dilution of the measured MIC."));
children.push(img("figures/07_mic.png",560,232));
children.push(cap("3","(A) MIC-prediction accuracy per drug: essential agreement (±1 dilution) and agreement within ±2 dilutions. (B) Predicted versus measured isoniazid MIC; the dashed line denotes equality and dotted lines ±1 dilution."));

children.push(H2("3.3 Pre-resistance and resistance amplification"));
children.push(P([
  t("To test whether genomic background anticipates resistance before the canonical mutation is present, each drug’s own target-gene features were withheld. The residual background predicted second-line resistance at AUROC 0.83–0.94 (Figure 4A). This figure must not be read as strain-intrinsic pre-resistance: a specificity analysis shows it is largely confounded by multidrug-resistant (MDR) status. Across six drugs and five random splits, a single binary indicator—whether the isolate is MDR (isoniazid- and rifampicin-resistant)—achieved a mean AUROC of 0.813 on its own, against 0.894 for the full 2,233-feature background model, so the entire high-dimensional model adds only 0.081 (95% CI 0.070–0.092), and roughly four-fifths of the above-chance discrimination is reproduced by that one variable. Stratifying to control for it, background discrimination falls to 0.820 within MDR isolates and 0.725 within non-MDR isolates. A residual signal above chance therefore remains, but the dominant driver is the well-established epidemiological observation that MDR strains are enriched for second-line resistance, not anticipation of resistance yet to arise. When the analysis was restricted to isolates carrying no target-gene mutation—those a catalogue would classify as susceptible—background features retained discrimination for the truly resistant minority (moxifloxacin AUROC 0.71, ethambutol 0.89), identifying phenotypic resistance not explained by catalogued target-gene variants. Consistent with established amplification dynamics, compensatory "),
  it("rpoC"), t(" mutations and efflux variants increased monotonically with resistance breadth, the former rising from 0.32 in pan-susceptible isolates to 0.64 in pre-extensively-resistant isolates (Figure 4B). The cohort also reproduced the reported ordering of acquisition: 74% of isoniazid-resistant isolates were rifampicin-resistant, against 5% of isoniazid-susceptible isolates."),
]));
children.push(img("figures/06_forecast.png",600,259));
children.push(cap("4","(A) Resistance anticipated from genomic background with the drug’s own target gene removed (orange) versus the full model (navy); red points denote prediction restricted to catalogue-susceptible isolates. (B) Prevalence of compensatory and efflux variants across escalating resistance categories."));

children.push(H2("3.4 Within-host resistance acquisition"));
children.push(P([
  t("The preceding analyses are cross-sectional, in which co-occurring resistance inflates apparent anticipation. We therefore tested a within-host, prospective-style design on the 226 patients with serial isolates. Predicting resistance acquisition from a susceptible isolate’s background (target gene excluded), evaluated with patient-grouped cross-validation, yielded an AUROC of 0.66 (95% confidence interval 0.60–0.73; 138 acquisition events across 226 patients; Figure 8). As expected, discrimination was substantially lower than the cross-sectional estimate (mean 0.87), consistent with the more difficult task of anticipating an event not yet manifest. The signal was attributable to the breadth of existing resistance rather than to specific evolvability markers: compensatory "),
  it("rpoC"), t(" mutations, mutator and efflux variants showed no individual enrichment in susceptible isolates that subsequently displayed resistance (Fisher exact p > 0.3 for all), a negative result we report without adjustment."),
]));
children.push(img("figures/08_longitudinal.png",600,260));
children.push(cap("8","Within-host prediction of resistance acquisition. (A) Receiver-operating characteristic for predicting, from a susceptible isolate’s background, whether the same patient harbours a resistant isolate (patient-grouped cross-validation). (B) The within-host (prospective) estimate is lower than the cross-sectional estimate, as expected; error bar denotes the patient-level bootstrap 95% confidence interval."));
children.push(H2("3.5 Resistance epidemiology and population structure"));
children.push(P([
  t("At the population level, resistance prevalences in the cohort (isoniazid 49%, rifampicin 39%) and a multidrug-resistance frequency of 27% among co-tested isolates are consistent with the enriched composition of CRyPTIC and with independent BV-BRC records (Figure 5). Co-resistance structure (Figure 6) and a genotype-distance phylogeny (Figure 7) illustrate the clustering of resistant lineages that underpins genomic surveillance."),
]));
children.push(img("figures/01_resistance_heatmap.png",380,450));
children.push(cap("5","Susceptibility profile of 90 representative isolates across nine drugs, ordered by resistance burden."));
children.push(img("figures/02_coresistance_map.png",440,380));
children.push(cap("6","Co-resistance: probability of resistance to each column drug given resistance to the row drug."));
children.push(img("figures/04_transmission_tree.png",390,425));
children.push(cap("7","Genotype-distance phylogeny of 60 isolates; multidrug-resistant isolates (isoniazid- and rifampicin-resistant) are indicated in red."));

// ===================== 4. DISCUSSION =====================
children.push(new Paragraph({pageBreakBefore:true,heading:HeadingLevel.HEADING_1,children:[new TextRun("4. Discussion")]}));
children.push(P([
  t("We present a reproducible, interpretable analysis of genome-based resistance prediction in "),
  it("M. tuberculosis"), t(" on the complete CRyPTIC cohort, with three findings: susceptibility prediction that recovers known biology at performance comparable to deployed tools; quantitative MIC prediction at "),
  t((micEA*100).toFixed(0)+"%"),
  t(" essential agreement; a background signal that predicts resistance independently of the target gene but which we show is largely attributable to multidrug-resistant status rather than to pre-resistance; and, in a within-host analysis of 226 patients with serial isolates, a modest prospective signal for resistance acquisition (AUROC 0.66, 95% CI 0.60–0.73)."),
]));
children.push(P([
  t("These results should be read in relation to existing work rather than as methodological novelty. Binary genotype-based susceptibility prediction for "),
  it("M. tuberculosis"), t(" is established and clinically deployed (TB-Profiler, Mykrobe, and the WHO catalogue), and our classifier reproduces that capability as a validated foundation. The existence of pre-resistance signatures has likewise been demonstrated previously by genome-wide survival analysis on longitudinally sampled isolates. The contribution here is therefore one of consolidation and extension: combining binary, quantitative, and pre-resistance analyses within a single open and interpretable pipeline on the full cohort, and quantifying the residual predictability of resistance among isolates that catalogue-based tools classify as susceptible—the population in which very-major errors occur. Such a residual signal is of direct surveillance value, flagging isolates for confirmatory DST that current pipelines would pass."),
]));
children.push(P([
  t("Several limitations qualify these findings. First, the principal analyses are cross-sectional, and co-occurring resistance inflates apparent anticipation to a degree we quantify explicitly: a single multidrug-resistance indicator reproduces roughly four-fifths of the above-chance discrimination attributed to the background model (Section 3.3). We report this because the confound is easy to miss and, left unexamined, would materially overstate the result; any future claim of pre-resistance should be accompanied by such a specificity analysis. The within-host analysis (Section 3.4) partially addresses the same issue and shows a modest prospective signal (AUROC 0.66), yet it too is constrained by the absence of collection dates and modest event numbers, so it cannot distinguish de novo acquisition from heteroresistance or mixed infection. A fully prospective forecast would require time-stamped, longitudinally sampled isolates with ancestral-state reconstruction and survival analysis, with external validation in regional cohorts such as available Thai serial-isolate collections. Second, features were defined at gene and codon resolution, with loss of function at diffuse loci ("),
  it("pncA"), t(", "), it("ethA"), t(", "), it("gidB"),
  t(") is encoded at gene level. We tested full per-allele integration against the WHO catalogue and found negligible benefit (mean AUROC +0.002); further gains would require cohorts enriched for rare resistance variants, not finer encoding of common ones. Third, pyrazinamide is absent from the CRyPTIC plate and was not modelled, and would require an external dataset with paired pncA genotype and pyrazinamide phenotype. Finally, the cohort is globally sourced; prospective re-validation on local isolates is required before any clinical application, and the tool is intended to augment, not replace, confirmatory DST."),
]));
children.push(P([
  t("Future work follows directly: assembly of time-stamped longitudinal and regional cohorts to convert the pre-resistance signature into a validated prospective forecast; a rare-variant-enriched cohort to quantify any benefit of per-allele over codon-level features; an external dataset to extend the model to pyrazinamide; and ingestion of long-read (nanopore) data to reduce turnaround further. Each is supported by the present pipeline without architectural change."),
]));

// ===================== 5. CONCLUSION =====================
children.push(H("5. Conclusion"));
children.push(P([
  t("Genome-based analysis of 12,245 "), it("M. tuberculosis"),
  t(" isolates reproduces deployed susceptibility prediction, extends it to quantitative MIC at "),
  t((micEA*100).toFixed(0)+"%"),
  t(" essential agreement, and recovers an interpretable pre-resistance signal that remains informative among catalogue-susceptible isolates. The work is fully reproducible on open data and is positioned as decision support and genomic surveillance for high-burden settings, with longitudinal prospective validation identified as the principal next step."),
]));

children.push(H("Data and code availability"));
children.push(P([ t("All data are openly available: the CRyPTIC compendium (EBI FTP, release June 2022) and the BV-BRC "),
  it("genome_amr"), t(" resource. All analysis code, feature definitions, figures, and trained-model parameters accompany this submission and reproduce every reported number.") ]));

// ===================== REFERENCES =====================
children.push(H("References"));
const refs=[
 "World Health Organization. Catalogue of mutations in Mycobacterium tuberculosis complex and their association with drug resistance, 2nd edn. Geneva: WHO; 2023.",
 "The CRyPTIC Consortium. A data compendium associating the genomes of 12,289 Mycobacterium tuberculosis isolates with quantitative resistance phenotypes to 13 antibiotics. PLoS Biol. 2022;20(8):e3001721.",
 "The CRyPTIC Consortium. Genomic signatures of pre-resistance in Mycobacterium tuberculosis. Nat Commun. 2021;12:7312.",
 "Walker TM, Kohl TA, Omar SV, et al. Whole-genome sequencing for prediction of Mycobacterium tuberculosis drug susceptibility and resistance: a retrospective cohort study. Lancet Infect Dis. 2015;15(10):1193-1202.",
 "Hunt M, Bradley P, Lapierre SG, et al. Antibiotic resistance prediction for Mycobacterium tuberculosis from genome sequence data with Mykrobe. Wellcome Open Res. 2019;4:191.",
 "Phelan JE, O'Sullivan DM, Machado D, et al. Integrating informatics tools and portable sequencing technology for rapid detection of resistance to anti-tuberculous drugs (TB-Profiler). Genome Med. 2019;11:41.",
 "Miotto P, Tessema B, Tagliani E, et al. A standardised method for interpreting the association between mutations and phenotypic drug resistance in Mycobacterium tuberculosis. Eur Respir J. 2017;50:1701354.",
 "Faksri K, Tan JH, Disratthakit A, et al. Whole-genome sequencing analysis of serially isolated multi-drug and extensively drug resistant Mycobacterium tuberculosis from Thai patients. PLoS One. 2016;11(8):e0160992.",
 "Olson RD, Assaf R, Brettin T, et al. Introducing the Bacterial and Viral Bioinformatics Resource Center (BV-BRC). Nucleic Acids Res. 2023;51(D1):D678-D689.",
 "World Health Organization. Global Tuberculosis Report 2024. Geneva: WHO; 2024.",
];
refs.forEach((r,i)=>children.push(new Paragraph({spacing:{after:60},
  children:[new TextRun({text:`${i+1}. `,bold:true,size:18}),new TextRun({text:r,size:18})]})));

// ===================== ASSEMBLE =====================
const doc=new Document({
  styles:{ default:{document:{run:{font:"Times New Roman",size:22}}},
    paragraphStyles:[
      {id:"Heading1",name:"Heading 1",basedOn:"Normal",next:"Normal",quickFormat:true,
       run:{size:26,bold:true,font:"Times New Roman",color:INK},
       paragraph:{spacing:{before:240,after:120},outlineLevel:0}},
      {id:"Heading2",name:"Heading 2",basedOn:"Normal",next:"Normal",quickFormat:true,
       run:{size:23,bold:true,italics:true,font:"Times New Roman",color:INK},
       paragraph:{spacing:{before:160,after:80},outlineLevel:1}},
    ]},
  sections:[{
    properties:{page:{size:{width:12240,height:15840},margin:{top:1440,right:1440,bottom:1440,left:1440}}},
    headers:{default:new Header({children:[new Paragraph({alignment:AlignmentType.RIGHT,
      children:[new TextRun({text:"Genome-based resistance prediction in M. tuberculosis",italics:true,size:16,color:"99a2b5"})]})]})},
    footers:{default:new Footer({children:[new Paragraph({alignment:AlignmentType.CENTER,
      children:[new TextRun({text:"",size:16,color:"99a2b5"}),new TextRun({children:[PageNumber.CURRENT],size:16,color:"99a2b5"})]})]})},
    children
  }]
});
Packer.toBuffer(doc).then(b=>{fs.writeFileSync("GenoTB_Manuscript.docx",b);
  console.log("wrote GenoTB_Manuscript.docx",b.length,"bytes");});
