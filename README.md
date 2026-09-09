# GenoTB — genome-based drug-susceptibility, MIC, and pre-resistance analysis in *Mycobacterium tuberculosis*

A reproducible analysis of the full **CRyPTIC** cohort (12,245 clinical isolates with matched genome and
96-well broth-microdilution phenotypes). Each isolate's variant calls are translated to amino-acid changes,
restricted to WHO-catalogue resistance loci, and used to predict (i) binary susceptibility, (ii) the
minimum inhibitory concentration (MIC), and (iii) a pre-resistance / resistance-acquisition signal.

## Scope
This is a rigorous **reproduction and consolidation**, not new methodology. Binary genome→resistance
prediction is established and deployed (TB-Profiler, Mykrobe, the WHO catalogue); the pre-resistance
signal is published (CRyPTIC Consortium, *Nat Commun* 2021). The value here is a single open, interpretable
pipeline that (a) benchmarks susceptibility prediction on the full cohort, (b) extends it to quantitative
MIC, and (c) tests a pre-resistance signal both cross-sectionally and within-host. See the manuscript
(`GenoTB_Manuscript.docx`) for full positioning and limitations.

## Results 
| Analysis | Result |
|---|---|
| Susceptibility, 9 drugs | mean AUROC **0.94**, specificity 97% |
| MIC regression | **79%** essential agreement (±1 dilution), 93% within ±2 |
| Pre-resistance (background, target gene removed) | AUROC 0.83–0.94; catalogue-susceptible MXF 0.71 |
| Within-host acquisition (226 serial-isolate patients) | AUROC **0.66** (95% CI 0.60–0.73) |
| Per-allele WHO catalogue vs codon-level | +0.002 AUROC (negligible — tested) |

## Pipeline
```bash
./dl.sh                        # download 12,245 CRyPTIC VCFs (resumable) — not in repo
python3 parse_v3.py            # variant → amino-acid features (aa_engine.py + catalogue.py)
python3 classifier_full_aa.py  # susceptibility classifiers
python3 mic_regression.py      # MIC regression
python3 forecast_emergence.py  # cross-sectional pre-resistance
python3 longitudinal_acquisition.py   # within-host acquisition
python3 parse_who.py && python3 compare_who.py   # per-allele WHO catalogue test
python3 viz_full.py; python3 viz_mic.py; python3 viz_forecast.py; python3 viz_longitudinal.py
node build_proposal.js         # -> GenoTB_Manuscript.docx
```

## Key files
- `aa_engine.py` — variant→amino-acid translation (H37Rv gene models)
- `catalogue.py` / `who_catalogue.py` — codon-restricted features / per-allele WHO matching
- `classifier_full_aa.py`, `mic_regression.py`, `forecast_emergence.py`, `longitudinal_acquisition.py`
- `GenoTB_Manuscript.docx`, `abstract.md`, `figures/`, `demo.html`

## Data
Raw genomes/phenotypes are **not** committed (410 MB of VCFs); they are openly available and re-downloaded
by `dl.sh` from the CRyPTIC compendium (EBI FTP, release June 2022). The WHO catalogue is fetched from the
TB-Profiler `tbdb` database. Only code, figures, the manuscript, and small result tables are tracked.

## Data availability & license
CRyPTIC compendium (EBI FTP), BV-BRC `genome_amr`, WHO catalogue v2. Research prototype — decision support
and surveillance, **not** a replacement for confirmatory DST.
