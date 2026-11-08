# Integrative Functional Impact Assessment of BRCA1 and BRCA2 Missense Variants

**A reproducible computational genomics study integrating ClinVar, Ensembl VEP, and five variant-effect prediction algorithms.**

## Overview

This project investigates how multiple computational predictors characterize missense variants in **BRCA1** and **BRCA2**, and how their predictions compare with clinical significance classifications recorded in NCBI ClinVar.

The analysis combines transcript-aware variant annotation, ClinVar evidence, and five computational predictors:

- SIFT
- PolyPhen-2
- AlphaMissense
- CADD PHRED
- REVEL

The study evaluates predictor discrimination using ROC-AUC and average precision, examines agreement between predictors using Spearman correlation, and identifies variants for which predictors produce substantially different relative rankings.

The intended application is **research-level variant characterization and prioritization**, not clinical diagnosis or independent clinical classification.

## Research question

How consistently do five computational predictors characterize BRCA1 and BRCA2 missense variants, how well do their scores distinguish ClinVar pathogenic-side and benign-side classifications, and where do their relative rankings disagree?

### Objectives

1. Construct a GRCh38 BRCA1/BRCA2 candidate SNV dataset from ClinVar.
2. Annotate variants using Ensembl Variant Effect Predictor (VEP).
3. Select target-gene MANE Select transcripts, retaining explicit exceptions for manual review.
4. Define a transcript-aware missense cohort.
5. Construct transparent binary ClinVar benchmark subsets.
6. Integrate five predictor scores.
7. Evaluate discrimination, predictor correlations, gene-specific performance, and relative ranking disagreement.
8. Produce reproducible tables, figures, reports, and provenance records.

## Study design and cohort

### Variant acquisition and annotation

ClinVar's variant summary archive was used to identify BRCA1 and BRCA2 variants on GRCh38. Candidate SNVs were identified using the VCF-specific allele fields and genomic position.

The initial extraction yielded:

| Processing stage | Variants |
|---|---:|
| BRCA1/BRCA2 GRCh38 candidate SNVs | 28,311 |
| Canonical A/C/G/T SNVs retained | 28,302 |
| Ambiguous or invalid allele records excluded | 9 |

The 28,302 retained candidates were annotated using Ensembl VEP. Full annotation returned 28,302 variants and 867,142 transcript consequence records, with no failed variants.

### Transcript selection

Transcript selection was restricted to the ClinVar target gene:

- MANE Select was preferred.
- Canonical transcript was the fallback only when no target-gene MANE Select transcript was available.
- Off-target gene transcripts were not accepted as primary transcripts.
- Variants without an appropriate target-gene transcript were retained as exceptions and flagged for manual review.

The transcript selection produced 28,301 target-gene MANE Select records and one manual-review exception.

The exception was excluded from transcript-dependent missense analysis.

### Missense cohort

| Category | Variants |
|---|---:|
| BRCA1 missense | 6,374 |
| BRCA2 missense | 9,819 |
| **Total missense cohort** | **16,193** |
| Non-missense | 12,108 |
| Manual-review exception | 1 |

All 16,193 missense variants had complete values for the five integrated predictors.

## ClinVar benchmark construction

ClinVar classifications were retained in the annotated dataset, but not every classification was suitable for binary evaluation.

The primary binary benchmark included:

- Pathogenic
- Likely pathogenic
- Pathogenic/Likely pathogenic
- Benign
- Likely benign
- Benign/Likely benign

Uncertain, conflicting, missing, unprovided, and other ambiguous classifications were excluded from primary binary evaluation. Their records and exclusion reasons were retained in the annotated dataset.

| Benchmark | Total | Pathogenic-side | Benign-side | Pathogenic prevalence |
|---|---:|---:|---:|---:|
| Full missense cohort | 16,193 | 342 eligible labels | 1,204 eligible labels | Not a complete binary benchmark |
| Primary benchmark | 1,546 | 342 | 1,204 | 22.12% |
| High-confidence benchmark | 656 | 252 | 404 | 38.41% |

The primary benchmark excludes 14,647 missense variants from binary evaluation.

The high-confidence benchmark is a stricter sensitivity analysis based on stronger ClinVar review evidence. It is not a replacement for the primary benchmark.

**Important:** ClinVar classifications are clinical assertions with varying evidence and review status. They are not an independent experimental gold standard.

## Predictors

| Predictor | Score used | Direction |
|---|---|---|
| SIFT | SIFT score | Lower = more damaging |
| PolyPhen-2 | PolyPhen-2 score | Higher = more damaging |
| AlphaMissense | Pathogenicity score | Higher = more pathogenic |
| CADD | PHRED score | Higher = more deleterious |
| REVEL | REVEL score | Higher = more pathogenic |

SIFT was direction-transformed as `1 - SIFT` for analyses requiring a common direction. This transformation does not imply that scores from different algorithms are calibrated or directly interchangeable.

CADD PHRED was used for primary predictor evaluation. CADD raw scores were retained separately and were not constrained to the interval [0, 1].

AlphaMissense matching was performed using exact genomic variant coordinates and alleles. Protein substitutions were also checked for consistency. Source transcript identifiers were retained because transcript conventions differ between resources.

## Main results

### Primary ClinVar benchmark

| Predictor | ROC-AUC | 95% CI | Average precision |
|---|---:|---:|---:|
| AlphaMissense | 0.953 | 0.937?0.967 | 0.910 |
| CADD PHRED | 0.951 | 0.938?0.961 | 0.838 |
| REVEL | 0.928 | 0.910?0.946 | 0.847 |
| SIFT | 0.859 | 0.843?0.874 | 0.513 |
| PolyPhen-2 | 0.750 | 0.720?0.779 | 0.490 |

### High-confidence ClinVar benchmark

| Predictor | ROC-AUC | 95% CI | Average precision |
|---|---:|---:|---:|
| AlphaMissense | 0.957 | 0.938?0.973 | 0.953 |
| CADD PHRED | 0.948 | 0.931?0.964 | 0.889 |
| REVEL | 0.927 | 0.904?0.949 | 0.907 |
| SIFT | 0.874 | 0.851?0.896 | 0.723 |
| PolyPhen-2 | 0.761 | 0.724?0.797 | 0.668 |

AlphaMissense and CADD PHRED showed the strongest ROC-AUC values in both benchmarks. AlphaMissense had the highest observed average precision in both.

These are descriptive comparisons. Overlapping confidence intervals do not establish statistical equivalence or superiority; formal paired comparisons would be needed to test differences between predictors.

Average precision depends on pathogenic prevalence, so values should be interpreted within their respective benchmark populations.

## Figures

### ROC curves

**Primary benchmark**

![ROC curves for the primary benchmark](results/figures/roc_curves_primary.png)

**High-confidence benchmark**

![ROC curves for the high-confidence benchmark](results/figures/roc_curves_high_confidence.png)

### Precision-recall curves

**Primary benchmark**

![Precision-recall curves for the primary benchmark](results/figures/pr_curves_primary.png)

**High-confidence benchmark**

![Precision-recall curves for the high-confidence benchmark](results/figures/pr_curves_high_confidence.png)

The horizontal reference in each precision-recall plot represents the pathogenic prevalence of that benchmark.

### Predictor agreement

![Spearman correlation heatmap](results/figures/predictor_correlation_heatmap.png)

Spearman correlations were calculated across the full 16,193-variant missense cohort using direction-aligned scores.

The strongest observed correlation was between SIFT and CADD PHRED (rho = 0.793). Correlation indicates ranking association, not agreement in calibration or biological interpretation.

### Relative predictor ranking disagreement

![Predictor ranking disagreement](results/figures/predictor_rank_disagreement.png)

The distribution represents the range of within-predictor percentile ranks assigned to each variant. It measures relative ranking spread, not binary classification disagreement or prediction correctness.

### Direction-aligned score distributions

![Direction-aligned percentile distributions](results/figures/predictor_score_distributions.png)

Scores are converted to within-predictor percentile ranks to permit visual comparison of relative ranking distributions. These percentiles are not calibrated probabilities and do not imply equivalent score scales.

## Statistical analysis

- ROC-AUC measures discrimination across classification thresholds.
- Average precision summarizes the precision-recall curve and is sensitive to class prevalence.
- Confidence intervals were estimated using 2,000 stratified bootstrap resamples, with random seed 42.
- Spearman rank correlation was used to assess monotonic association between predictor scores.
- Gene-specific metrics were calculated separately for BRCA1 and BRCA2.
- Relative ranking disagreement was quantified using the range and standard deviation of within-predictor percentile ranks.

The analysis preserves predictor missingness rather than converting missing values into benign or zero scores. The integrated predictor dataset used in the final evaluation contained complete values for all five predictors.

## Repository structure

```text
ncbi-functional-variant-analysis/
??? config/
??? data/
?   ??? raw/                 # Downloaded source data (not committed)
?   ??? processed/           # Derived datasets (not committed)
?   ??? metadata/            # Provenance, manifests, audit reports
??? docs/
?   ??? references.md
?   ??? study_summary.md
?   ??? research_writeup.md
??? notebooks/
??? results/
?   ??? figures/
?   ??? reports/
?   ??? tables/
??? src/
?   ??? clinvar_downloader.py
?   ??? vep_full_annotation.py
?   ??? select_primary_transcripts.py
?   ??? classify_primary_missense.py
?   ??? build_clinvar_benchmark.py
?   ??? alphamissense_downloader.py
?   ??? extract_alphamissense_cohort.py
?   ??? integrate_alphamissense.py
?   ??? vep_full_predictor_annotation.py
?   ??? integrate_all_predictors.py
?   ??? run_predictor_evaluation.py
??? tests/
??? workflows/
??? .gitignore
??? README.md
??? requirements.txt
```

The repository also contains audit and pilot scripts. The structure above highlights the main processing stages rather than every utility.

## Reproducibility

### Environment

The project was developed with Python 3.10.7 in a virtual environment.

From the repository root, create and activate the environment:

```powershell
C:\Python310\python.exe -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

### Run the evaluation

The final evaluation consumes the integrated processed dataset:

```text
data/processed/brca1_brca2_all_predictors_integrated.csv
```

Run:

```powershell
python .\src\run_predictor_evaluation.py
```

The evaluation script generates the performance tables, correlation outputs, ranking-disagreement outputs, figures, and evaluation reports under `results/` and `data/metadata/`.

The upstream acquisition and annotation stages require their source files and, where applicable, access to external services. See the scripts and provenance manifests under `data/metadata/` for source details.

### Reproducibility considerations

- Random seed: 42.
- Bootstrap iterations: 2,000.
- Genome assembly: GRCh38.
- Variant matching uses genomic alleles and coordinates.
- Transcript selection is restricted to the target gene.
- Raw and processed data are excluded from Git by design.
- Source checksums and acquisition details are recorded in metadata manifests where available.

## Data availability and provenance

The repository tracks scripts, documentation, metadata, tables, figures, and reports. Large raw and processed datasets are intentionally excluded from version control.

Key source resources:

- [NCBI ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [Ensembl Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html)
- [AlphaMissense source repository](https://github.com/google-deepmind/alphamissense)
- [AlphaMissense Zenodo record](https://zenodo.org/records/8208688)

The AlphaMissense source file header and a newer Zenodo record contain differing licence declarations. This discrepancy is documented in the download manifest and remains unresolved. Do not assume commercial-use rights without checking the applicable terms with the relevant institutional or legal team.

## Limitations

1. **ClinVar evidence heterogeneity:** Clinical classifications differ in review status, evidence quality, and submission history.
2. **Benchmark coverage:** Only 1,546 of 16,193 missense variants entered the primary binary benchmark.
3. **Evidence dependence:** Predictor training datasets may overlap with ClinVar or related clinical evidence, creating potential dependence or leakage.
4. **Transcript representation:** Transcript and protein conventions differ between resources; genomic allele matching was therefore central to integration.
5. **Predictor calibration:** Scores from different algorithms are not directly comparable as probabilities.
6. **Gene subgroup differences:** BRCA1 and BRCA2 have different benchmark sizes and pathogenic prevalence; subgroup metrics are descriptive.
7. **No independent experimental validation:** This study evaluates agreement with ClinVar assertions, not experimentally measured functional effects.
8. **No clinical validation:** The analysis does not establish clinical utility or support patient-level diagnostic decisions.

## Responsible interpretation

The results support comparative research assessment and prioritization of variants for further investigation. They do not establish a clinical classification, replace expert review, or demonstrate that any predictor is universally superior.

## References

See [docs/references.md](docs/references.md) for the bibliography and [docs/research_writeup.md](docs/research_writeup.md) for the manuscript-style study draft.

## License

A repository licence has not been specified. Source datasets and third-party prediction resources remain subject to their own terms. Check those terms before redistribution or commercial use.
