# Study Summary

## Title

Integrative Evaluation of Computational Missense Variant Effect Predictors
in BRCA1 and BRCA2 Using Transcript-Aware Annotation and ClinVar Evidence

## Objective

To construct a reproducible GRCh38 BRCA1/BRCA2 missense variant cohort,
integrate five computational prediction tools, and evaluate their
discrimination and ranking relationships against ClinVar-derived
pathogenic and benign classifications.

## Cohort construction

- ClinVar BRCA1/BRCA2 GRCh38 candidate SNVs: 28,302.
- Variants with target-gene MANE Select transcript: 28,301.
- Variants requiring manual transcript review: 1.
- Primary MANE-transcript missense cohort: 16,193.
- BRCA1 missense variants: 6,374.
- BRCA2 missense variants: 9,819.

Transcript selection was restricted to the ClinVar target gene. A
non-target-gene transcript was not substituted when a target-gene
transcript was unavailable.

## Predictor integration

Five predictors were integrated:

- SIFT
- PolyPhen-2
- AlphaMissense
- CADD PHRED
- REVEL

All five scores were available for all 16,193 missense variants.

SIFT scores were direction-aligned as 1 minus the original score so that
higher values consistently represented greater predicted damage.

## Benchmark construction

The primary binary benchmark contained 1,546 variants:

- Pathogenic-side labels: 342.
- Benign-side labels: 1,204.

The high-confidence sensitivity benchmark contained 656 variants:

- Pathogenic-side labels: 252.
- Benign-side labels: 404.

The remaining 14,647 missense variants were not included in primary
binary evaluation because their ClinVar classifications were uncertain,
conflicting, missing, or otherwise unsuitable for unambiguous binary
labelling.

## Main performance results

| Predictor | Primary ROC-AUC | Primary AP | High-confidence ROC-AUC | High-confidence AP |
|---|---:|---:|---:|---:|
| AlphaMissense | 0.953 | 0.910 | 0.957 | 0.953 |
| CADD PHRED | 0.951 | 0.838 | 0.948 | 0.889 |
| REVEL | 0.928 | 0.847 | 0.927 | 0.907 |
| SIFT | 0.859 | 0.513 | 0.874 | 0.723 |
| PolyPhen-2 | 0.750 | 0.490 | 0.761 | 0.668 |

AlphaMissense had the highest observed ROC-AUC and average precision in
both benchmark cohorts. CADD PHRED showed similarly strong ROC-AUC.
These are descriptive comparisons; formal paired comparisons are still
required before asserting statistically significant differences.

## Correlation findings

The strongest direction-aligned Spearman relationship was between SIFT
and CADD PHRED (rho=0.793). AlphaMissense and CADD PHRED were also
moderately correlated (rho=0.666). The predictors were not perfectly
correlated, indicating differences in their variant rankings.

Correlation p-values that numerically underflowed to zero must not be
reported as literal zero probabilities.

## Gene-specific observations

Predictor performance varied between BRCA1 and BRCA2. For example,
PolyPhen-2 primary benchmark ROC-AUC was approximately 0.663 in BRCA1
and 0.915 in BRCA2. SIFT average precision was approximately 0.657 in
BRCA1 and 0.367 in BRCA2.

These subgroup observations are descriptive. They may reflect cohort
composition, label prevalence, evidence distribution, and sample size.
They do not by themselves establish intrinsic gene-specific model
performance.

## Limitations

1. ClinVar aggregate classifications are not independent experimental
   ground truth.
2. Only 1,546 of 16,193 missense variants had unambiguous binary labels.
3. The high-confidence cohort is a sensitivity analysis, not a
   replacement for the primary benchmark.
4. Predictor training datasets may overlap with ClinVar or related
   clinical evidence.
5. Average precision depends on pathogenic prevalence.
6. Gene-specific subgroup comparisons are descriptive.
7. Independent experimental validation and external-cohort evaluation
   are needed before clinical use.
8. Predictor score disagreement does not establish which prediction is
   biologically correct.

## Conclusion

This analysis establishes a reproducible, transcript-aware framework
for integrating five missense variant effect predictors across BRCA1
and BRCA2. AlphaMissense and CADD PHRED showed the strongest overall
discrimination against the selected ClinVar benchmark, while
gene-specific analyses revealed meaningful heterogeneity among tools.

The framework is suitable for research prioritization and further
validation. It is not a clinical diagnostic or variant-classification
system.
