# Integrative Functional Impact Assessment of BRCA1 and BRCA2 Missense Variants Using Five Computational Predictors

## Abstract

### Background

Computational prediction algorithms are widely used to prioritize missense variants for further investigation. However, different algorithms use distinct training data, features, score conventions, and modelling strategies. Their predictions may therefore differ, and their performance depends on the reference evidence used for evaluation.

### Objective

This study aimed to construct a transcript-aware BRCA1 and BRCA2 missense variant cohort, integrate five computational predictors, evaluate their discrimination against selected ClinVar clinical significance assertions, and characterize agreement and relative ranking disagreement between predictors.

### Methods

ClinVar GRCh38 records were used to identify candidate BRCA1 and BRCA2 single-nucleotide variants. Variants were annotated with Ensembl Variant Effect Predictor (VEP), retaining transcript-level consequences and selecting target-gene MANE Select transcripts where available. A missense cohort of 16,193 variants was assembled. Scores from SIFT, PolyPhen-2, AlphaMissense, CADD PHRED, and REVEL were integrated. A primary binary benchmark included pathogenic-side and benign-side ClinVar classifications, while uncertain, conflicting, missing, and other ambiguous classifications were excluded from binary evaluation. A stricter high-confidence benchmark was used for sensitivity analysis. Predictor discrimination was assessed using ROC-AUC and average precision, with 2,000 stratified bootstrap resamples for confidence intervals. Spearman correlation and within-predictor percentile ranks were used to assess predictor agreement and ranking disagreement.

### Results

The final missense cohort contained 6,374 BRCA1 and 9,819 BRCA2 variants. The primary benchmark included 1,546 variants (342 pathogenic-side and 1,204 benign-side); the high-confidence benchmark included 656 variants (252 pathogenic-side and 404 benign-side). In the primary benchmark, AlphaMissense achieved an ROC-AUC of 0.953 and average precision of 0.910, while CADD PHRED achieved 0.951 and 0.838, respectively. REVEL achieved 0.928 and 0.847, SIFT 0.859 and 0.513, and PolyPhen-2 0.750 and 0.490. The high-confidence benchmark showed a similar overall pattern. The strongest pairwise Spearman correlation was observed between SIFT and CADD PHRED (rho = 0.793).

### Conclusion

The five predictors showed different levels of discrimination against the selected ClinVar benchmark, and their relative rankings were not interchangeable. AlphaMissense and CADD PHRED showed the strongest observed ROC-AUC values in the evaluated cohorts. These findings should be interpreted in the context of ClinVar evidence heterogeneity, benchmark selection, possible training-data dependence, and the absence of independent experimental validation. The results are intended for research prioritization rather than clinical decision-making.

**Keywords:** BRCA1; BRCA2; missense variants; ClinVar; Ensembl VEP; SIFT; PolyPhen-2; AlphaMissense; CADD; REVEL; variant prioritization.

---

## 1. Introduction

Missense variants alter protein-coding sequences and may affect protein structure, stability, interactions, or function. Their interpretation is challenging because the presence of an amino-acid substitution does not by itself establish a functional or clinical consequence.

Computational prediction methods provide complementary approaches for estimating the potential impact of amino-acid substitutions. These methods differ in their input features, evolutionary information, training strategies, and output scales. Consequently, predictions from different algorithms may agree for some variants and diverge for others.

BRCA1 and BRCA2 are clinically important cancer susceptibility genes with extensive variant records in ClinVar. Their large variant collections provide an opportunity to examine how computational predictors characterize missense variation in a defined genomic context.

This study integrates five prediction methods?SIFT, PolyPhen-2, AlphaMissense, CADD, and REVEL?with ClinVar and transcript-aware VEP annotation. Rather than treating one algorithm as a universal reference, the analysis evaluates their discrimination against selected ClinVar assertions and quantifies their agreement and relative ranking differences.

### 1.1 Objectives

The objectives were to:

1. Construct a GRCh38 BRCA1/BRCA2 candidate SNV collection from ClinVar.
2. Annotate variants and transcript consequences using Ensembl VEP.
3. Select target-gene MANE Select transcripts for primary transcript-level analysis.
4. Derive a missense variant cohort.
5. Define primary and high-confidence ClinVar binary benchmark subsets.
6. Integrate five computational predictor scores.
7. Compare predictor discrimination using ROC-AUC and average precision.
8. Examine predictor correlations and relative ranking disagreement.
9. Describe BRCA1 and BRCA2 subgroup performance.

---

## 2. Materials and methods

### 2.1 ClinVar variant acquisition

The ClinVar variant summary archive was used as the source for candidate BRCA1 and BRCA2 variants. Records were filtered to the GRCh38 assembly and candidate single-nucleotide substitutions.

A data-field issue was identified during extraction: the generic `ReferenceAllele` and `AlternateAllele` columns contained unavailable values for the candidate records. The extraction was corrected to use the VCF-specific allele fields, `ReferenceAlleleVCF` and `AlternateAlleleVCF`, together with `PositionVCF`.

Genomic variant identity was represented using chromosome, position, reference allele, and alternate allele. Only canonical A, C, G, and T alleles were retained for the primary candidate SNV collection.

This process yielded 28,302 canonical candidate SNVs from 28,311 candidate records. Nine ambiguous or invalid allele records were excluded.

### 2.2 VEP annotation and transcript selection

The candidate variants were annotated using Ensembl VEP against GRCh38. The annotation included transcript consequences and available transcript, protein, and prediction fields.

Full annotation returned 28,302 variants and 867,142 transcript consequence records, with no failed variants.

Transcript selection was restricted to the ClinVar target gene. MANE Select was preferred. A canonical target-gene transcript was permitted as a fallback only where a target-gene MANE Select transcript was unavailable. Off-target transcripts were excluded from primary transcript selection.

One variant lacked an appropriate target-gene transcript and was retained as a manual-review exception. It was excluded from transcript-dependent missense classification.

### 2.3 Missense cohort definition

Transcript consequence annotations were used to identify missense variants in the selected primary transcript.

The resulting cohort comprised:

- 6,374 BRCA1 missense variants.
- 9,819 BRCA2 missense variants.
- 16,193 missense variants in total.

The remaining 12,108 variants were classified as non-missense within the primary transcript analysis. One additional variant remained flagged for manual review.

### 2.4 ClinVar benchmark construction

ClinVar clinical significance labels were retained in the annotated cohort. For binary evaluation, the following pathogenic-side classifications were grouped together:

- Pathogenic.
- Likely pathogenic.
- Pathogenic/Likely pathogenic.

The following benign-side classifications were grouped together:

- Benign.
- Likely benign.
- Benign/Likely benign.

Uncertain significance, conflicting classifications, missing classifications, not-provided classifications, and other ambiguous categories were excluded from primary binary evaluation. Exclusion reasons were preserved.

The primary binary benchmark contained 1,546 variants: 342 pathogenic-side and 1,204 benign-side.

A high-confidence benchmark of 656 variants (252 pathogenic-side and 404 benign-side) was used as a sensitivity analysis based on stronger ClinVar review evidence.

The primary benchmark represented a selected subset of the missense cohort. A total of 14,647 missense variants were not included in the primary binary benchmark.

ClinVar assertions were treated as reference clinical evidence for comparative evaluation, not as independent experimental measurements of variant function.

### 2.5 Predictor integration

Five computational predictors were integrated:

| Predictor | Score | Interpretation direction |
|---|---|---|
| SIFT | SIFT score | Lower values indicate greater predicted impact |
| PolyPhen-2 | PolyPhen-2 score | Higher values indicate greater predicted impact |
| AlphaMissense | Pathogenicity score | Higher values indicate greater predicted pathogenicity |
| CADD | PHRED score | Higher values indicate greater predicted deleteriousness |
| REVEL | REVEL score | Higher values indicate greater predicted pathogenicity |

SIFT scores were transformed using `1 - score` for analyses requiring a consistent direction. Higher transformed values therefore represented greater predicted damaging impact.

The primary CADD analysis used PHRED scores. Raw CADD scores were retained separately.

AlphaMissense records were integrated using genomic coordinates and reference/alternate alleles. Protein substitutions were checked for consistency. Differences in transcript identifiers between resources were retained as provenance rather than treated as evidence of genomic mismatch.

All five predictors had complete values for the 16,193-variant integrated missense cohort.

### 2.6 Predictor performance evaluation

Each predictor was evaluated separately in the primary and high-confidence benchmarks.

The principal metrics were:

- **ROC-AUC:** discrimination across false-positive and true-positive rate thresholds.
- **Average precision:** summary of precision-recall performance.

Ninety-five percent percentile confidence intervals were estimated using 2,000 stratified bootstrap resamples. Pathogenic and benign classes were resampled separately to preserve class counts in each replicate. The random seed was 42.

Average precision was interpreted in relation to the pathogenic prevalence of each benchmark.

### 2.7 Predictor agreement

Spearman rank correlations were calculated between direction-aligned predictor scores across the full 16,193-variant missense cohort.

Pairwise correlation coefficients were used to characterize monotonic ranking association. They were not interpreted as evidence of calibration, equivalent biological meaning, or clinical agreement.

### 2.8 Relative ranking disagreement

For each predictor, variants were assigned within-predictor percentile ranks. The range between the highest and lowest percentile rank across the five predictors was calculated for each variant.

This measure represents the spread of relative predictor rankings. It does not represent binary classification disagreement, prediction correctness, or experimentally measured functional variation.

### 2.9 Gene-specific analyses

Performance metrics were calculated separately for BRCA1 and BRCA2 within both benchmark subsets.

These subgroup analyses were descriptive. Differences in benchmark size and pathogenic prevalence were considered when interpreting average precision and subgroup results.

---

## 3. Results

### 3.1 Cohort construction

The ClinVar extraction identified 28,311 BRCA1/BRCA2 GRCh38 candidate SNV records. After excluding nine ambiguous or invalid allele records, 28,302 canonical SNVs were retained.

VEP successfully annotated all 28,302 variants, producing 867,142 transcript consequence records. Target-gene MANE Select transcript selection was possible for 28,301 variants. One variant was retained as a manual-review exception because no appropriate target-gene transcript was available.

The final transcript-aware missense cohort contained 16,193 variants, comprising 6,374 BRCA1 and 9,819 BRCA2 missense variants.

### 3.2 Benchmark composition

The primary benchmark included 1,546 variants, with 342 pathogenic-side and 1,204 benign-side assertions. Pathogenic prevalence was 22.12%.

The high-confidence benchmark included 656 variants, with 252 pathogenic-side and 404 benign-side assertions. Pathogenic prevalence was 38.41%.

The difference in prevalence is important when interpreting average precision between benchmarks.

### 3.3 Primary benchmark performance

| Predictor | ROC-AUC | 95% confidence interval | Average precision |
|---|---:|---:|---:|
| AlphaMissense | 0.9525 | 0.9365?0.9672 | 0.9101 |
| CADD PHRED | 0.9507 | 0.9382?0.9614 | 0.8379 |
| REVEL | 0.9279 | 0.9100?0.9463 | 0.8472 |
| SIFT | 0.8587 | 0.8430?0.8737 | 0.5133 |
| PolyPhen-2 | 0.7499 | 0.7204?0.7792 | 0.4900 |

AlphaMissense and CADD PHRED had the highest observed ROC-AUC values. AlphaMissense also had the highest average precision. REVEL showed strong discrimination, while SIFT and PolyPhen-2 had lower average precision in this benchmark.

The confidence intervals for AlphaMissense and CADD PHRED overlap. The observed difference in ROC-AUC should therefore not be presented as proof of statistical superiority.

### 3.4 High-confidence benchmark performance

| Predictor | ROC-AUC | 95% confidence interval | Average precision |
|---|---:|---:|---:|
| AlphaMissense | 0.9569 | 0.9377?0.9732 | 0.9532 |
| CADD PHRED | 0.9478 | 0.9308?0.9637 | 0.8890 |
| REVEL | 0.9268 | 0.9036?0.9487 | 0.9071 |
| SIFT | 0.8736 | 0.8505?0.8961 | 0.7227 |
| PolyPhen-2 | 0.7612 | 0.7238?0.7972 | 0.6676 |

The overall pattern was consistent with the primary benchmark. AlphaMissense had the highest observed ROC-AUC and average precision, followed by CADD PHRED and REVEL in ROC-AUC.

The increase in average precision relative to the primary benchmark must be interpreted in light of the higher pathogenic prevalence in the high-confidence subset.

### 3.5 Predictor correlations

The strongest observed pairwise Spearman correlation was between SIFT and CADD PHRED (rho = 0.7928).

Other notable correlations included:

| Predictor pair | Spearman rho |
|---|---:|
| SIFT ? CADD PHRED | 0.7928 |
| AlphaMissense ? CADD PHRED | 0.6656 |
| PolyPhen-2 ? CADD PHRED | 0.6642 |
| SIFT ? AlphaMissense | 0.6349 |
| SIFT ? PolyPhen-2 | 0.5753 |
| CADD PHRED ? REVEL | 0.5167 |

The correlations indicate varying degrees of monotonic association. They do not establish that predictors make identical variant classifications or provide equivalent functional measurements.

### 3.6 Gene-specific performance

In the primary benchmark, BRCA1 contained 689 variants (224 pathogenic-side; prevalence 32.51%), while BRCA2 contained 857 variants (118 pathogenic-side; prevalence 13.77%).

The high-confidence gene-specific subsets contained 347 BRCA1 variants (158 pathogenic-side) and 309 BRCA2 variants (94 pathogenic-side).

The predictor metrics differed between gene subsets. These differences should be interpreted descriptively because the sample sizes and class prevalences differed. They do not establish intrinsic gene-specific superiority of one predictor.

### 3.7 Relative ranking disagreement

The within-predictor percentile analysis demonstrated that variants can occupy substantially different relative positions across prediction algorithms.

This result supports the use of multiple predictors as complementary evidence sources. It does not show that the most discordant variants are necessarily incorrectly predicted or functionally important.

---

## 4. Discussion

### 4.1 Comparative predictor performance

The evaluation showed that the five predictors did not exhibit identical performance against the selected ClinVar benchmarks. AlphaMissense and CADD PHRED had the strongest observed ROC-AUC values, while AlphaMissense achieved the highest average precision in both benchmark subsets.

REVEL also showed strong discrimination. SIFT and PolyPhen-2 had lower performance in the primary benchmark, although their average precision increased in the high-confidence subset.

These observations are specific to the selected BRCA1/BRCA2 missense cohort, the ClinVar benchmark definitions, and the evaluation design. They should not be generalized to all genes, variant types, or clinical settings.

### 4.2 Interpretation of predictor agreement

The correlation analysis showed that predictor scores were associated to varying degrees. The relatively strong SIFT?CADD correlation suggests that these methods often assign similar relative rankings in this cohort, despite their different score definitions.

Moderate or weaker correlations between other pairs indicate that predictor rankings can diverge. Such disagreement may arise from differences in model design, input information, training data, and score calibration.

A combined interpretation framework may therefore be more informative than reliance on a single predictor. However, a consensus score should not be constructed without a separately justified integration method and validation strategy.

### 4.3 Primary and high-confidence benchmarks

The high-confidence subset was used to assess whether observed performance patterns persisted under stricter ClinVar evidence criteria.

The broadly similar ROC-AUC ordering provides descriptive evidence of consistency across the two benchmark definitions. The different pathogenic prevalence means average precision values cannot be compared without considering the benchmark composition.

Neither benchmark eliminates the limitations of using ClinVar assertions as reference labels.

### 4.4 Potential evidence dependence

Computational predictors may have been trained or calibrated using datasets that overlap with ClinVar or related clinical evidence resources. This creates the possibility of evidence dependence or information leakage.

The present analysis does not establish the absence of such overlap. Future work should examine predictor training-data provenance and, where possible, evaluate performance on temporally or evidentially independent variant sets.

### 4.5 Research prioritization implications

The observed differences in performance and relative rankings support treating computational predictions as complementary evidence rather than interchangeable outputs.

Variants with strong agreement across methods may be useful candidates for further investigation. Variants with marked ranking disagreement may warrant closer examination of transcript context, model assumptions, clinical evidence, and experimental literature.

Neither pattern is sufficient to establish a clinical classification or a measured functional effect.

---

## 5. Limitations

Several limitations should be considered.

1. **ClinVar is not experimental ground truth.** The benchmark reflects clinical assertions with heterogeneous evidence and review status.
2. **Benchmark coverage is limited.** Only 1,546 of the 16,193 missense variants entered the primary binary benchmark.
3. **Potential training-data overlap.** Predictor training datasets may overlap with ClinVar or related resources.
4. **Transcript representation differs.** Transcript identifiers and protein conventions vary across annotation resources, even where genomic alleles match.
5. **Scores are not calibrated equivalently.** Raw predictor values cannot be compared directly across algorithms.
6. **Average precision is prevalence-dependent.** Different pathogenic prevalence complicates comparisons between benchmark subsets.
7. **Subgroup analyses are descriptive.** BRCA1 and BRCA2 have different sample sizes and class distributions.
8. **No independent experimental validation was performed.** The study does not measure protein function or clinical outcomes.
9. **No paired statistical comparison was performed.** The reported confidence intervals do not by themselves establish statistically significant differences between predictors.
10. **No clinical utility evaluation was performed.** The study does not assess clinical decision-making or patient outcomes.

---

## 6. Conclusion

This study established a transcript-aware BRCA1/BRCA2 missense variant cohort and integrated five computational prediction algorithms with ClinVar evidence.

AlphaMissense and CADD PHRED showed the strongest observed ROC-AUC values in the primary and high-confidence benchmarks, while AlphaMissense achieved the highest observed average precision. Predictor correlations and relative ranking analyses demonstrated that the five methods provide overlapping but non-identical information.

The findings support comparative, multi-predictor research prioritization while emphasizing the limitations of ClinVar-derived benchmarks, potential evidence dependence, and the absence of independent experimental validation. The outputs should not be interpreted as clinical classifications.

---

## 7. Data and code availability

The analysis scripts, documentation, result tables, figures, and provenance metadata are maintained in this repository.

Large raw and processed datasets are excluded from version control. Source acquisition details and checksums are documented in the metadata manifests where available.

The AlphaMissense source file contains a licence declaration that differs from a newer Zenodo record. This discrepancy remains unresolved and requires review before commercial use or redistribution.

## 8. References

See [`references.md`](references.md) for the bibliography, including the primary publications for VEP, ClinVar, SIFT, PolyPhen-2, REVEL, CADD, and AlphaMissense.
