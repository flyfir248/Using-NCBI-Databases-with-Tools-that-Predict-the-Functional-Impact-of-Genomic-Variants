# Integrative Functional Impact Assessment

## Study cohort

- Total missense variants: **16,193**
- Primary ClinVar benchmark: **1,546**
- High-confidence benchmark: **656**

The primary benchmark includes unambiguous benign-side and pathogenic-side ClinVar classifications. The high-confidence benchmark is a review-status-restricted sensitivity analysis.

## Predictor evaluation

Higher direction-aligned scores indicate greater predicted deleteriousness or pathogenicity. SIFT was transformed as `1 - SIFT score`; the other four scores were retained in their original direction.

ROC-AUC and average precision are reported with stratified bootstrap 95% confidence intervals.

### Primary benchmark

| cohort  | predictor     | n    | pathogenic_count | benign_count | pathogenic_prevalence | roc_auc | roc_auc_ci_lower | roc_auc_ci_upper | average_precision_ci_lower | average_precision_ci_upper | average_precision |\n| ------- | ------------- | ---- | ---------------- | ------------ | --------------------- | ------- | ---------------- | ---------------- | -------------------------- | -------------------------- | ----------------- |\n| primary | SIFT          | 1546 | 342              | 1204         | 0.2212                | 0.8587  | 0.8430           | 0.8737           | 0.4864                     | 0.5405                     | 0.5133            |\n| primary | PolyPhen-2    | 1546 | 342              | 1204         | 0.2212                | 0.7499  | 0.7204           | 0.7792           | 0.4457                     | 0.5381                     | 0.4900            |\n| primary | AlphaMissense | 1546 | 342              | 1204         | 0.2212                | 0.9525  | 0.9365           | 0.9672           | 0.8851                     | 0.9349                     | 0.9101            |\n| primary | CADD PHRED    | 1546 | 342              | 1204         | 0.2212                | 0.9507  | 0.9382           | 0.9614           | 0.8014                     | 0.8735                     | 0.8379            |\n| primary | REVEL         | 1546 | 342              | 1204         | 0.2212                | 0.9279  | 0.9100           | 0.9463           | 0.8101                     | 0.8822                     | 0.8472            |

### High-confidence benchmark

| cohort          | predictor     | n   | pathogenic_count | benign_count | pathogenic_prevalence | roc_auc | roc_auc_ci_lower | roc_auc_ci_upper | average_precision_ci_lower | average_precision_ci_upper | average_precision |\n| --------------- | ------------- | --- | ---------------- | ------------ | --------------------- | ------- | ---------------- | ---------------- | -------------------------- | -------------------------- | ----------------- |\n| high_confidence | SIFT          | 656 | 252              | 404          | 0.3841                | 0.8736  | 0.8505           | 0.8961           | 0.6872                     | 0.7625                     | 0.7227            |\n| high_confidence | PolyPhen-2    | 656 | 252              | 404          | 0.3841                | 0.7612  | 0.7238           | 0.7972           | 0.6162                     | 0.7180                     | 0.6676            |\n| high_confidence | AlphaMissense | 656 | 252              | 404          | 0.3841                | 0.9569  | 0.9377           | 0.9732           | 0.9327                     | 0.9701                     | 0.9532            |\n| high_confidence | CADD PHRED    | 656 | 252              | 404          | 0.3841                | 0.9478  | 0.9308           | 0.9637           | 0.8508                     | 0.9284                     | 0.8890            |\n| high_confidence | REVEL         | 656 | 252              | 404          | 0.3841                | 0.9268  | 0.9036           | 0.9487           | 0.8781                     | 0.9341                     | 0.9071            |

### Gene-specific performance

| cohort          | gene  | predictor     | n   | pathogenic_count | benign_count | roc_auc | average_precision |\n| --------------- | ----- | ------------- | --- | ---------------- | ------------ | ------- | ----------------- |\n| primary         | BRCA1 | SIFT          | 689 | 224              | 465          | 0.8666  | 0.6572            |\n| primary         | BRCA1 | PolyPhen-2    | 689 | 224              | 465          | 0.6632  | 0.4193            |\n| primary         | BRCA1 | AlphaMissense | 689 | 224              | 465          | 0.9518  | 0.9382            |\n| primary         | BRCA1 | CADD PHRED    | 689 | 224              | 465          | 0.9492  | 0.8964            |\n| primary         | BRCA1 | REVEL         | 689 | 224              | 465          | 0.9273  | 0.8790            |\n| primary         | BRCA2 | SIFT          | 857 | 118              | 739          | 0.8575  | 0.3669            |\n| primary         | BRCA2 | PolyPhen-2    | 857 | 118              | 739          | 0.9145  | 0.6983            |\n| primary         | BRCA2 | AlphaMissense | 857 | 118              | 739          | 0.9505  | 0.8498            |\n| primary         | BRCA2 | CADD PHRED    | 857 | 118              | 739          | 0.9630  | 0.7989            |\n| primary         | BRCA2 | REVEL         | 857 | 118              | 739          | 0.9184  | 0.8002            |\n| high_confidence | BRCA1 | SIFT          | 347 | 158              | 189          | 0.8817  | 0.7917            |\n| high_confidence | BRCA1 | PolyPhen-2    | 347 | 158              | 189          | 0.6871  | 0.5679            |\n| high_confidence | BRCA1 | AlphaMissense | 347 | 158              | 189          | 0.9608  | 0.9649            |\n| high_confidence | BRCA1 | CADD PHRED    | 347 | 158              | 189          | 0.9465  | 0.9173            |\n| high_confidence | BRCA1 | REVEL         | 347 | 158              | 189          | 0.9354  | 0.9310            |\n| high_confidence | BRCA2 | SIFT          | 309 | 94               | 215          | 0.8689  | 0.6336            |\n| high_confidence | BRCA2 | PolyPhen-2    | 309 | 94               | 215          | 0.9063  | 0.8312            |\n| high_confidence | BRCA2 | AlphaMissense | 309 | 94               | 215          | 0.9492  | 0.9297            |\n| high_confidence | BRCA2 | CADD PHRED    | 309 | 94               | 215          | 0.9564  | 0.8631            |\n| high_confidence | BRCA2 | REVEL         | 309 | 94               | 215          | 0.9225  | 0.8746            |

## Predictor correlations

Spearman correlations were calculated across all 16,193 missense variants using direction-aligned scores.

| predictor_1   | predictor_2   | spearman_rho | p_value | n     |\n| ------------- | ------------- | ------------ | ------- | ----- |\n| SIFT          | PolyPhen-2    | 0.5753       | 0.0000  | 16193 |\n| SIFT          | AlphaMissense | 0.6349       | 0.0000  | 16193 |\n| SIFT          | CADD PHRED    | 0.7928       | 0.0000  | 16193 |\n| SIFT          | REVEL         | 0.3950       | 0.0000  | 16193 |\n| PolyPhen-2    | AlphaMissense | 0.4951       | 0.0000  | 16193 |\n| PolyPhen-2    | CADD PHRED    | 0.6642       | 0.0000  | 16193 |\n| PolyPhen-2    | REVEL         | 0.3411       | 0.0000  | 16193 |\n| AlphaMissense | CADD PHRED    | 0.6656       | 0.0000  | 16193 |\n| AlphaMissense | REVEL         | 0.4568       | 0.0000  | 16193 |\n| CADD PHRED    | REVEL         | 0.5167       | 0.0000  | 16193 |

## Predictor ranking disagreement

Disagreement is quantified as the range of within-predictor percentile ranks for each variant. It describes relative ranking spread; it is not a measure of binary classification agreement or biological truth.

- Mean rank range: 0.4302
- Median rank range: 0.4230

The 100 variants with the greatest ranking spread are listed in `results/tables/top_100_predictor_rank_disagreements.csv`.

## Limitations

- ClinVar is not an independent experimental ground truth.
- The primary benchmark represents only the subset with unambiguous benign or pathogenic classifications.
- Conflicting classifications and VUS records are excluded from primary binary performance evaluation.
- The high-confidence subset is a sensitivity analysis, not a replacement for the primary benchmark.
- Predictor training data and clinical evidence may overlap. Potential evidence dependence and data leakage must be considered when interpreting performance.
- Average precision is prevalence-dependent.
- Correlation and ranking disagreement do not establish which predictor is correct for an individual variant.

## Reproducibility

- Random seed: 42
- Bootstrap iterations: 2000
- Bootstrap: stratified percentile bootstrap by benchmark class
- Input SHA256: `7c8eea9224091f9fb68a449e147753fcf31a35bbd9d44fcdccd08dab78de5797`
- Generated: 2026-11-08T18:12:19.033751+00:00


## Additional statistical reporting notes

The correlation analysis was conducted across the complete missense cohort
(n=16,193). Some Spearman correlation p-values were returned as 0.0 due to
floating-point underflow. These values must not be interpreted as literal
zero probabilities. Correlation magnitude is reported as Spearman's rho.

The all-missense cohort contains 16,193 variants, of which 1,546 have
unambiguous binary ClinVar labels used in the primary benchmark. The
remaining 14,647 variants are not part of that binary benchmark.

Gene-specific performance metrics are descriptive subgroup analyses.
Differences between BRCA1 and BRCA2 should not be interpreted as formal
evidence of gene-specific superiority without additional statistical
testing and consideration of class prevalence and cohort composition.

The bootstrap confidence intervals are calculated independently for each
predictor. Overlap or non-overlap of these intervals is not a substitute
for a paired statistical comparison between predictors.
