from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.stats import spearmanr
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)


ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    ROOT / "data/processed/"
    "brca1_brca2_all_predictors_integrated.csv"
)

TABLE_DIR = ROOT / "results/tables"
FIGURE_DIR = ROOT / "results/figures"
REPORT_DIR = ROOT / "results/reports"

JSON_REPORT = ROOT / "data/metadata/predictor_evaluation_report.json"
MARKDOWN_REPORT = REPORT_DIR / "predictor_evaluation_report.md"

RANDOM_SEED = 42
BOOTSTRAP_ITERATIONS = 2000

PREDICTORS = {
    "SIFT": {
        "column": "sift_score",
        "direction": "lower_is_more_damaging",
    },
    "PolyPhen-2": {
        "column": "polyphen_score",
        "direction": "higher_is_more_damaging",
    },
    "AlphaMissense": {
        "column": "alphamissense_pathogenicity",
        "direction": "higher_is_more_pathogenic",
    },
    "CADD PHRED": {
        "column": "cadd_phred",
        "direction": "higher_is_more_damaging",
    },
    "REVEL": {
        "column": "revel",
        "direction": "higher_is_more_pathogenic",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "_No rows available._"

    columns = [str(column) for column in frame.columns]

    def format_value(value):
        if pd.isna(value):
            return ""
        if isinstance(value, (float, np.floating)):
            return f"{value:.4f}"
        return str(value).replace("|", "\\|").replace("\\n", " ")

    rows = [
        [format_value(value) for value in row]
        for row in frame.itertuples(index=False, name=None)
    ]

    widths = [
        max(len(columns[i]), *(len(row[i]) for row in rows))
        for i in range(len(columns))
    ]

    header = "| " + " | ".join(
        columns[i].ljust(widths[i]) for i in range(len(columns))
    ) + " |"

    separator = "| " + " | ".join(
        "-" * widths[i] for i in range(len(columns))
    ) + " |"

    body = [
        "| " + " | ".join(
            row[i].ljust(widths[i]) for i in range(len(columns))
        ) + " |"
        for row in rows
    ]

    return "\\n".join([header, separator, *body])


def parse_boolean(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)

    normalized = series.astype("string").str.strip().str.lower()
    return normalized.isin(["true", "1", "yes", "y"])


def prepare_scores(df: pd.DataFrame) -> pd.DataFrame:
    scores = pd.DataFrame(index=df.index)

    for name, config in PREDICTORS.items():
        values = pd.to_numeric(
            df[config["column"]], errors="coerce"
        )

        if values.isna().any():
            raise ValueError(
                f"{name} contains missing or non-numeric scores."
            )

        if not np.isfinite(values.to_numpy(dtype=float)).all():
            raise ValueError(f"{name} contains non-finite scores.")

        if config["direction"] == "lower_is_more_damaging":
            values = 1.0 - values

        scores[name] = values.astype(float)

    return scores


def bootstrap_metrics(
    y: np.ndarray,
    scores: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
) -> dict:
    positive = np.flatnonzero(y == 1)
    negative = np.flatnonzero(y == 0)

    if len(positive) == 0 or len(negative) == 0:
        raise ValueError("Both benchmark classes are required.")

    auc_values = []
    ap_values = []

    for _ in range(iterations):
        sampled_positive = rng.choice(
            positive, size=len(positive), replace=True
        )
        sampled_negative = rng.choice(
            negative, size=len(negative), replace=True
        )

        indices = np.concatenate(
            [sampled_positive, sampled_negative]
        )
        rng.shuffle(indices)

        sampled_y = y[indices]
        sampled_scores = scores[indices]

        auc_values.append(
            roc_auc_score(sampled_y, sampled_scores)
        )
        ap_values.append(
            average_precision_score(sampled_y, sampled_scores)
        )

    return {
        "roc_auc_ci_lower": float(np.percentile(auc_values, 2.5)),
        "roc_auc_ci_upper": float(np.percentile(auc_values, 97.5)),
        "average_precision_ci_lower": float(
            np.percentile(ap_values, 2.5)
        ),
        "average_precision_ci_upper": float(
            np.percentile(ap_values, 97.5)
        ),
    }


def evaluate_cohort(
    df: pd.DataFrame,
    scores: pd.DataFrame,
    cohort_name: str,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, dict]:
    y = pd.to_numeric(
        df["clinvar_binary_label"], errors="coerce"
    ).to_numpy(dtype=float)

    if not np.isin(y, [0, 1]).all():
        raise ValueError(
            f"{cohort_name} contains invalid benchmark labels."
        )

    y = y.astype(int)

    if len(np.unique(y)) != 2:
        raise ValueError(
            f"{cohort_name} does not contain both benchmark classes."
        )

    rows = []
    curves = {}

    for name in PREDICTORS:
        values = scores.loc[df.index, name].to_numpy(dtype=float)

        auc_value = roc_auc_score(y, values)
        ap_value = average_precision_score(y, values)

        confidence = bootstrap_metrics(
            y,
            values,
            BOOTSTRAP_ITERATIONS,
            rng,
        )

        rows.append({
            "cohort": cohort_name,
            "predictor": name,
            "n": len(y),
            "pathogenic_count": int((y == 1).sum()),
            "benign_count": int((y == 0).sum()),
            "pathogenic_prevalence": float(y.mean()),
            "roc_auc": float(auc_value),
            **confidence,
            "average_precision": float(ap_value),
        })

        fpr, tpr, _ = roc_curve(y, values)
        precision, recall, _ = precision_recall_curve(y, values)

        curves[name] = {
            "fpr": fpr,
            "tpr": tpr,
            "precision": precision,
            "recall": recall,
            "roc_auc": auc_value,
            "average_precision": ap_value,
        }

    return pd.DataFrame(rows), curves


def save_roc_plot(curves: dict, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 7))

    for name, curve in curves.items():
        ax.plot(
            curve["fpr"],
            curve["tpr"],
            label=f"{name} (AUC={curve['roc_auc']:.3f})",
        )

    ax.plot([0, 1], [0, 1], linestyle="--", label="Random baseline")
    ax.set(
        xlabel="False positive rate",
        ylabel="True positive rate",
        title=title,
        xlim=(0, 1),
        ylim=(0, 1.02),
    )
    ax.legend(loc="lower right")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def save_pr_plot(curves: dict, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 7))

    for name, curve in curves.items():
        ax.plot(
            curve["recall"],
            curve["precision"],
            label=f"{name} (AP={curve['average_precision']:.3f})",
        )

    ax.set(
        xlabel="Recall",
        ylabel="Precision",
        title=title,
        xlim=(0, 1),
        ylim=(0, 1.02),
    )
    ax.legend(loc="best")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main() -> None:
    for directory in [TABLE_DIR, FIGURE_DIR, REPORT_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    print("Loading five-predictor integrated dataset...")
    df = pd.read_csv(INPUT_FILE, low_memory=False)

    required = [
        "variant_key",
        "GeneSymbol",
        "clinvar_binary_label",
        "benchmark_eligible",
        "high_confidence_eligible",
        *[config["column"] for config in PREDICTORS.values()],
    ]

    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if df["variant_key"].duplicated().any():
        raise ValueError("Duplicate variant keys detected.")

    if len(df) != 16193:
        raise ValueError(
            f"Expected 16,193 variants; found {len(df)}."
        )

    scores = prepare_scores(df)

    primary_mask = (
        parse_boolean(df["benchmark_eligible"])
        & df["clinvar_binary_label"].notna()
    )

    high_mask = (
        parse_boolean(df["high_confidence_eligible"])
        & df["clinvar_binary_label"].notna()
    )

    primary = df.loc[primary_mask].copy()
    high = df.loc[high_mask].copy()

    print("Full missense cohort:", len(df))
    print("Primary benchmark:", len(primary))
    print("High-confidence benchmark:", len(high))

    rng = np.random.default_rng(RANDOM_SEED)

    primary_metrics, primary_curves = evaluate_cohort(
        primary, scores, "primary", rng
    )

    high_metrics, high_curves = evaluate_cohort(
        high, scores, "high_confidence", rng
    )

    primary_metrics.to_csv(
        TABLE_DIR / "predictor_performance_primary.csv",
        index=False,
    )

    high_metrics.to_csv(
        TABLE_DIR / "predictor_performance_high_confidence.csv",
        index=False,
    )

    # Gene-specific performance, without bootstrap CIs.
    gene_rows = []

    for cohort_name, cohort in [
        ("primary", primary),
        ("high_confidence", high),
    ]:
        for gene in ["BRCA1", "BRCA2"]:
            subset = cohort[cohort["GeneSymbol"] == gene]

            if subset["clinvar_binary_label"].nunique() != 2:
                continue

            y = subset["clinvar_binary_label"].astype(int).to_numpy()

            for name in PREDICTORS:
                values = scores.loc[subset.index, name].to_numpy()

                gene_rows.append({
                    "cohort": cohort_name,
                    "gene": gene,
                    "predictor": name,
                    "n": len(subset),
                    "pathogenic_count": int((y == 1).sum()),
                    "benign_count": int((y == 0).sum()),
                    "roc_auc": float(roc_auc_score(y, values)),
                    "average_precision": float(
                        average_precision_score(y, values)
                    ),
                })

    gene_metrics = pd.DataFrame(gene_rows)
    gene_metrics.to_csv(
        TABLE_DIR / "predictor_performance_by_gene.csv",
        index=False,
    )

    # Spearman correlations use direction-aligned scores.
    correlation = scores.corr(method="spearman")
    correlation.to_csv(
        TABLE_DIR / "predictor_spearman_correlations.csv"
    )

    correlation_rows = []
    for i, first in enumerate(PREDICTORS):
        for second in list(PREDICTORS)[i + 1:]:
            result = spearmanr(scores[first], scores[second])

            correlation_rows.append({
                "predictor_1": first,
                "predictor_2": second,
                "spearman_rho": float(result.statistic),
                "p_value": float(result.pvalue),
                "n": len(scores),
            })

    correlation_pairs = pd.DataFrame(correlation_rows)
    correlation_pairs.to_csv(
        TABLE_DIR / "predictor_correlation_pairs.csv",
        index=False,
    )

    # Relative-rank disagreement across predictors.
    # This measures ranking spread, not binary classification disagreement.
    percentile_ranks = scores.rank(pct=True)

    disagreement = df[
        ["variant_key", "GeneSymbol", "ClinicalSignificance"]
    ].copy()

    disagreement["rank_percentile_range"] = (
        percentile_ranks.max(axis=1)
        - percentile_ranks.min(axis=1)
    )

    disagreement["rank_percentile_sd"] = percentile_ranks.std(axis=1)

    disagreement = disagreement.sort_values(
        "rank_percentile_range",
        ascending=False,
    )

    disagreement.to_csv(
        TABLE_DIR / "predictor_rank_disagreement_all_variants.csv",
        index=False,
    )

    disagreement.head(100).to_csv(
        TABLE_DIR / "top_100_predictor_rank_disagreements.csv",
        index=False,
    )

    # Cohort accounting.
    cohort_summary = pd.DataFrame([
        {
            "cohort": "all_missense",
            "n": len(df),
            "pathogenic": int(
                (df["clinvar_binary_label"] == 1).sum()
            ),
            "benign": int(
                (df["clinvar_binary_label"] == 0).sum()
            ),
        },
        {
            "cohort": "primary",
            "n": len(primary),
            "pathogenic": int(
                (primary["clinvar_binary_label"] == 1).sum()
            ),
            "benign": int(
                (primary["clinvar_binary_label"] == 0).sum()
            ),
        },
        {
            "cohort": "high_confidence",
            "n": len(high),
            "pathogenic": int(
                (high["clinvar_binary_label"] == 1).sum()
            ),
            "benign": int(
                (high["clinvar_binary_label"] == 0).sum()
            ),
        },
    ])

    cohort_summary.to_csv(
        TABLE_DIR / "benchmark_cohort_summary.csv",
        index=False,
    )

    # Figures.
    save_roc_plot(
        primary_curves,
        FIGURE_DIR / "roc_curves_primary.png",
        "Five-predictor ROC curves: primary ClinVar benchmark",
    )

    save_pr_plot(
        primary_curves,
        FIGURE_DIR / "pr_curves_primary.png",
        "Five-predictor precision-recall curves: primary benchmark",
    )

    save_roc_plot(
        high_curves,
        FIGURE_DIR / "roc_curves_high_confidence.png",
        "Five-predictor ROC curves: high-confidence benchmark",
    )

    save_pr_plot(
        high_curves,
        FIGURE_DIR / "pr_curves_high_confidence.png",
        "Five-predictor precision-recall curves: high-confidence benchmark",
    )

    fig, ax = plt.subplots(figsize=(8, 7))
    image = ax.imshow(correlation.to_numpy(), vmin=-1, vmax=1)
    ax.set_xticks(range(len(correlation.columns)))
    ax.set_xticklabels(correlation.columns, rotation=35, ha="right")
    ax.set_yticks(range(len(correlation.index)))
    ax.set_yticklabels(correlation.index)
    ax.set_title("Spearman correlation between predictors")
    fig.colorbar(image, ax=ax, label="Spearman correlation")
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "predictor_correlation_heatmap.png",
        dpi=300,
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.boxplot(
        [scores[name].to_numpy() for name in PREDICTORS],
        tick_labels=list(PREDICTORS.keys()),
        showfliers=False,
    )
    ax.set_ylabel("Direction-aligned predictor score")
    ax.set_title("Predictor score distributions across missense cohort")
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "predictor_score_distributions.png",
        dpi=300,
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(
        disagreement["rank_percentile_range"],
        bins=40,
    )
    ax.set(
        xlabel="Range of within-predictor percentile ranks",
        ylabel="Number of variants",
        title="Relative predictor ranking disagreement",
    )
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "predictor_rank_disagreement.png",
        dpi=300,
    )
    plt.close(fig)

    # Consolidated JSON report.
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(INPUT_FILE.relative_to(ROOT)),
        "input_sha256": sha256(INPUT_FILE),
        "total_missense_variants": len(df),
        "primary_benchmark": {
            "n": len(primary),
            "pathogenic": int(
                (primary["clinvar_binary_label"] == 1).sum()
            ),
            "benign": int(
                (primary["clinvar_binary_label"] == 0).sum()
            ),
            "metrics": primary_metrics.to_dict(orient="records"),
        },
        "high_confidence_benchmark": {
            "n": len(high),
            "pathogenic": int(
                (high["clinvar_binary_label"] == 1).sum()
            ),
            "benign": int(
                (high["clinvar_binary_label"] == 0).sum()
            ),
            "metrics": high_metrics.to_dict(orient="records"),
        },
        "predictors": PREDICTORS,
        "correlation_method": "Spearman rank correlation",
        "correlation_score_direction": (
            "All scores direction-aligned so higher means more damaging."
        ),
        "disagreement_method": (
            "Per-variant range and standard deviation of predictor "
            "percentile ranks across the full missense cohort. "
            "This is relative ranking disagreement, not a binary call."
        ),
        "bootstrap": {
            "iterations": BOOTSTRAP_ITERATIONS,
            "seed": RANDOM_SEED,
            "method": "Stratified percentile bootstrap by benchmark class",
            "confidence_interval": "95%",
        },
        "limitations": [
            "ClinVar classifications are an evidence-based benchmark, "
            "not an independent experimental ground truth.",
            "The primary benchmark contains only unambiguous benign-side "
            "and pathogenic-side ClinVar classifications.",
            "The high-confidence benchmark is a sensitivity analysis "
            "restricted by ClinVar review status.",
            "Predictor training data may overlap with ClinVar or related "
            "clinical evidence; apparent performance may therefore be "
            "affected by data leakage or evidence dependence.",
            "Average precision depends on pathogenic prevalence and "
            "should be interpreted alongside ROC-AUC and class counts.",
            "The relative-rank disagreement analysis does not establish "
            "which predictor is correct for an individual variant.",
        ],
        "outputs": {
            "primary_metrics": "results/tables/predictor_performance_primary.csv",
            "high_confidence_metrics": (
                "results/tables/predictor_performance_high_confidence.csv"
            ),
            "gene_metrics": (
                "results/tables/predictor_performance_by_gene.csv"
            ),
            "correlations": (
                "results/tables/predictor_correlation_pairs.csv"
            ),
            "rank_disagreement": (
                "results/tables/predictor_rank_disagreement_all_variants.csv"
            ),
            "report": "results/reports/predictor_evaluation_report.md",
        },
    }

    JSON_REPORT.parent.mkdir(parents=True, exist_ok=True)
    JSON_REPORT.write_text(
        json.dumps(report, indent=2, default=str),
        encoding="utf-8",
    )

    # Human-readable report.
    lines = [
        "# Integrative Functional Impact Assessment",
        "",
        "## Study cohort",
        "",
        f"- Total missense variants: **{len(df):,}**",
        f"- Primary ClinVar benchmark: **{len(primary):,}**",
        f"- High-confidence benchmark: **{len(high):,}**",
        "",
        "The primary benchmark includes unambiguous benign-side and "
        "pathogenic-side ClinVar classifications. The high-confidence "
        "benchmark is a review-status-restricted sensitivity analysis.",
        "",
        "## Predictor evaluation",
        "",
        "Higher direction-aligned scores indicate greater predicted "
        "deleteriousness or pathogenicity. SIFT was transformed as "
        "`1 - SIFT score`; the other four scores were retained in "
        "their original direction.",
        "",
        "ROC-AUC and average precision are reported with stratified "
        "bootstrap 95% confidence intervals.",
        "",
        "### Primary benchmark",
        "",
        markdown_table(primary_metrics),
        "",
        "### High-confidence benchmark",
        "",
        markdown_table(high_metrics),
        "",
        "### Gene-specific performance",
        "",
        markdown_table(gene_metrics)
        if not gene_metrics.empty else "No eligible gene subgroups.",
        "",
        "## Predictor correlations",
        "",
        "Spearman correlations were calculated across all 16,193 "
        "missense variants using direction-aligned scores.",
        "",
        markdown_table(correlation_pairs),
        "",
        "## Predictor ranking disagreement",
        "",
        "Disagreement is quantified as the range of within-predictor "
        "percentile ranks for each variant. It describes relative "
        "ranking spread; it is not a measure of binary classification "
        "agreement or biological truth.",
        "",
        f"- Mean rank range: "
        f"{disagreement['rank_percentile_range'].mean():.4f}",
        f"- Median rank range: "
        f"{disagreement['rank_percentile_range'].median():.4f}",
        "",
        "The 100 variants with the greatest ranking spread are listed "
        "in `results/tables/top_100_predictor_rank_disagreements.csv`.",
        "",
        "## Limitations",
        "",
        "- ClinVar is not an independent experimental ground truth.",
        "- The primary benchmark represents only the subset with "
        "unambiguous benign or pathogenic classifications.",
        "- Conflicting classifications and VUS records are excluded "
        "from primary binary performance evaluation.",
        "- The high-confidence subset is a sensitivity analysis, "
        "not a replacement for the primary benchmark.",
        "- Predictor training data and clinical evidence may overlap. "
        "Potential evidence dependence and data leakage must be "
        "considered when interpreting performance.",
        "- Average precision is prevalence-dependent.",
        "- Correlation and ranking disagreement do not establish "
        "which predictor is correct for an individual variant.",
        "",
        "## Reproducibility",
        "",
        f"- Random seed: {RANDOM_SEED}",
        f"- Bootstrap iterations: {BOOTSTRAP_ITERATIONS}",
        "- Bootstrap: stratified percentile bootstrap by benchmark class",
        f"- Input SHA256: `{sha256(INPUT_FILE)}`",
        f"- Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
    ]

    MARKDOWN_REPORT.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print("\nPREDICTOR EVALUATION COMPLETED")
    print("\nPrimary benchmark:")
    print(primary_metrics.to_string(index=False))

    print("\nHigh-confidence benchmark:")
    print(high_metrics.to_string(index=False))

    print("\nOutputs:")
    print("Tables:", TABLE_DIR)
    print("Figures:", FIGURE_DIR)
    print("Markdown report:", MARKDOWN_REPORT)
    print("JSON report:", JSON_REPORT)


if __name__ == "__main__":
    main()
