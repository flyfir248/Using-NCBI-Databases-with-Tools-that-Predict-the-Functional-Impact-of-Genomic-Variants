import json
import logging
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    ROOT / "data/processed/"
    "brca1_brca2_primary_missense_cohort.csv"
)

OUTPUT_DIR = ROOT / "data/processed"
METADATA_DIR = ROOT / "data/metadata"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

logger = logging.getLogger(__name__)

PATHOGENIC_LABELS = {
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
}

BENIGN_LABELS = {
    "Benign",
    "Likely benign",
    "Benign/Likely benign",
}

HIGH_CONFIDENCE_STATUSES = {
    "reviewed by expert panel",
    "criteria provided, multiple submitters, no conflicts",
}

# Explicit mapping of review status to evidence category.
REVIEW_STATUS_GROUPS = {
    "reviewed by expert panel": "expert_panel",
    "criteria provided, multiple submitters, no conflicts":
        "multiple_submitters_no_conflicts",
    "criteria provided, single submitter": "single_submitter",
    "criteria provided, conflicting classifications":
        "conflicting_submissions",
    "no assertion criteria provided": "no_assertion_criteria",
    "no classification provided": "no_classification_provided",
    "-": "missing_review_status",
}


def main():

    logger.info("Loading primary missense cohort...")

    df = pd.read_csv(INPUT_FILE, low_memory=False)

    required_columns = {
        "#AlleleID",
        "GeneSymbol",
        "ClinicalSignificance",
        "ReviewStatus",
    }

    missing = required_columns - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}"
        )

    if not df["#AlleleID"].is_unique:
        raise ValueError("Duplicate ClinVar AlleleIDs detected.")

    if len(df) != 16193:
        raise ValueError(
            f"Unexpected missense cohort size: {len(df)}"
        )

    # Preserve source fields; add independent benchmark fields.
    df["benchmark_group"] = "excluded"
    df["clinvar_binary_label"] = pd.Series(
        pd.NA, index=df.index, dtype="Int64"
    )

    df["benchmark_eligible"] = False
    df["benchmark_exclusion_reason"] = ""

    # Normalize missing values for classification logic only.
    clinical = (
        df["ClinicalSignificance"]
        .fillna("<MISSING>")
        .astype(str)
        .str.strip()
    )

    review = (
        df["ReviewStatus"]
        .fillna("<MISSING>")
        .astype(str)
        .str.strip()
    )

    # Assign pathogenic labels.
    pathogenic_mask = clinical.isin(PATHOGENIC_LABELS)

    df.loc[pathogenic_mask, "benchmark_group"] = "pathogenic"
    df.loc[pathogenic_mask, "clinvar_binary_label"] = 1

    # Assign benign labels.
    benign_mask = clinical.isin(BENIGN_LABELS)

    df.loc[benign_mask, "benchmark_group"] = "benign"
    df.loc[benign_mask, "clinvar_binary_label"] = 0

    eligible_mask = pathogenic_mask | benign_mask

    df["benchmark_eligible"] = eligible_mask

    # Explicit exclusion reasons.
    df.loc[
        clinical == "Uncertain significance",
        "benchmark_exclusion_reason"
    ] = "uncertain_significance"

    df.loc[
        clinical == "Conflicting classifications of pathogenicity",
        "benchmark_exclusion_reason"
    ] = "conflicting_clinvar_classifications"

    df.loc[
        clinical == "-",
        "benchmark_exclusion_reason"
    ] = "missing_clinical_significance"

    df.loc[
        clinical == "not provided",
        "benchmark_exclusion_reason"
    ] = "clinical_significance_not_provided"

    known_values = (
        PATHOGENIC_LABELS
        | BENIGN_LABELS
        | {
            "Uncertain significance",
            "Conflicting classifications of pathogenicity",
            "-",
            "not provided",
        }
    )

    unknown_mask = ~clinical.isin(known_values)

    df.loc[
        unknown_mask,
        "benchmark_exclusion_reason"
    ] = "other_or_unrecognized_clinical_significance"

    # Review-status evidence categories.
    df["review_confidence_group"] = (
        review.map(REVIEW_STATUS_GROUPS)
        .fillna("other_review_status")
    )

    # High-confidence sensitivity benchmark.
    df["high_confidence_eligible"] = (
        eligible_mask
        & review.isin(HIGH_CONFIDENCE_STATUSES)
    )

    # Audit consistency.
    if int(pathogenic_mask.sum()) != 342:
        raise ValueError("Unexpected pathogenic label count.")

    if int(benign_mask.sum()) != 1204:
        raise ValueError("Unexpected benign label count.")

    if int(df["high_confidence_eligible"].sum()) != 656:
        raise ValueError("Unexpected high-confidence count.")

    if int(df["benchmark_eligible"].sum()) != 1546:
        raise ValueError("Unexpected primary benchmark count.")

    if df.loc[
        ~eligible_mask, "clinvar_binary_label"
    ].notna().any():
        raise ValueError(
            "Excluded variants unexpectedly have binary labels."
        )

    # Create benchmark subsets.
    primary = df[df["benchmark_eligible"]].copy()

    high_confidence = df[
        df["high_confidence_eligible"]
    ].copy()

    # Verify class counts.
    primary_counts = (
        primary["benchmark_group"].value_counts().to_dict()
    )

    high_confidence_counts = (
        high_confidence["benchmark_group"]
        .value_counts()
        .to_dict()
    )

    expected_primary = {
        "pathogenic": 342,
        "benign": 1204,
    }

    expected_high_confidence = {
        "benign": 404,
        "pathogenic": 252,
    }

    if primary_counts != expected_primary:
        raise ValueError(
            f"Unexpected primary counts: {primary_counts}"
        )

    if high_confidence_counts != expected_high_confidence:
        raise ValueError(
            "Unexpected high-confidence counts: "
            f"{high_confidence_counts}"
        )

    # Output files.
    annotated_path = (
        OUTPUT_DIR /
        "brca1_brca2_missense_benchmark_annotated.csv"
    )

    primary_path = (
        OUTPUT_DIR /
        "brca1_brca2_primary_benchmark.csv"
    )

    high_confidence_path = (
        OUTPUT_DIR /
        "brca1_brca2_high_confidence_benchmark.csv"
    )

    annotated_path.parent.mkdir(parents=True, exist_ok=True)
    METADATA_DIR.mkdir(parents=True, exist_ok=True)

    df.to_csv(annotated_path, index=False)
    primary.to_csv(primary_path, index=False)
    high_confidence.to_csv(high_confidence_path, index=False)

    # Detailed summaries.
    gene_summary = {}

    for gene in ["BRCA1", "BRCA2"]:

        gene_data = df[df["GeneSymbol"] == gene]
        gene_primary = primary[
            primary["GeneSymbol"] == gene
        ]
        gene_high = high_confidence[
            high_confidence["GeneSymbol"] == gene
        ]

        gene_summary[gene] = {
            "total_missense_variants": int(len(gene_data)),
            "primary_benchmark": {
                "pathogenic": int(
                    (gene_primary["benchmark_group"] == "pathogenic").sum()
                ),
                "benign": int(
                    (gene_primary["benchmark_group"] == "benign").sum()
                ),
            },
            "high_confidence_benchmark": {
                "pathogenic": int(
                    (gene_high["benchmark_group"] == "pathogenic").sum()
                ),
                "benign": int(
                    (gene_high["benchmark_group"] == "benign").sum()
                ),
            },
        }

    review_summary = (
        pd.crosstab(
            df["ReviewStatus"].fillna("<MISSING>"),
            df["benchmark_group"],
        )
        .to_dict()
    )

    exclusion_summary = (
        df.loc[~eligible_mask, "benchmark_exclusion_reason"]
        .value_counts()
        .to_dict()
    )

    report = {
        "input_file": str(INPUT_FILE),
        "total_missense_variants": int(len(df)),
        "primary_benchmark_size": int(len(primary)),
        "high_confidence_benchmark_size": int(
            len(high_confidence)
        ),
        "primary_class_counts": {
            str(k): int(v)
            for k, v in primary_counts.items()
        },
        "high_confidence_class_counts": {
            str(k): int(v)
            for k, v in high_confidence_counts.items()
        },
        "excluded_variants": int((~eligible_mask).sum()),
        "exclusion_reasons": {
            str(k): int(v)
            for k, v in exclusion_summary.items()
        },
        "gene_summary": gene_summary,
        "review_status_classification": {
            str(status): {
                str(group): int(count)
                for group, count in groups.items()
            }
            for status, groups in review_summary.items()
        },
        "high_confidence_review_statuses": sorted(
            HIGH_CONFIDENCE_STATUSES
        ),
        "label_policy": {
            "pathogenic": sorted(PATHOGENIC_LABELS),
            "benign": sorted(BENIGN_LABELS),
            "excluded": [
                "Uncertain significance",
                "Conflicting classifications of pathogenicity",
                "-",
                "not provided",
                "other or unrecognized values",
            ],
        },
        "validation": {
            "unique_allele_ids": bool(df["#AlleleID"].is_unique),
            "primary_count_correct": len(primary) == 1546,
            "high_confidence_count_correct": len(high_confidence) == 656,
            "excluded_records_have_no_binary_label": bool(
                df.loc[~eligible_mask, "clinvar_binary_label"]
                .isna().all()
            ),
            "all_primary_labels_assigned": bool(
                primary["clinvar_binary_label"].notna().all()
            ),
        },
    }

    report["all_checks_passed"] = all(
        report["validation"].values()
    )

    report_path = (
        METADATA_DIR / "clinvar_benchmark_report.json"
    )

    report_path.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )

    logger.info("=" * 60)
    logger.info("CLINVAR BENCHMARK CONSTRUCTION COMPLETED")
    logger.info("=" * 60)

    logger.info("Total missense variants: %s", len(df))
    logger.info("Primary benchmark: %s", len(primary))
    logger.info("High-confidence benchmark: %s", len(high_confidence))
    logger.info("Excluded: %s", int((~eligible_mask).sum()))

    logger.info("Primary class counts: %s", primary_counts)
    logger.info("High-confidence class counts: %s", high_confidence_counts)

    logger.info("Gene summary: %s", gene_summary)

    logger.info("Annotated dataset: %s", annotated_path)
    logger.info("Primary benchmark: %s", primary_path)
    logger.info("High-confidence benchmark: %s", high_confidence_path)
    logger.info("Report: %s", report_path)

    if not report["all_checks_passed"]:
        raise SystemExit("BENCHMARK VALIDATION FAILED.")

    logger.info("ALL BENCHMARK VALIDATIONS PASSED.")


if __name__ == "__main__":
    main()
