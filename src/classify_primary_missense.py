
import json
import logging
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    ROOT / "data/processed/vep_primary_transcripts.csv"
)

ALL_OUTPUT = (
    ROOT / "data/processed/vep_primary_missense_classification.csv"
)

MISSENSE_OUTPUT = (
    ROOT / "data/processed/brca1_brca2_primary_missense_cohort.csv"
)

NON_MISSENSE_OUTPUT = (
    ROOT / "data/processed/brca1_brca2_primary_non_missense.csv"
)

REPORT_OUTPUT = (
    ROOT / "data/metadata/vep_primary_missense_report.json"
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

logger = logging.getLogger(__name__)


def parse_consequence_terms(value):
    """Convert VEP pipe-separated consequence terms into a set."""

    if pd.isna(value):
        return set()

    return {
        term.strip()
        for term in str(value).split("|")
        if term.strip()
    }


def main():

    logger.info("Loading primary transcript dataset...")

    data = pd.read_csv(
        INPUT_FILE,
        low_memory=False,
    )

    required_columns = {
        "#AlleleID",
        "GeneSymbol",
        "transcript_selection_method",
        "transcript_id",
        "gene_symbol",
        "biotype",
        "consequence_terms",
    }

    missing_columns = required_columns - set(data.columns)

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {sorted(missing_columns)}"
        )

    if data["#AlleleID"].duplicated().any():
        raise ValueError("Duplicate ClinVar AlleleIDs detected.")

    if len(data) != 28302:
        logger.warning(
            "Expected 28302 variants; observed %s.",
            len(data),
        )

    # Confirm selected transcript belongs to the ClinVar target gene.
    selected = data[
        data["transcript_selection_method"].isin(
            ["MANE_Select", "canonical_fallback"]
        )
    ].copy()

    selected["target_gene_match"] = (
        selected["GeneSymbol"].astype("string").str.upper()
        == selected["gene_symbol"].astype("string").str.upper()
    )

    if not selected["target_gene_match"].all():
        raise ValueError(
            "A selected transcript does not match its ClinVar target gene."
        )

    # Parse consequence terms from the selected transcript only.
    data["primary_consequence_term_set"] = (
        data["consequence_terms"].apply(parse_consequence_terms)
    )

    data["is_primary_missense"] = (
        data["primary_consequence_term_set"].apply(
            lambda terms: "missense_variant" in terms
        )
    )

    # Classify every original variant without silently dropping records.
    def classify(row):

        method = row["transcript_selection_method"]

        if method == "manual_review":
            return "manual_review_no_selected_target_gene_transcript"

        if not row["primary_consequence_term_set"]:
            return "manual_review_missing_consequence"

        if row["is_primary_missense"]:
            return "missense"

        return "non_missense"

    data["missense_classification"] = data.apply(
        classify,
        axis=1,
    )

    # Check selected transcripts are protein-coding for missense records.
    missense = data[
        data["missense_classification"] == "missense"
    ].copy()

    invalid_biotype = (
        missense["biotype"].astype("string").str.lower()
        != "protein_coding"
    )

    if invalid_biotype.any():
        raise ValueError(
            "Missense records include non-protein-coding transcripts."
        )

    non_missense = data[
        data["missense_classification"] == "non_missense"
    ].copy()

    # Remove the internal Python set column before CSV serialization.
    output_data = data.drop(
        columns=["primary_consequence_term_set"]
    )

    output_missense = missense.drop(
        columns=["primary_consequence_term_set"]
    )

    output_non_missense = non_missense.drop(
        columns=["primary_consequence_term_set"]
    )

    # Save all variants and the separate analysis cohorts.
    output_data.to_csv(ALL_OUTPUT, index=False)
    output_missense.to_csv(MISSENSE_OUTPUT, index=False)
    output_non_missense.to_csv(NON_MISSENSE_OUTPUT, index=False)

    # Build transparent summary statistics.
    classification_counts = (
        data["missense_classification"]
        .value_counts()
        .to_dict()
    )

    gene_summary = {}

    for gene, group in data.groupby("GeneSymbol"):

        gene_summary[str(gene)] = {
            "total_candidates": int(len(group)),
            "missense": int(
                (group["missense_classification"] == "missense").sum()
            ),
            "non_missense": int(
                (group["missense_classification"] == "non_missense").sum()
            ),
            "manual_review": int(
                group["missense_classification"].str.startswith(
                    "manual_review"
                ).sum()
            ),
        }

    consequence_distribution = (
        data.loc[
            data["missense_classification"] == "non_missense",
            "consequence_terms",
        ]
        .value_counts()
        .head(30)
        .to_dict()
    )

    report = {
        "input_variant_count": int(len(data)),
        "missense_cohort_count": int(len(missense)),
        "non_missense_count": int(len(non_missense)),
        "manual_review_count": int(
            data["missense_classification"]
            .str.startswith("manual_review")
            .sum()
        ),
        "classification_distribution": {
            str(key): int(value)
            for key, value in classification_counts.items()
        },
        "gene_summary": gene_summary,
        "top_non_missense_consequences": {
            str(key): int(value)
            for key, value in consequence_distribution.items()
        },
        "classification_policy": [
            "Use selected target-gene MANE Select transcript.",
            "Use target-gene canonical transcript only as fallback.",
            "Missense requires missense_variant in selected transcript consequence_terms.",
            "Do not use variant-level most_severe_consequence for cohort selection.",
            "Retain manual-review and non-missense variants in the full classification dataset.",
            "Exclude manual-review cases from the missense cohort.",
        ],
        "outputs": {
            "all_classifications": str(ALL_OUTPUT),
            "missense_cohort": str(MISSENSE_OUTPUT),
            "non_missense": str(NON_MISSENSE_OUTPUT),
        },
    }

    REPORT_OUTPUT.write_text(
        json.dumps(report, indent=4),
        encoding="utf-8",
    )

    # Final integrity checks.
    if len(output_data) != len(missense) + len(non_missense) + report["manual_review_count"]:
        raise ValueError("Classification counts do not reconcile.")

    logger.info("==============================================")
    logger.info("PRIMARY MISSENSE CLASSIFICATION COMPLETED")
    logger.info("==============================================")
    logger.info("Total candidates: %s", len(data))
    logger.info("Missense cohort: %s", len(missense))
    logger.info("Non-missense: %s", len(non_missense))
    logger.info("Manual review: %s", report["manual_review_count"])
    logger.info("Gene summary: %s", gene_summary)
    logger.info("Missense output: %s", MISSENSE_OUTPUT)
    logger.info("Report: %s", REPORT_OUTPUT)


if __name__ == "__main__":
    main()
