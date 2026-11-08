
import json
import logging
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

VARIANTS_FILE = (
    ROOT / "data/processed/vep_full_variants.csv"
)
TRANSCRIPTS_FILE = (
    ROOT / "data/processed/vep_full_transcript_consequences.csv"
)

OUTPUT_FILE = (
    ROOT / "data/processed/vep_primary_transcripts.csv"
)
EXCEPTIONS_FILE = (
    ROOT / "data/metadata/vep_transcript_selection_exceptions.csv"
)
REPORT_FILE = (
    ROOT / "data/metadata/vep_transcript_selection_report.json"
)

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(message)s"

logging.basicConfig(
    level=logging.INFO,
    format=LOG_FORMAT,
)

logger = logging.getLogger(__name__)


def main():

    logger.info("Loading VEP datasets...")

    variants = pd.read_csv(VARIANTS_FILE)

    transcripts = pd.read_csv(
        TRANSCRIPTS_FILE,
        low_memory=False,
    )

    required_variant_columns = {
        "#AlleleID",
        "vep_identity_key",
        "GeneSymbol",
    }

    required_transcript_columns = {
        "#AlleleID",
        "transcript_id",
        "mane_select",
        "canonical",
        "consequence_terms",
    }

    missing_variant_columns = (
        required_variant_columns - set(variants.columns)
    )

    missing_transcript_columns = (
        required_transcript_columns - set(transcripts.columns)
    )

    if missing_variant_columns:
        raise ValueError(
            f"Missing variant columns: {missing_variant_columns}"
        )

    if missing_transcript_columns:
        raise ValueError(
            f"Missing transcript columns: {missing_transcript_columns}"
        )

    if variants["#AlleleID"].duplicated().any():
        raise ValueError("Duplicate AlleleIDs in VEP variants.")

    if transcripts["transcript_id"].isna().any():
        raise ValueError("Missing transcript IDs detected.")

    logger.info("Variants loaded: %s", len(variants))
    logger.info("Transcript records loaded: %s", len(transcripts))

    # Identify MANE Select and canonical transcript records.
    mane_mask = (
        transcripts["mane_select"].notna()
        & transcripts["mane_select"].astype(str).str.strip().ne("")
    )

    canonical_values = pd.to_numeric(
        transcripts["canonical"],
        errors="coerce",
    )

    canonical_mask = canonical_values.eq(1)

    decisions = []
    selected_indices = []
    exceptions = []

    logger.info("Selecting primary transcripts...")

    for allele_id, group in transcripts.groupby(
        "#AlleleID",
        sort=False,
    ):

        mane_rows = group.loc[mane_mask.loc[group.index]]
        canonical_rows = group.loc[canonical_mask.loc[group.index]]

        selected_index = None

        if len(mane_rows) == 1:

            selected_index = mane_rows.index[0]
            method = "MANE_Select"

        elif len(mane_rows) > 1:

            method = "manual_review"
            reason = "Multiple MANE Select transcript records"

            exceptions.append({
                "#AlleleID": allele_id,
                "reason": reason,
                "mane_count": len(mane_rows),
                "canonical_count": len(canonical_rows),
            })

        elif len(canonical_rows) == 1:

            selected_index = canonical_rows.index[0]
            method = "canonical_fallback"

        elif len(canonical_rows) > 1:

            method = "manual_review"
            reason = "Multiple canonical transcript records"

            exceptions.append({
                "#AlleleID": allele_id,
                "reason": reason,
                "mane_count": 0,
                "canonical_count": len(canonical_rows),
            })

        else:

            method = "manual_review"
            reason = "Neither MANE Select nor canonical transcript available"

            exceptions.append({
                "#AlleleID": allele_id,
                "reason": reason,
                "mane_count": 0,
                "canonical_count": 0,
            })

        if selected_index is not None:
            selected_indices.append(selected_index)

        decisions.append({
            "#AlleleID": allele_id,
            "transcript_selection_method": method,
        })

    # Keep only the transcript-specific fields when joining.
    base_columns = set(variants.columns)

    transcript_specific_columns = [
        column
        for column in transcripts.columns
        if column not in base_columns
        and column != "#AlleleID"
    ]

    selected_records = transcripts.loc[
        selected_indices,
        ["#AlleleID"] + transcript_specific_columns,
    ].copy()

    if selected_records["#AlleleID"].duplicated().any():
        raise ValueError(
            "More than one selected transcript per AlleleID."
        )

    # Retain every original VEP variant, including manual-review cases.
    primary = variants.merge(
        selected_records,
        on="#AlleleID",
        how="left",
        validate="one_to_one",
    )

    decision_df = pd.DataFrame(decisions)

    primary = primary.merge(
        decision_df,
        on="#AlleleID",
        how="left",
        validate="one_to_one",
    )

    # Verify output integrity.
    if len(primary) != len(variants):
        raise ValueError(
            "Primary transcript table does not match variant count."
        )

    if primary["#AlleleID"].duplicated().any():
        raise ValueError(
            "Duplicate AlleleIDs in primary transcript output."
        )

    # Save exceptions for manual inspection.
    exception_df = pd.DataFrame(
        exceptions,
        columns=[
            "#AlleleID",
            "reason",
            "mane_count",
            "canonical_count",
        ],
    )

    exception_df.to_csv(EXCEPTIONS_FILE, index=False)

    primary.to_csv(OUTPUT_FILE, index=False)

    # Generate selection summary.
    selection_counts = (
        primary["transcript_selection_method"]
        .value_counts()
        .to_dict()
    )

    gene_counts = {}

    for gene, group in primary.groupby("GeneSymbol"):
        gene_counts[str(gene)] = {
            "total_variants": int(len(group)),
            "mane_select": int(
                (group["transcript_selection_method"] == "MANE_Select").sum()
            ),
            "canonical_fallback": int(
                (
                    group["transcript_selection_method"]
                    == "canonical_fallback"
                ).sum()
            ),
            "manual_review": int(
                (
                    group["transcript_selection_method"]
                    == "manual_review"
                ).sum()
            ),
        }

    report = {
        "input_variant_count": int(len(variants)),
        "input_transcript_record_count": int(len(transcripts)),
        "primary_transcript_output_count": int(len(primary)),
        "selection_method_distribution": {
            key: int(value)
            for key, value in selection_counts.items()
        },
        "gene_distribution": gene_counts,
        "manual_review_count": int(len(exception_df)),
        "selection_policy": [
            "MANE Select preferred when exactly one is available",
            "Canonical fallback only when MANE Select is absent",
            "Ambiguous or unavailable transcript selection flagged for manual review",
            "All original variants retained",
            "All original transcript records preserved separately",
        ],
        "output_file": str(OUTPUT_FILE),
        "exceptions_file": str(EXCEPTIONS_FILE),
    }

    REPORT_FILE.write_text(
        json.dumps(report, indent=4),
        encoding="utf-8",
    )

    logger.info("==============================================")
    logger.info("PRIMARY TRANSCRIPT SELECTION COMPLETED")
    logger.info("==============================================")
    logger.info("Variants retained: %s", len(primary))
    logger.info("Selection distribution: %s", selection_counts)
    logger.info("Manual review cases: %s", len(exception_df))
    logger.info("Output: %s", OUTPUT_FILE)
    logger.info("Report: %s", REPORT_FILE)
    logger.info("Exceptions: %s", EXCEPTIONS_FILE)


if __name__ == "__main__":
    main()
