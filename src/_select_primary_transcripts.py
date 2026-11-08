
import json
import logging
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

VARIANTS_FILE = ROOT / "data/processed/vep_full_variants.csv"
TRANSCRIPTS_FILE = ROOT / "data/processed/vep_full_transcript_consequences.csv"

OUTPUT_FILE = ROOT / "data/processed/vep_primary_transcripts.csv"
EXCEPTIONS_FILE = ROOT / "data/metadata/vep_transcript_selection_exceptions.csv"
REPORT_FILE = ROOT / "data/metadata/vep_transcript_selection_report.json"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


def normalize_gene(series):
    return series.astype("string").str.strip().str.upper()


def main():

    logger.info("Loading VEP datasets...")

    variants = pd.read_csv(VARIANTS_FILE, low_memory=False)
    transcripts = pd.read_csv(TRANSCRIPTS_FILE, low_memory=False)

    required_variant_columns = {
        "#AlleleID",
        "vep_identity_key",
        "GeneSymbol",
    }

    required_transcript_columns = {
        "#AlleleID",
        "transcript_id",
        "gene_symbol",
        "mane_select",
        "canonical",
        "consequence_terms",
    }

    missing = required_variant_columns - set(variants.columns)
    if missing:
        raise ValueError(f"Missing variant columns: {missing}")

    missing = required_transcript_columns - set(transcripts.columns)
    if missing:
        raise ValueError(f"Missing transcript columns: {missing}")

    if variants["#AlleleID"].duplicated().any():
        raise ValueError("Duplicate AlleleIDs in variant dataset.")

    if transcripts["transcript_id"].isna().any():
        raise ValueError("Missing transcript IDs detected.")

    # Map every variant to its ClinVar target gene.
    gene_map = variants.set_index("#AlleleID")["GeneSymbol"]

    transcripts["target_gene"] = transcripts["#AlleleID"].map(gene_map)

    if transcripts["target_gene"].isna().any():
        raise ValueError("Transcript records contain unknown ClinVar AlleleIDs.")

    # Match transcript gene to the ClinVar target gene.
    transcripts["target_gene_normalized"] = normalize_gene(
        transcripts["target_gene"]
    )

    transcripts["transcript_gene_normalized"] = normalize_gene(
        transcripts["gene_symbol"]
    )

    transcripts["is_target_gene"] = (
        transcripts["target_gene_normalized"].notna()
        & transcripts["transcript_gene_normalized"].notna()
        & (
            transcripts["target_gene_normalized"]
            == transcripts["transcript_gene_normalized"]
        )
    )

    # Identify transcript annotations.
    mane_mask = (
        transcripts["mane_select"].notna()
        & transcripts["mane_select"].astype(str).str.strip().ne("")
    )

    canonical_mask = pd.to_numeric(
        transcripts["canonical"],
        errors="coerce",
    ).eq(1)

    target_gene_mask = transcripts["is_target_gene"]

    decisions = []
    selected_indices = []
    exceptions = []

    logger.info("Selecting transcripts belonging to ClinVar target genes...")

    for allele_id, group in transcripts.groupby("#AlleleID", sort=False):

        target_group = group.loc[group["is_target_gene"]]

        target_mane = target_group.loc[
            mane_mask.loc[target_group.index]
        ]

        target_canonical = target_group.loc[
            canonical_mask.loc[target_group.index]
        ]

        selected_index = None
        reason = ""

        if len(target_mane) == 1:

            selected_index = target_mane.index[0]
            method = "MANE_Select"

        elif len(target_mane) > 1:

            method = "manual_review"
            reason = "Multiple MANE Select transcripts in target gene"

        elif len(target_canonical) == 1:

            selected_index = target_canonical.index[0]
            method = "canonical_fallback"

        elif len(target_canonical) > 1:

            method = "manual_review"
            reason = "Multiple canonical transcripts in target gene"

        else:

            method = "manual_review"
            reason = "No MANE Select or canonical transcript in target gene"

        if selected_index is not None:
            selected_indices.append(selected_index)

        decisions.append({
            "#AlleleID": allele_id,
            "transcript_selection_method": method,
            "target_gene_transcript_count": int(len(target_group)),
            "target_gene_mane_count": int(len(target_mane)),
            "target_gene_canonical_count": int(len(target_canonical)),
        })

        if reason:
            exceptions.append({
                "#AlleleID": allele_id,
                "target_gene": group["target_gene"].iloc[0],
                "reason": reason,
                "target_gene_transcript_count": len(target_group),
                "target_gene_mane_count": len(target_mane),
                "target_gene_canonical_count": len(target_canonical),
            })

    # Extract the selected transcript rows.
    base_columns = set(variants.columns)

    transcript_specific_columns = [
        column
        for column in transcripts.columns
        if column not in base_columns
        and column not in {
            "target_gene",
            "target_gene_normalized",
            "transcript_gene_normalized",
            "is_target_gene",
        }
        and column != "#AlleleID"
    ]

    selected_records = transcripts.loc[
        selected_indices,
        ["#AlleleID"] + transcript_specific_columns,
    ].copy()

    if selected_records["#AlleleID"].duplicated().any():
        raise ValueError("More than one selected transcript per variant.")

    # Retain every original variant, including manual-review cases.
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

    # Integrity checks.
    if len(primary) != len(variants):
        raise ValueError("Variant count changed during transcript selection.")

    if primary["#AlleleID"].duplicated().any():
        raise ValueError("Duplicate AlleleIDs in primary output.")

    selected = primary[
        primary["transcript_selection_method"].isin(
            ["MANE_Select", "canonical_fallback"]
        )
    ]

    if not (
        normalize_gene(selected["gene_symbol"])
        == normalize_gene(selected["GeneSymbol"])
    ).all():
        raise ValueError(
            "A selected transcript does not belong to its ClinVar target gene."
        )

    # Save output and exceptions.
    exception_columns = [
        "#AlleleID",
        "target_gene",
        "reason",
        "target_gene_transcript_count",
        "target_gene_mane_count",
        "target_gene_canonical_count",
    ]

    exception_df = pd.DataFrame(
        exceptions,
        columns=exception_columns,
    )

    primary.to_csv(OUTPUT_FILE, index=False)
    exception_df.to_csv(EXCEPTIONS_FILE, index=False)

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
            str(key): int(value)
            for key, value in selection_counts.items()
        },
        "gene_distribution": gene_counts,
        "manual_review_count": int(len(exception_df)),
        "mane_records_from_other_genes_excluded": int(
            (mane_mask & ~target_gene_mask).sum()
        ),
        "selection_policy": [
            "Only transcripts belonging to the ClinVar target gene are eligible",
            "MANE Select preferred when exactly one target-gene record is available",
            "Canonical fallback only when target-gene MANE Select is absent",
            "Ambiguous or unavailable target-gene transcripts flagged for manual review",
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
    logger.info("TARGET-GENE TRANSCRIPT SELECTION COMPLETED")
    logger.info("==============================================")
    logger.info("Variants retained: %s", len(primary))
    logger.info("Selection distribution: %s", selection_counts)
    logger.info(
        "MANE records from other genes excluded: %s",
        report["mane_records_from_other_genes_excluded"],
    )
    logger.info("Manual review cases: %s", len(exception_df))
    logger.info("Output: %s", OUTPUT_FILE)
    logger.info("Report: %s", REPORT_FILE)
    logger.info("Exceptions: %s", EXCEPTIONS_FILE)


if __name__ == "__main__":
    main()
