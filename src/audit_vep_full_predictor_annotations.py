from __future__ import annotations

import hashlib
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

ANNOTATION_FILE = (
    ROOT / "data/processed/vep_full_predictor_annotations.csv"
)
COHORT_FILE = (
    ROOT / "data/processed/brca1_brca2_primary_missense_cohort.csv"
)
PRIMARY_FILE = (
    ROOT / "data/processed/vep_primary_transcripts.csv"
)
FAILURE_FILE = (
    ROOT / "data/processed/vep_full_predictor_failures.csv"
)
PRODUCTION_REPORT_FILE = (
    ROOT / "data/metadata/vep_full_predictor_report.json"
)
OUTPUT_REPORT = (
    ROOT / "data/metadata/vep_full_predictor_integrity_audit.json"
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def normalize_id(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .str.strip()
        .str.replace(r"\.0$", "", regex=True)
    )


def normalize_text(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip()


def record_check(
    checks: list[dict],
    name: str,
    passed: bool,
    details: dict | str,
) -> None:
    checks.append(
        {
            "check": name,
            "passed": bool(passed),
            "details": details,
        }
    )

    status = "PASS" if passed else "FAIL"
    logger.info("%s | %s | %s", status, name, details)


def main() -> None:
    logger.info("Loading annotation and source datasets.")

    annotations = pd.read_csv(
        ANNOTATION_FILE,
        low_memory=False,
    )
    cohort = pd.read_csv(
        COHORT_FILE,
        low_memory=False,
    )
    primary = pd.read_csv(
        PRIMARY_FILE,
        low_memory=False,
    )
    failures = pd.read_csv(
        FAILURE_FILE,
        low_memory=False,
    )

    with PRODUCTION_REPORT_FILE.open(encoding="utf-8") as handle:
        production_report = json.load(handle)

    checks: list[dict] = []

    required_annotation_columns = [
        "#AlleleID",
        "GeneSymbol",
        "variant_key",
        "vep_identity_key",
        "Chromosome",
        "PositionVCF",
        "ReferenceAlleleVCF",
        "AlternateAlleleVCF",
        "transcript_id",
        "mane_select",
        "hgvsp",
        "cadd_raw",
        "cadd_phred",
        "revel",
        "cadd_raw_available",
        "cadd_phred_available",
        "revel_available",
        "mane_target_gene_validated",
        "assembly_name",
    ]

    required_cohort_columns = [
        "#AlleleID",
        "GeneSymbol",
        "variant_key",
        "transcript_id",
        "hgvsp",
    ]

    required_primary_columns = [
        "#AlleleID",
        "GeneSymbol",
        "variant_key",
        "transcript_id",
        "hgvsp",
    ]

    missing_annotation_columns = sorted(
        set(required_annotation_columns) - set(annotations.columns)
    )
    missing_cohort_columns = sorted(
        set(required_cohort_columns) - set(cohort.columns)
    )
    missing_primary_columns = sorted(
        set(required_primary_columns) - set(primary.columns)
    )

    record_check(
        checks,
        "required_columns_present",
        not (
            missing_annotation_columns
            or missing_cohort_columns
            or missing_primary_columns
        ),
        {
            "missing_annotation_columns": missing_annotation_columns,
            "missing_cohort_columns": missing_cohort_columns,
            "missing_primary_columns": missing_primary_columns,
        },
    )

    if any(
        [
            missing_annotation_columns,
            missing_cohort_columns,
            missing_primary_columns,
        ]
    ):
        raise RuntimeError(
            "Required columns are missing. Audit cannot continue."
        )

    expected_count = len(cohort)

    logger.info(
        "Cohort=%s | annotations=%s | primary transcripts=%s",
        len(cohort),
        len(annotations),
        len(primary),
    )

    # ---------------------------------------------------------
    # 1. Dataset sizes and uniqueness
    # ---------------------------------------------------------

    record_check(
        checks,
        "annotation_count_matches_cohort",
        len(annotations) == expected_count,
        {
            "cohort_count": expected_count,
            "annotation_count": len(annotations),
        },
    )

    record_check(
        checks,
        "cohort_variant_keys_unique",
        cohort["variant_key"].is_unique,
        {
            "duplicate_count": int(
                cohort["variant_key"].duplicated().sum()
            )
        },
    )

    record_check(
        checks,
        "annotation_variant_keys_unique",
        annotations["variant_key"].is_unique,
        {
            "duplicate_count": int(
                annotations["variant_key"].duplicated().sum()
            )
        },
    )

    record_check(
        checks,
        "annotation_allele_ids_unique",
        normalize_id(annotations["#AlleleID"]).is_unique,
        {
            "duplicate_count": int(
                normalize_id(annotations["#AlleleID"]).duplicated().sum()
            )
        },
    )

    # ---------------------------------------------------------
    # 2. Exact cohort-to-annotation reconciliation
    # ---------------------------------------------------------

    cohort_keys = set(cohort["variant_key"].dropna())
    annotation_keys = set(annotations["variant_key"].dropna())

    missing_keys = sorted(cohort_keys - annotation_keys)
    unexpected_keys = sorted(annotation_keys - cohort_keys)

    record_check(
        checks,
        "variant_key_membership_matches",
        not missing_keys and not unexpected_keys,
        {
            "missing_count": len(missing_keys),
            "unexpected_count": len(unexpected_keys),
            "missing_examples": missing_keys[:10],
            "unexpected_examples": unexpected_keys[:10],
        },
    )

    # Merge with explicit one-to-one validation.
    merged = cohort[
        [
            "#AlleleID",
            "GeneSymbol",
            "variant_key",
            "transcript_id",
            "hgvsp",
        ]
    ].merge(
        annotations,
        on="variant_key",
        how="outer",
        suffixes=("_cohort", "_annotation"),
        indicator=True,
        validate="one_to_one",
    )

    record_check(
        checks,
        "cohort_annotation_one_to_one_merge",
        bool((merged["_merge"] == "both").all()),
        merged["_merge"].value_counts().to_dict(),
    )

    matched = merged.loc[merged["_merge"] == "both"].copy()

    # ---------------------------------------------------------
    # 3. ClinVar identity and gene reconciliation
    # ---------------------------------------------------------

    cohort_ids = normalize_id(matched["#AlleleID_cohort"])
    annotation_ids = normalize_id(matched["#AlleleID_annotation"])

    allele_id_mismatches = int(
        (cohort_ids != annotation_ids).sum()
    )

    gene_mismatches = int(
        (
            normalize_text(matched["GeneSymbol_cohort"])
            != normalize_text(matched["GeneSymbol_annotation"])
        ).sum()
    )

    record_check(
        checks,
        "clinvar_allele_ids_match",
        allele_id_mismatches == 0,
        {"mismatch_count": allele_id_mismatches},
    )

    record_check(
        checks,
        "gene_assignments_match",
        gene_mismatches == 0,
        {"mismatch_count": gene_mismatches},
    )

    # ---------------------------------------------------------
    # 4. Genomic identity reconstruction
    # ---------------------------------------------------------

    reconstructed_keys = (
        "GRCh38:"
        + annotations["Chromosome"].astype("string").str.replace(
            r"^chr", "", regex=True, case=False
        )
        + ":"
        + pd.to_numeric(
            annotations["PositionVCF"], errors="coerce"
        ).astype("Int64").astype("string")
        + ":"
        + annotations["ReferenceAlleleVCF"].astype("string").str.upper()
        + ">"
        + annotations["AlternateAlleleVCF"].astype("string").str.upper()
    )

    identity_mismatches = int(
        (
            reconstructed_keys.reset_index(drop=True)
            != annotations["variant_key"].astype("string").reset_index(
                drop=True
            )
        ).sum()
    )

    record_check(
        checks,
        "genomic_variant_keys_reconstruct",
        identity_mismatches == 0,
        {"mismatch_count": identity_mismatches},
    )

    expected_vep_identity = (
        annotations["Chromosome"].astype("string").str.replace(
            r"^chr", "", regex=True, case=False
        )
        + ":"
        + pd.to_numeric(
            annotations["PositionVCF"], errors="coerce"
        ).astype("Int64").astype("string")
        + ":"
        + annotations["ReferenceAlleleVCF"].astype("string").str.upper()
        + ":"
        + annotations["AlternateAlleleVCF"].astype("string").str.upper()
    )

    vep_identity_mismatches = int(
        (
            expected_vep_identity.reset_index(drop=True)
            != annotations["vep_identity_key"].astype("string").reset_index(
                drop=True
            )
        ).sum()
    )

    record_check(
        checks,
        "vep_identity_keys_reconstruct",
        vep_identity_mismatches == 0,
        {"mismatch_count": vep_identity_mismatches},
    )

    # ---------------------------------------------------------
    # 5. Transcript consistency against cohort
    # ---------------------------------------------------------

    transcript_mismatches = int(
        (
            normalize_text(matched["transcript_id_cohort"])
            != normalize_text(matched["transcript_id_annotation"])
        ).sum()
    )

    hgvsp_mismatches = int(
        (
            normalize_text(matched["hgvsp_cohort"])
            != normalize_text(matched["hgvsp_annotation"])
        ).sum()
    )

    record_check(
        checks,
        "transcript_ids_match_cohort",
        transcript_mismatches == 0,
        {"mismatch_count": transcript_mismatches},
    )

    record_check(
        checks,
        "protein_changes_match_cohort",
        hgvsp_mismatches == 0,
        {"mismatch_count": hgvsp_mismatches},
    )

    # Independently compare against the selected primary transcript file.
    primary_missense = primary[
        primary["variant_key"].isin(cohort_keys)
    ][
        [
            "#AlleleID",
            "GeneSymbol",
            "variant_key",
            "transcript_id",
            "hgvsp",
        ]
    ].copy()

    record_check(
        checks,
        "primary_transcript_subset_count_matches",
        len(primary_missense) == expected_count,
        {
            "primary_subset_count": len(primary_missense),
            "cohort_count": expected_count,
        },
    )

    primary_merge = cohort[
        [
            "#AlleleID",
            "GeneSymbol",
            "variant_key",
            "transcript_id",
            "hgvsp",
        ]
    ].merge(
        primary_missense,
        on="variant_key",
        how="outer",
        suffixes=("_cohort", "_primary"),
        indicator=True,
        validate="one_to_one",
    )

    primary_matched = primary_merge.loc[
        primary_merge["_merge"] == "both"
    ].copy()

    primary_transcript_mismatches = int(
        (
            normalize_text(primary_matched["transcript_id_cohort"])
            != normalize_text(primary_matched["transcript_id_primary"])
        ).sum()
    )

    primary_protein_mismatches = int(
        (
            normalize_text(primary_matched["hgvsp_cohort"])
            != normalize_text(primary_matched["hgvsp_primary"])
        ).sum()
    )

    record_check(
        checks,
        "primary_transcript_ids_match_cohort",
        primary_transcript_mismatches == 0
        and len(primary_matched) == expected_count,
        {
            "mismatch_count": primary_transcript_mismatches,
            "matched_count": len(primary_matched),
        },
    )

    record_check(
        checks,
        "primary_protein_changes_match_cohort",
        primary_protein_mismatches == 0,
        {"mismatch_count": primary_protein_mismatches},
    )

    # ---------------------------------------------------------
    # 6. MANE and assembly validation
    # ---------------------------------------------------------

    mane_flags = (
        annotations["mane_target_gene_validated"]
        .astype("string")
        .str.lower()
        .isin(["true", "1"])
    )

    record_check(
        checks,
        "all_annotations_target_gene_mane_validated",
        bool(mane_flags.all()),
        {
            "validated_count": int(mane_flags.sum()),
            "total": len(annotations),
        },
    )

    record_check(
        checks,
        "all_annotations_use_grch38",
        bool(
            annotations["assembly_name"]
            .astype("string")
            .eq("GRCh38")
            .all()
        ),
        annotations["assembly_name"].value_counts(
            dropna=False
        ).to_dict(),
    )

    record_check(
        checks,
        "transcript_ids_present",
        bool(
            annotations["transcript_id"].notna().all()
            and annotations["transcript_id"].astype("string")
            .str.strip().ne("").all()
        ),
        {
            "missing_count": int(
                annotations["transcript_id"].isna().sum()
            )
        },
    )

    # ---------------------------------------------------------
    # 7. Predictor score validation
    # ---------------------------------------------------------

    predictor_columns = [
        "cadd_raw",
        "cadd_phred",
        "revel",
    ]

    predictor_summary = {}

    for column in predictor_columns:
        values = pd.to_numeric(
            annotations[column],
            errors="coerce",
        )

        finite_mask = np.isfinite(values.to_numpy(dtype=float))
        missing_count = int(values.isna().sum())
        nonfinite_count = int((~finite_mask & values.notna()).sum())

        predictor_summary[column] = {
            "count": int(values.notna().sum()),
            "missing_count": missing_count,
            "nonfinite_count": nonfinite_count,
            "minimum": (
                float(values.min()) if values.notna().any() else None
            ),
            "median": (
                float(values.median()) if values.notna().any() else None
            ),
            "mean": (
                float(values.mean()) if values.notna().any() else None
            ),
            "maximum": (
                float(values.max()) if values.notna().any() else None
            ),
        }

        record_check(
            checks,
            f"{column}_complete",
            missing_count == 0,
            predictor_summary[column],
        )

        record_check(
            checks,
            f"{column}_finite",
            nonfinite_count == 0,
            {"nonfinite_count": nonfinite_count},
        )

    revel_values = pd.to_numeric(
        annotations["revel"], errors="coerce"
    )

    revel_out_of_range = int(
        (
            revel_values.notna()
            & ~revel_values.between(0, 1)
        ).sum()
    )

    record_check(
        checks,
        "revel_within_zero_and_one",
        revel_out_of_range == 0,
        {"out_of_range_count": revel_out_of_range},
    )

    # CADD raw is intentionally not restricted to [0, 1].
    cadd_raw_values = pd.to_numeric(
        annotations["cadd_raw"], errors="coerce"
    )

    record_check(
        checks,
        "cadd_raw_range_reported_without_zero_one_constraint",
        bool(np.isfinite(cadd_raw_values.dropna()).all()),
        {
            "minimum": float(cadd_raw_values.min()),
            "maximum": float(cadd_raw_values.max()),
            "negative_values_are_permitted": True,
        },
    )

    # ---------------------------------------------------------
    # 8. Availability flags agree with score availability
    # ---------------------------------------------------------

    availability_summary = {}

    for column in predictor_columns:
        flag_column = f"{column}_available"

        flags = (
            annotations[flag_column]
            .astype("string")
            .str.lower()
            .isin(["true", "1"])
        )

        scores_present = annotations[column].notna()

        mismatch_count = int(
            (flags != scores_present).sum()
        )

        availability_summary[flag_column] = {
            "true_count": int(flags.sum()),
            "false_count": int((~flags).sum()),
            "flag_score_mismatches": mismatch_count,
        }

        record_check(
            checks,
            f"{flag_column}_matches_score_availability",
            mismatch_count == 0,
            availability_summary[flag_column],
        )

    # ---------------------------------------------------------
    # 9. Gene-specific score coverage
    # ---------------------------------------------------------

    gene_coverage = {}

    for gene, group in annotations.groupby("GeneSymbol"):
        gene_coverage[gene] = {
            "variant_count": int(len(group)),
            "cadd_raw_available": int(group["cadd_raw"].notna().sum()),
            "cadd_phred_available": int(
                group["cadd_phred"].notna().sum()
            ),
            "revel_available": int(group["revel"].notna().sum()),
        }

    # ---------------------------------------------------------
    # 10. Production report and failure file consistency
    # ---------------------------------------------------------

    report_consistency = {
        "input_variant_count": production_report.get(
            "input_variant_count"
        ),
        "successful_batches": production_report.get(
            "successful_batches"
        ),
        "failed_batches": production_report.get(
            "failed_batches"
        ),
        "selected_target_gene_mane_count": production_report.get(
            "selected_target_gene_mane_count"
        ),
        "failed_variant_count": production_report.get(
            "failed_variant_count"
        ),
    }

    production_report_pass = (
        production_report.get("input_variant_count") == expected_count
        and production_report.get("successful_batches") == 81
        and production_report.get("failed_batches") == 0
        and production_report.get(
            "selected_target_gene_mane_count"
        ) == expected_count
        and production_report.get("failed_variant_count") == 0
    )

    record_check(
        checks,
        "production_report_matches_expected_completion",
        production_report_pass,
        report_consistency,
    )

    record_check(
        checks,
        "failure_file_empty",
        len(failures) == 0,
        {"failure_rows": len(failures)},
    )

    # ---------------------------------------------------------
    # 11. Provenance and checksums
    # ---------------------------------------------------------

    input_hashes = {
        "annotations": sha256_file(ANNOTATION_FILE),
        "cohort": sha256_file(COHORT_FILE),
        "primary_transcripts": sha256_file(PRIMARY_FILE),
        "failures": sha256_file(FAILURE_FILE),
        "production_report": sha256_file(PRODUCTION_REPORT_FILE),
    }

    report = {
        "audit_name": "VEP full predictor annotation integrity audit",
        "created_at_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "assembly": "GRCh38",
        "expected_variant_count": expected_count,
        "annotation_variant_count": len(annotations),
        "primary_transcript_subset_count": len(primary_missense),
        "predictor_summary": predictor_summary,
        "availability_summary": availability_summary,
        "gene_coverage": gene_coverage,
        "production_report_summary": report_consistency,
        "checks": checks,
        "total_checks": len(checks),
        "passed_checks": sum(c["passed"] for c in checks),
        "failed_checks": sum(not c["passed"] for c in checks),
        "all_checks_passed": all(c["passed"] for c in checks),
        "sha256": input_hashes,
    }

    OUTPUT_REPORT.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_REPORT.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)

    logger.info("=" * 72)
    logger.info("VEP PREDICTOR INTEGRITY AUDIT SUMMARY")
    logger.info("=" * 72)
    logger.info("Total checks: %s", report["total_checks"])
    logger.info("Passed: %s", report["passed_checks"])
    logger.info("Failed: %s", report["failed_checks"])
    logger.info("Report: %s", OUTPUT_REPORT)

    if report["all_checks_passed"]:
        logger.info("ALL VEP PREDICTOR INTEGRITY CHECKS PASSED.")
    else:
        logger.error("VEP PREDICTOR INTEGRITY AUDIT FAILED.")

        for check in checks:
            if not check["passed"]:
                logger.error(
                    "FAILED CHECK: %s | %s",
                    check["check"],
                    check["details"],
                )

        raise SystemExit(1)


if __name__ == "__main__":
    main()
