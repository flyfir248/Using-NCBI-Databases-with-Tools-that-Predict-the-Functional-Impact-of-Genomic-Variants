
import json
import logging
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

CANDIDATES_FILE = ROOT / "data/processed/brca1_brca2_grch38_candidate_snvs.csv"
VEP_VARIANTS_FILE = ROOT / "data/processed/vep_full_variants.csv"
VEP_TRANSCRIPTS_FILE = ROOT / "data/processed/vep_full_transcript_consequences.csv"
VEP_FAILURES_FILE = ROOT / "data/processed/vep_full_failures.csv"
COVERAGE_FILE = ROOT / "data/metadata/vep_full_transcript_coverage.csv"

REPORT_FILE = ROOT / "data/metadata/vep_full_integrity_audit.json"
MANE_MISSING_FILE = ROOT / "data/metadata/vep_variants_without_mane_select.csv"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


def normalize_chromosome(value):
    value = str(value).strip()
    return value.removeprefix("chr")


def candidate_identity(row):
    return (
        f"{normalize_chromosome(row['Chromosome'])}:"
        f"{int(row['PositionVCF'])}:"
        f"{str(row['ReferenceAlleleVCF']).upper()}:"
        f"{str(row['AlternateAlleleVCF']).upper()}"
    )


def vep_identity(row):
    alleles = str(row["allele_string"]).split("/")

    if len(alleles) != 2:
        return None

    return (
        f"{normalize_chromosome(row['seq_region_name'])}:"
        f"{int(row['start'])}:"
        f"{alleles[0].upper()}:"
        f"{alleles[1].upper()}"
    )


def main():
    logger.info("Loading all datasets...")

    candidates = pd.read_csv(
        CANDIDATES_FILE,
        dtype={"#AlleleID": "string"},
    )
    variants = pd.read_csv(
        VEP_VARIANTS_FILE,
        dtype={"#AlleleID": "string"},
    )
    transcripts = pd.read_csv(
        VEP_TRANSCRIPTS_FILE,
        dtype={"#AlleleID": "string", "transcript_id": "string"},
    )
    failures = pd.read_csv(
        VEP_FAILURES_FILE,
        dtype={"#AlleleID": "string"},
    )
    coverage = pd.read_csv(
        COVERAGE_FILE,
        dtype={"#AlleleID": "string"},
    )

    report = {
        "input_candidate_count": len(candidates),
        "vep_variant_count": len(variants),
        "transcript_row_count": len(transcripts),
        "failure_row_count": len(failures),
        "checks": {},
        "errors": [],
    }

    def record_check(name, passed, details=None):
        report["checks"][name] = {
            "passed": bool(passed),
            "details": details or {},
        }
        logger.info(
            "%s: %s",
            "PASS" if passed else "FAIL",
            name,
        )
        if not passed:
            report["errors"].append(name)

    # ---------------------------------------------------------
    # 1. Candidate and VEP record counts
    # ---------------------------------------------------------

    record_check(
        "candidate_count_matches_vep_count",
        len(candidates) == len(variants),
        {
            "candidates": len(candidates),
            "vep_variants": len(variants),
        },
    )

    record_check(
        "no_failed_variants",
        len(failures) == 0,
        {"failure_rows": len(failures)},
    )

    # ---------------------------------------------------------
    # 2. ClinVar AlleleID uniqueness
    # ---------------------------------------------------------

    candidate_ids = set(candidates["#AlleleID"].dropna())
    vep_ids = set(variants["#AlleleID"].dropna())

    record_check(
        "candidate_allele_ids_unique",
        candidates["#AlleleID"].is_unique,
        {"unique_ids": len(candidate_ids)},
    )

    record_check(
        "vep_allele_ids_unique",
        variants["#AlleleID"].is_unique,
        {"unique_ids": len(vep_ids)},
    )

    record_check(
        "all_candidate_allele_ids_present_in_vep",
        candidate_ids == vep_ids,
        {
            "missing_from_vep": len(candidate_ids - vep_ids),
            "unexpected_in_vep": len(vep_ids - candidate_ids),
        },
    )

    # ---------------------------------------------------------
    # 3. Independent genomic identity reconstruction
    # ---------------------------------------------------------

    candidates["audit_identity_key"] = candidates.apply(
        candidate_identity,
        axis=1,
    )

    variants["audit_identity_key"] = variants.apply(
        vep_identity,
        axis=1,
    )

    record_check(
        "candidate_genomic_keys_unique",
        candidates["audit_identity_key"].is_unique,
        {
            "unique_keys": candidates["audit_identity_key"].nunique()
        },
    )

    record_check(
        "vep_genomic_keys_unique",
        variants["audit_identity_key"].is_unique,
        {
            "unique_keys": variants["audit_identity_key"].nunique()
        },
    )

    candidate_keys = set(candidates["audit_identity_key"])
    vep_keys = set(variants["audit_identity_key"].dropna())

    record_check(
        "all_candidate_genomic_keys_present_in_vep",
        candidate_keys == vep_keys,
        {
            "missing_from_vep": len(candidate_keys - vep_keys),
            "unexpected_in_vep": len(vep_keys - candidate_keys),
        },
    )

    record_check(
        "stored_vep_identity_key_matches_reconstructed_key",
        (
            variants["vep_identity_key"].astype(str)
            == variants["audit_identity_key"].astype(str)
        ).all(),
        {
            "mismatches": int(
                (
                    variants["vep_identity_key"].astype(str)
                    != variants["audit_identity_key"].astype(str)
                ).sum()
            )
        },
    )

    # ---------------------------------------------------------
    # 4. ClinVar-to-VEP identity mapping by AlleleID
    # ---------------------------------------------------------

    identity_comparison = candidates[
        ["#AlleleID", "audit_identity_key"]
    ].merge(
        variants[["#AlleleID", "audit_identity_key"]],
        on="#AlleleID",
        how="outer",
        suffixes=("_clinvar", "_vep"),
        indicator=True,
    )

    mapping_mismatches = identity_comparison[
        (identity_comparison["_merge"] != "both")
        | (
            identity_comparison["audit_identity_key_clinvar"]
            != identity_comparison["audit_identity_key_vep"]
        )
    ]

    record_check(
        "allele_id_to_genomic_identity_mapping",
        len(mapping_mismatches) == 0,
        {"mapping_mismatches": len(mapping_mismatches)},
    )

    # ---------------------------------------------------------
    # 5. Transcript linkage and transcript counts
    # ---------------------------------------------------------

    valid_variant_ids = set(variants["#AlleleID"].dropna())
    transcript_ids = set(transcripts["#AlleleID"].dropna())

    orphan_transcripts = transcript_ids - valid_variant_ids

    record_check(
        "no_orphan_transcript_records",
        len(orphan_transcripts) == 0,
        {"orphan_allele_ids": len(orphan_transcripts)},
    )

    transcript_counts = (
        transcripts.groupby("#AlleleID")
        .size()
        .rename("observed_transcript_count")
    )

    count_comparison = variants[
        ["#AlleleID", "transcript_consequence_count"]
    ].merge(
        transcript_counts,
        on="#AlleleID",
        how="left",
    )

    count_comparison["observed_transcript_count"] = (
        count_comparison["observed_transcript_count"].fillna(0)
    )

    count_mismatches = count_comparison[
        count_comparison["transcript_consequence_count"]
        != count_comparison["observed_transcript_count"]
    ]

    record_check(
        "transcript_counts_match_variant_records",
        len(count_mismatches) == 0,
        {"count_mismatches": len(count_mismatches)},
    )

    record_check(
        "no_missing_transcript_ids",
        transcripts["transcript_id"].notna().all(),
        {
            "missing_transcript_ids": int(
                transcripts["transcript_id"].isna().sum()
            )
        },
    )

    # ---------------------------------------------------------
    # 6. MANE Select and canonical transcript coverage
    # ---------------------------------------------------------

    mane_missing = coverage[
        ~coverage["has_mane_select"].astype(bool)
    ].copy()

    canonical_missing = coverage[
        ~coverage["has_canonical"].astype(bool)
    ].copy()

    record_check(
        "coverage_table_has_one_row_per_variant",
        (
            len(coverage) == len(variants)
            and coverage["#AlleleID"].is_unique
        ),
        {
            "coverage_rows": len(coverage),
            "unique_coverage_allele_ids": coverage["#AlleleID"].nunique(),
        },
    )

    record_check(
        "every_variant_has_mane_or_canonical",
        (
            coverage["has_mane_select"].astype(bool)
            | coverage["has_canonical"].astype(bool)
        ).all(),
        {
            "without_mane": len(mane_missing),
            "without_canonical": len(canonical_missing),
        },
    )

    mane_missing.to_csv(MANE_MISSING_FILE, index=False)

    report["transcript_coverage"] = {
        "variants_with_mane_select": int(
            coverage["has_mane_select"].astype(bool).sum()
        ),
        "variants_without_mane_select": len(mane_missing),
        "variants_with_canonical": int(
            coverage["has_canonical"].astype(bool).sum()
        ),
        "variants_without_canonical": len(canonical_missing),
        "mane_missing_variants": mane_missing[
            ["#AlleleID", "GeneSymbol", "transcript_count",
             "canonical_count"]
        ].to_dict(orient="records"),
    }

    # ---------------------------------------------------------
    # 7. Dataset summaries
    # ---------------------------------------------------------

    report["gene_distribution"] = (
        variants["GeneSymbol"].value_counts().to_dict()
    )

    report["most_severe_consequence_distribution"] = (
        variants["most_severe_consequence"]
        .value_counts(dropna=False)
        .to_dict()
    )

    report["transcript_biotype_distribution"] = (
        transcripts["biotype"]
        .value_counts(dropna=False)
        .to_dict()
    )

    report["transcript_consequence_distribution"] = (
        transcripts["consequence_terms"]
        .value_counts(dropna=False)
        .to_dict()
    )

    report["unique_transcript_ids"] = int(
        transcripts["transcript_id"].nunique()
    )

    report["unique_protein_ids"] = int(
        transcripts["protein_id"].nunique()
    )

    report["audit_passed"] = len(report["errors"]) == 0

    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)

    with open(REPORT_FILE, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, default=str)

    logger.info("=" * 65)
    logger.info("VEP INTEGRITY AUDIT %s",
                "PASSED" if report["audit_passed"] else "FAILED")
    logger.info("Report: %s", REPORT_FILE)
    logger.info("MANE-missing variants: %s", MANE_MISSING_FILE)
    logger.info("=" * 65)

    if report["errors"]:
        logger.error("Failed checks: %s", report["errors"])
        raise SystemExit(1)


if __name__ == "__main__":
    main()
