from pathlib import Path
import hashlib
import json
import logging
import re
from datetime import datetime, timezone

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

COHORT_FILE = (
    ROOT / "data/processed/"
    "brca1_brca2_primary_missense_cohort.csv"
)

ALPHA_FILE = (
    ROOT / "data/processed/"
    "alphamissense_cohort_source_matches.csv"
)

OUTPUT_FILE = (
    ROOT / "data/processed/"
    "brca1_brca2_missense_predictor_integrated.csv"
)

REPORT_FILE = (
    ROOT / "data/metadata/"
    "alphamissense_integration_report.json"
)

EXPECTED_VARIANTS = 16193

EXPECTED_CLASSES = {
    "likely_benign",
    "ambiguous",
    "likely_pathogenic",
}

AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D",
    "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)


def sha256_file(path):
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def parse_hgvsp(value):
    if pd.isna(value):
        return None

    match = re.search(
        r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})",
        str(value),
    )

    if not match:
        return None

    ref3, position, alt3 = match.groups()

    if ref3 not in AA3_TO_AA1 or alt3 not in AA3_TO_AA1:
        return None

    return f"{AA3_TO_AA1[ref3]}{position}{AA3_TO_AA1[alt3]}"


def parse_alpha_variant(value):
    if pd.isna(value):
        return None

    match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", str(value))

    if not match:
        return None

    return "".join(match.groups())


def main():
    logging.info("Loading MANE missense cohort.")
    cohort = pd.read_csv(COHORT_FILE, low_memory=False)

    logging.info("Loading validated AlphaMissense source matches.")
    alpha = pd.read_csv(ALPHA_FILE, low_memory=False)

    required_cohort = {
        "variant_key",
        "GeneSymbol",
        "transcript_id",
        "hgvsp",
    }

    required_alpha = {
        "source_variant_key",
        "transcript_id",
        "protein_variant",
        "uniprot_id",
        "am_pathogenicity",
        "am_class",
    }

    if not required_cohort.issubset(cohort.columns):
        raise ValueError(
            f"Missing cohort columns: "
            f"{required_cohort - set(cohort.columns)}"
        )

    if not required_alpha.issubset(alpha.columns):
        raise ValueError(
            f"Missing AlphaMissense columns: "
            f"{required_alpha - set(alpha.columns)}"
        )

    if len(cohort) != EXPECTED_VARIANTS:
        raise ValueError(
            f"Unexpected cohort size: {len(cohort)}"
        )

    if cohort["variant_key"].isna().any():
        raise ValueError("Missing cohort genomic keys.")

    if alpha["source_variant_key"].isna().any():
        raise ValueError("Missing AlphaMissense genomic keys.")

    if cohort["variant_key"].duplicated().any():
        raise ValueError("Duplicate cohort genomic keys.")

    if alpha["source_variant_key"].duplicated().any():
        raise ValueError("Duplicate AlphaMissense genomic keys.")

    alpha = alpha.rename(
        columns={
            "source_variant_key": "alphamissense_source_variant_key",
            "transcript_id": "alphamissense_transcript_id",
            "protein_variant": "alphamissense_protein_variant",
            "uniprot_id": "alphamissense_uniprot_id",
            "am_pathogenicity": "alphamissense_pathogenicity",
            "am_class": "alphamissense_class",
        }
    )

    alpha_columns = [
        "alphamissense_source_variant_key",
        "alphamissense_transcript_id",
        "alphamissense_protein_variant",
        "alphamissense_uniprot_id",
        "alphamissense_pathogenicity",
        "alphamissense_class",
    ]

    logging.info("Performing one-to-one genomic integration.")

    integrated = cohort.merge(
        alpha[alpha_columns],
        left_on="variant_key",
        right_on="alphamissense_source_variant_key",
        how="outer",
        validate="one_to_one",
        indicator=True,
    )

    merge_counts = {
        str(key): int(value)
        for key, value in integrated["_merge"].value_counts().items()
    }

    if not (integrated["_merge"] == "both").all():
        raise RuntimeError(
            f"Genomic integration mismatch: {merge_counts}"
        )

    integrated["alphamissense_genomic_match"] = True

    # Transcript accession comparisons.
    cohort_base = (
        integrated["transcript_id"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    alpha_base = (
        integrated["alphamissense_transcript_id"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    integrated["alphamissense_transcript_exact_match"] = (
        integrated["transcript_id"]
        == integrated["alphamissense_transcript_id"]
    )

    integrated["alphamissense_transcript_base_match"] = (
        cohort_base == alpha_base
    )

    # Protein substitution comparisons.
    integrated["mane_protein_change_normalized"] = (
        integrated["hgvsp"].apply(parse_hgvsp)
    )

    integrated["alphamissense_protein_change_normalized"] = (
        integrated["alphamissense_protein_variant"]
        .apply(parse_alpha_variant)
    )

    integrated["alphamissense_protein_change_comparable"] = (
        integrated["mane_protein_change_normalized"].notna()
        & integrated["alphamissense_protein_change_normalized"].notna()
    )

    integrated["alphamissense_protein_change_match"] = (
        integrated["mane_protein_change_normalized"]
        == integrated["alphamissense_protein_change_normalized"]
    )

    # Score and category validation.
    integrated["alphamissense_pathogenicity"] = pd.to_numeric(
        integrated["alphamissense_pathogenicity"],
        errors="coerce",
    )

    scores = integrated["alphamissense_pathogenicity"]
    classes = set(integrated["alphamissense_class"].dropna().unique())

    checks = {
        "expected_cohort_size": len(integrated) == EXPECTED_VARIANTS,
        "one_to_one_genomic_merge": (
            merge_counts.get("both", 0) == EXPECTED_VARIANTS
            and merge_counts.get("left_only", 0) == 0
            and merge_counts.get("right_only", 0) == 0
        ),
        "unique_genomic_keys": (
            not integrated["variant_key"].duplicated().any()
        ),
        "no_missing_scores": not scores.isna().any(),
        "finite_scores": bool(np.isfinite(scores).all()),
        "scores_within_zero_and_one": bool(
            scores.between(0, 1).all()
        ),
        "no_missing_classes": (
            not integrated["alphamissense_class"].isna().any()
        ),
        "expected_prediction_classes": classes == EXPECTED_CLASSES,
        "all_protein_changes_comparable": bool(
            integrated["alphamissense_protein_change_comparable"].all()
        ),
        "all_protein_changes_match": bool(
            integrated["alphamissense_protein_change_match"].all()
        ),
        "all_genomic_matches_validated": bool(
            integrated["alphamissense_genomic_match"].all()
        ),
    }

    if not all(checks.values()):
        raise RuntimeError(
            "AlphaMissense integration validation failed: "
            + json.dumps(checks, indent=2)
        )

    # Remove temporary merge indicator and helper columns.
    integrated = integrated.drop(
        columns=[
            "_merge",
            "mane_protein_change_normalized",
            "alphamissense_protein_change_normalized",
        ]
    )

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)

    integrated.to_csv(OUTPUT_FILE, index=False)

    # Verify the saved file, not just the in-memory dataframe.
    saved = pd.read_csv(OUTPUT_FILE, low_memory=False)

    if len(saved) != EXPECTED_VARIANTS:
        raise RuntimeError("Saved integration row count is incorrect.")

    if saved["variant_key"].duplicated().any():
        raise RuntimeError("Saved integration contains duplicate keys.")

    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_cohort": str(COHORT_FILE.relative_to(ROOT)),
        "input_alphamissense_matches": str(ALPHA_FILE.relative_to(ROOT)),
        "output_file": str(OUTPUT_FILE.relative_to(ROOT)),
        "input_cohort_sha256": sha256_file(COHORT_FILE),
        "input_alphamissense_matches_sha256": sha256_file(ALPHA_FILE),
        "output_sha256": sha256_file(OUTPUT_FILE),
        "total_variants": len(integrated),
        "genomic_merge_counts": merge_counts,
        "transcript_concordance": {
            "exact_accession_matches": int(
                integrated["alphamissense_transcript_exact_match"].sum()
            ),
            "version_independent_matches": int(
                integrated["alphamissense_transcript_base_match"].sum()
            ),
        },
        "protein_concordance": {
            "comparable": int(
                integrated["alphamissense_protein_change_comparable"].sum()
            ),
            "matching": int(
                integrated["alphamissense_protein_change_match"].sum()
            ),
            "mismatching": int(
                (
                    integrated["alphamissense_protein_change_comparable"]
                    & ~integrated["alphamissense_protein_change_match"]
                ).sum()
            ),
        },
        "score_summary": {
            "minimum": float(scores.min()),
            "median": float(scores.median()),
            "mean": float(scores.mean()),
            "maximum": float(scores.max()),
        },
        "class_distribution": {
            str(key): int(value)
            for key, value in integrated[
                "alphamissense_class"
            ].value_counts().items()
        },
        "validation_checks": checks,
        "all_checks_passed": True,
    }

    with open(REPORT_FILE, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=4)

    logging.info("ALPHAMISSENSE INTEGRATION COMPLETED.")
    logging.info("Integrated variants: %s", len(integrated))
    logging.info(
        "Protein substitutions matching: %s",
        int(integrated["alphamissense_protein_change_match"].sum()),
    )
    logging.info("Output: %s", OUTPUT_FILE)
    logging.info("Report: %s", REPORT_FILE)


if __name__ == "__main__":
    main()

