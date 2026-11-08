from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

BASE_FILE = ROOT / "data/processed/brca1_brca2_missense_predictor_integrated.csv"
VEP_FILE = ROOT / "data/processed/vep_full_predictor_annotations.csv"
BENCHMARK_FILE = ROOT / "data/processed/brca1_brca2_missense_benchmark_annotated.csv"

OUTPUT_FILE = ROOT / "data/processed/brca1_brca2_all_predictors_integrated.csv"
REPORT_FILE = ROOT / "data/metadata/all_predictors_integration_report.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    print("Loading existing integrated dataset...")
    base = pd.read_csv(BASE_FILE, low_memory=False)

    print("Loading CADD/REVEL annotations...")
    vep = pd.read_csv(VEP_FILE, low_memory=False)

    print("Loading ClinVar benchmark annotations...")
    benchmark = pd.read_csv(BENCHMARK_FILE, low_memory=False)

    for name, frame in [
        ("base", base),
        ("VEP", vep),
        ("benchmark", benchmark),
    ]:
        if "variant_key" not in frame.columns:
            raise ValueError(f"{name} dataset has no variant_key.")
        if frame["variant_key"].isna().any():
            raise ValueError(f"{name} dataset contains missing variant keys.")
        if frame["variant_key"].duplicated().any():
            raise ValueError(f"{name} dataset has duplicate variant keys.")

    expected = set(base["variant_key"])
    if set(vep["variant_key"]) != expected:
        raise ValueError("VEP variant membership differs from base dataset.")
    if set(benchmark["variant_key"]) != expected:
        raise ValueError("Benchmark variant membership differs from base dataset.")

    # Bring in only new VEP predictor values and their availability flags.
    vep_columns = [
        "variant_key",
        "cadd_raw",
        "cadd_phred",
        "revel",
        "cadd_raw_available",
        "cadd_phred_available",
        "revel_available",
    ]

    existing_collisions = sorted(
        (set(vep_columns) - {"variant_key"}) & set(base.columns)
    )
    if existing_collisions:
        raise ValueError(
            f"Base dataset already contains these fields: {existing_collisions}"
        )

    combined = base.merge(
        vep[vep_columns],
        on="variant_key",
        how="left",
        validate="one_to_one",
        indicator=True,
    )

    if not combined["_merge"].eq("both").all():
        raise ValueError("Not all base variants received VEP predictor scores.")

    combined.drop(columns="_merge", inplace=True)

    # Add ClinVar benchmark fields, without duplicating fields already present.
    benchmark_fields = [
        "variant_key",
        "benchmark_group",
        "clinvar_binary_label",
        "benchmark_eligible",
        "benchmark_exclusion_reason",
        "review_confidence_group",
        "high_confidence_eligible",
    ]
    benchmark_fields = [
        c for c in benchmark_fields if c in benchmark.columns
    ]

    new_benchmark_fields = [
        c for c in benchmark_fields
        if c == "variant_key" or c not in combined.columns
    ]

    combined = combined.merge(
        benchmark[new_benchmark_fields],
        on="variant_key",
        how="left",
        validate="one_to_one",
        indicator="_benchmark_merge",
    )

    if not combined["_benchmark_merge"].eq("both").all():
        raise ValueError("Benchmark annotation merge was incomplete.")

    combined.drop(columns="_benchmark_merge", inplace=True)

    # Validate five expected predictor score fields.
    predictors = [
        "sift_score",
        "polyphen_score",
        "alphamissense_pathogenicity",
        "cadd_phred",
        "revel",
    ]

    missing_predictors = [
        c for c in predictors if c not in combined.columns
    ]
    if missing_predictors:
        raise ValueError(
            f"Expected predictor columns are missing: {missing_predictors}"
        )

    score_summary = {}
    for column in predictors:
        values = pd.to_numeric(combined[column], errors="coerce")
        finite = np.isfinite(values.dropna().to_numpy(dtype=float))

        score_summary[column] = {
            "available": int(values.notna().sum()),
            "missing": int(values.isna().sum()),
            "finite_nonmissing": int(finite.sum()),
            "minimum": float(values.min()) if values.notna().any() else None,
            "median": float(values.median()) if values.notna().any() else None,
            "maximum": float(values.max()) if values.notna().any() else None,
        }

        if not finite.all():
            raise ValueError(f"{column} contains non-finite values.")

    # Predictor-specific domain checks.
    for column in [
        "sift_score",
        "polyphen_score",
        "alphamissense_pathogenicity",
        "revel",
    ]:
        values = pd.to_numeric(combined[column], errors="coerce")
        if not values.dropna().between(0, 1).all():
            raise ValueError(f"{column} has values outside [0, 1].")

    if combined["variant_key"].duplicated().any():
        raise ValueError("Duplicate variant keys after integration.")

    combined.to_csv(OUTPUT_FILE, index=False)

    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "variant_count": len(combined),
        "gene_distribution": combined["GeneSymbol"].value_counts().to_dict(),
        "predictors": predictors,
        "predictor_summary": score_summary,
        "benchmark_eligible_count": int(
            combined["benchmark_eligible"].fillna(False).astype(bool).sum()
        ) if "benchmark_eligible" in combined.columns else None,
        "high_confidence_eligible_count": int(
            combined["high_confidence_eligible"].fillna(False).astype(bool).sum()
        ) if "high_confidence_eligible" in combined.columns else None,
        "output": str(OUTPUT_FILE.relative_to(ROOT)),
        "sha256": {
            "base": sha256(BASE_FILE),
            "vep": sha256(VEP_FILE),
            "benchmark": sha256(BENCHMARK_FILE),
            "output": sha256(OUTPUT_FILE),
        },
    }

    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)
    REPORT_FILE.write_text(
        json.dumps(report, indent=2, default=str),
        encoding="utf-8",
    )

    print("\nFIVE-PREDICTOR INTEGRATION COMPLETED")
    print("Variants:", len(combined))
    print("Predictors:", ", ".join(predictors))
    print("Output:", OUTPUT_FILE)
    print("Report:", REPORT_FILE)
    print("\nPredictor coverage:")
    for name, values in score_summary.items():
        print(
            f"{name}: {values['available']}/{len(combined)} "
            f"(missing={values['missing']})"
        )


if __name__ == "__main__":
    main()
