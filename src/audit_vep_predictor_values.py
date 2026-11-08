import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

INPUT = (
    ROOT / "data/processed/"
    "brca1_brca2_primary_missense_cohort.csv"
)

OUTPUT = (
    ROOT / "data/metadata/"
    "vep_predictor_value_audit.json"
)

PREDICTORS = {
    "SIFT": {
        "score": "sift_score",
        "prediction": "sift_prediction",
        "expected_categories": {
            "deleterious",
            "tolerated",
            "deleterious_low_confidence",
            "tolerated_low_confidence",
        },
    },
    "PolyPhen-2": {
        "score": "polyphen_score",
        "prediction": "polyphen_prediction",
        "expected_categories": {
            "benign",
            "possibly_damaging",
            "probably_damaging",
            "unknown",
        },
    },
}


def main():

    df = pd.read_csv(INPUT, low_memory=False)

    if len(df) != 16193:
        raise ValueError(f"Unexpected cohort size: {len(df)}")

    report = {
        "input_file": str(INPUT),
        "cohort_size": len(df),
        "predictors": {},
        "all_checks_passed": True,
    }

    for name, config in PREDICTORS.items():

        score_col = config["score"]
        prediction_col = config["prediction"]

        scores = pd.to_numeric(
            df[score_col], errors="coerce"
        )

        predictions = (
            df[prediction_col]
            .fillna("<MISSING>")
            .astype(str)
            .str.strip()
        )

        category_counts = {
            str(k): int(v)
            for k, v in predictions.value_counts().items()
        }

        unexpected_categories = sorted(
            set(predictions.unique())
            - config["expected_categories"]
        )

        gene_coverage = {}

        for gene, subset in df.groupby("GeneSymbol"):

            gene_scores = pd.to_numeric(
                subset[score_col], errors="coerce"
            )

            gene_coverage[gene] = {
                "total": int(len(subset)),
                "available": int(gene_scores.notna().sum()),
                "missing": int(gene_scores.isna().sum()),
            }

        checks = {
            "all_scores_present": bool(scores.notna().all()),
            "all_scores_finite": bool(
                np.isfinite(scores.dropna()).all()
            ),
            "scores_within_zero_and_one": bool(
                scores.dropna().between(0, 1).all()
            ),
            "all_predictions_present": bool(
                predictions.ne("<MISSING>").all()
            ),
            "no_unexpected_prediction_categories": (
                len(unexpected_categories) == 0
            ),
            "prediction_count_matches_cohort": (
                len(predictions) == len(df)
            ),
        }

        report["predictors"][name] = {
            "score_column": score_col,
            "prediction_column": prediction_col,
            "score_statistics": {
                "count": int(scores.count()),
                "minimum": float(scores.min()),
                "median": float(scores.median()),
                "maximum": float(scores.max()),
            },
            "prediction_distribution": category_counts,
            "unexpected_categories": unexpected_categories,
            "gene_coverage": gene_coverage,
            "checks": checks,
            "all_checks_passed": all(checks.values()),
        }

        if not all(checks.values()):
            report["all_checks_passed"] = False

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    OUTPUT.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(report, indent=2))

    if not report["all_checks_passed"]:
        raise SystemExit("PREDICTOR VALUE AUDIT FAILED.")

    print("\nALL PREDICTOR VALUE CHECKS PASSED.")


if __name__ == "__main__":
    main()
