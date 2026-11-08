from pathlib import Path
import json

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

INPUT = (
    ROOT / "data/processed/"
    "brca1_brca2_primary_missense_cohort.csv"
)

OUTPUT = (
    ROOT / "data/metadata/"
    "vep_predictor_coverage_report.json"
)

df = pd.read_csv(INPUT, low_memory=False)

predictors = {
    "SIFT": ["sift_score", "sift_prediction"],
    "PolyPhen": ["polyphen_score", "polyphen_prediction"],
}

report = {
    "input_variants": len(df),
    "predictors": {},
}

for predictor, columns in predictors.items():

    print(f"\n{predictor}")
    print("=" * 65)

    predictor_report = {}

    for column in columns:

        if column not in df.columns:
            print(f"{column}: COLUMN MISSING")
            predictor_report[column] = {
                "column_present": False
            }
            continue

        values = df[column]
        available = values.notna() & values.astype(str).str.strip().ne("")

        print(f"Available: {available.sum()} / {len(df)}")
        print(f"Missing: {(~available).sum()}")

        predictor_report[column] = {
            "column_present": True,
            "available": int(available.sum()),
            "missing": int((~available).sum()),
            "coverage_percent": round(
                available.mean() * 100, 3
            ),
            "unique_values": int(values.nunique()),
        }

        if column.endswith("_score"):
            numeric = pd.to_numeric(values, errors="coerce")

            predictor_report[column].update({
                "numeric_available": int(numeric.notna().sum()),
                "minimum": (
                    float(numeric.min())
                    if numeric.notna().any() else None
                ),
                "maximum": (
                    float(numeric.max())
                    if numeric.notna().any() else None
                ),
                "median": (
                    float(numeric.median())
                    if numeric.notna().any() else None
                ),
            })

        else:
            predictor_report[column]["prediction_distribution"] = {
                str(k): int(v)
                for k, v in values.value_counts(
                    dropna=False
                ).items()
            }

    report["predictors"][predictor] = predictor_report

print("\nGENE-LEVEL SCORE COVERAGE")
print("=" * 65)

for gene, subset in df.groupby("GeneSymbol"):

    print(f"\n{gene}: {len(subset)} variants")

    for predictor, columns in predictors.items():

        score_column = columns[0]

        if score_column in subset.columns:
            available = pd.to_numeric(
                subset[score_column], errors="coerce"
            ).notna().sum()

            print(
                f"{predictor}: {available}/{len(subset)}"
            )

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

OUTPUT.write_text(
    json.dumps(report, indent=2),
    encoding="utf-8",
)

print("\nReport saved:", OUTPUT)