from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FILE = ROOT / "data/processed/brca1_brca2_primary_missense_cohort.csv"

df = pd.read_csv(FILE, low_memory=False)

positive = [
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
]

negative = [
    "Benign",
    "Likely benign",
    "Benign/Likely benign",
]

df["benchmark_group"] = "excluded"

df.loc[
    df["ClinicalSignificance"].isin(positive),
    "benchmark_group"
] = "pathogenic"

df.loc[
    df["ClinicalSignificance"].isin(negative),
    "benchmark_group"
] = "benign"

print("\nBENCHMARK LABELS BY REVIEW STATUS")
print("=" * 85)

table = pd.crosstab(
    df["ReviewStatus"],
    df["benchmark_group"],
    margins=True
)

print(table.to_string())

print("\nBENIGN/PATHOGENIC RECORDS BY GENE AND REVIEW STATUS")
print("=" * 85)

eligible = df[df["benchmark_group"] != "excluded"]

gene_review = pd.crosstab(
    [eligible["GeneSymbol"], eligible["ReviewStatus"]],
    eligible["benchmark_group"],
    margins=False
)

print(gene_review.to_string())

print("\nHIGH-CONFIDENCE SUBSET COUNTS")
print("=" * 85)

high_confidence_statuses = [
    "reviewed by expert panel",
    "criteria provided, multiple submitters, no conflicts",
]

high_confidence = eligible[
    eligible["ReviewStatus"].isin(high_confidence_statuses)
]

print("High-confidence total:", len(high_confidence))
print("Pathogenic:", (high_confidence["benchmark_group"] == "pathogenic").sum())
print("Benign:", (high_confidence["benchmark_group"] == "benign").sum())

print("\nHIGH-CONFIDENCE COUNTS BY GENE")
print(
    pd.crosstab(
        high_confidence["GeneSymbol"],
        high_confidence["benchmark_group"]
    ).to_string()
)
