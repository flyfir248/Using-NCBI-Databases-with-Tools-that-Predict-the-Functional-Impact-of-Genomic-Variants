from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FILE = ROOT / "data/processed/brca1_brca2_primary_missense_cohort.csv"

df = pd.read_csv(FILE, low_memory=False)

print("\nCLINVAR CLASSIFICATIONS IN THE MISSENSE COHORT")
print("=" * 65)

counts = (
    df["ClinicalSignificance"]
    .fillna("<MISSING>")
    .astype(str)
    .value_counts(dropna=False)
)

for classification, count in counts.items():
    percentage = count / len(df) * 100
    print(f"{classification:<55} {count:>6} ({percentage:>6.2f}%)")

print("\nGENE-LEVEL CLASSIFICATION DISTRIBUTION")
print("=" * 65)

gene_counts = pd.crosstab(
    df["GeneSymbol"],
    df["ClinicalSignificance"].fillna("<MISSING>")
)

print(gene_counts.to_string())

print("\nREVIEW STATUS DISTRIBUTION")
print("=" * 65)

print(
    df["ReviewStatus"]
    .fillna("<MISSING>")
    .value_counts(dropna=False)
    .to_string()
)

print("\nTotal missense records:", len(df))
print("Unique AlleleIDs:", df["#AlleleID"].nunique())
