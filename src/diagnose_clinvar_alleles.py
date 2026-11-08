
import pandas as pd

FILE = "data/raw/clinvar_variant_summary.txt.gz"

columns = [
    "GeneSymbol",
    "Assembly",
    "Type",
    "VariationID",
    "Chromosome",
    "Start",
    "Stop",
    "ReferenceAllele",
    "AlternateAllele",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
]

df = pd.read_csv(
    FILE,
    sep="\t",
    compression="gzip",
    dtype=str,
    usecols=columns,
    keep_default_na=False,
)

mask = (
    df["GeneSymbol"].isin(["BRCA1", "BRCA2"])
    & df["Assembly"].eq("GRCh38")
    & df["Type"].eq("single nucleotide variant")
)

snvs = df.loc[mask].copy()

print("Total BRCA1/BRCA2 GRCh38 SNVs:", len(snvs))

fields = [
    "ReferenceAllele",
    "AlternateAllele",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
]

print("\nMissing-value distributions:")

for field in fields:
    print(f"\n{field}:")
    print(snvs[field].value_counts().head(10))

print("\nRepresentative VCF records:")

print(
    snvs[
        [
            "GeneSymbol",
            "VariationID",
            "Chromosome",
            "Start",
            "PositionVCF",
            "ReferenceAlleleVCF",
            "AlternateAlleleVCF",
        ]
    ].head(20).to_string(index=False)
)

print("\nUnique VCF reference alleles:")
print(snvs["ReferenceAlleleVCF"].value_counts().head(10))

print("\nUnique VCF alternate alleles:")
print(snvs["AlternateAlleleVCF"].value_counts().head(10))
