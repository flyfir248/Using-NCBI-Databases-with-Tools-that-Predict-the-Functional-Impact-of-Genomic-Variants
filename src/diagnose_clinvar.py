
import pandas as pd

FILE = "data/raw/clinvar_variant_summary.txt.gz"

columns = [
    "GeneSymbol",
    "Assembly",
    "Type",
    "Start",
    "Stop",
    "ReferenceAllele",
    "AlternateAllele",
]

print("Loading relevant ClinVar columns...")

df = pd.read_csv(
    FILE,
    sep="\t",
    compression="gzip",
    dtype=str,
    usecols=columns,
    keep_default_na=False,
)

print("\nTotal records:", len(df))

# Stage 1: Exact gene symbols
gene = df["GeneSymbol"].isin(["BRCA1", "BRCA2"])
print("Exact BRCA1/BRCA2:", gene.sum())

# Stage 2: Genome assembly
assembly = df["Assembly"].eq("GRCh38")
print("GRCh38:", (gene & assembly).sum())

# Stage 3: Variant type
snv = df["Type"].str.lower().eq("single nucleotide variant")
print("SNVs:", (gene & assembly & snv).sum())

# Stage 4: Valid alleles
ref = df["ReferenceAllele"].str.strip()
alt = df["AlternateAllele"].str.strip()

valid_bases = {"A", "C", "G", "T"}

valid_alleles = (
    ref.isin(valid_bases)
    & alt.isin(valid_bases)
    & ref.ne(alt)
)

print(
    "Valid SNV alleles:",
    (gene & assembly & snv & valid_alleles).sum()
)

# Stage 5: Matching genomic coordinates
start = pd.to_numeric(df["Start"], errors="coerce")
stop = pd.to_numeric(df["Stop"], errors="coerce")

valid_coordinates = (
    start.notna()
    & stop.notna()
    & start.eq(stop)
)

print(
    "Valid coordinates:",
    (
        gene
        & assembly
        & snv
        & valid_alleles
        & valid_coordinates
    ).sum()
)

# Inspect representative BRCA1/BRCA2 GRCh38 SNVs
sample = df.loc[
    gene & assembly & snv,
    [
        "GeneSymbol",
        "Start",
        "Stop",
        "ReferenceAllele",
        "AlternateAllele",
    ],
].head(20)

print("\nRepresentative records:\n")
print(sample.to_string(index=False))

print("\nAllele value distribution:\n")
print(df.loc[gene & assembly & snv, "ReferenceAllele"].value_counts().head(10))

print("\nCoordinate equality:\n")
print(
    df.loc[gene & assembly & snv, ["Start", "Stop"]]
    .head(20)
    .to_string(index=False)
)
