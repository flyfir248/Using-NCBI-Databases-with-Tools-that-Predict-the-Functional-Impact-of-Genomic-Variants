import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

ALL_FILE = ROOT / "data/processed/vep_primary_missense_classification.csv"
MISSENSE_FILE = ROOT / "data/processed/brca1_brca2_primary_missense_cohort.csv"
NON_MISSENSE_FILE = ROOT / "data/processed/brca1_brca2_primary_non_missense.csv"

all_variants = pd.read_csv(ALL_FILE, low_memory=False)
missense = pd.read_csv(MISSENSE_FILE, low_memory=False)
non_missense = pd.read_csv(NON_MISSENSE_FILE, low_memory=False)

checks = {}

# 1. Total dataset integrity
checks["total_candidate_count"] = len(all_variants) == 28302
checks["unique_allele_ids"] = all_variants["#AlleleID"].is_unique

# 2. Classification reconciliation
distribution = (
    all_variants["missense_classification"]
    .value_counts(dropna=False)
    .to_dict()
)

checks["missense_count"] = len(missense) == 16193
checks["non_missense_count"] = len(non_missense) == 12108
checks["manual_review_count"] = (
    all_variants["missense_classification"]
    .astype(str)
    .str.startswith("manual_review")
    .sum() == 1
)
checks["classification_counts_reconcile"] = (
    len(missense) + len(non_missense)
    + all_variants["missense_classification"]
      .astype(str).str.startswith("manual_review").sum()
    == len(all_variants)
)

# 3. Output membership
all_ids = set(all_variants["#AlleleID"].astype(str))
missense_ids = set(missense["#AlleleID"].astype(str))
non_missense_ids = set(non_missense["#AlleleID"].astype(str))

checks["missense_ids_match_classification"] = (
    missense_ids == set(
        all_variants.loc[
            all_variants["missense_classification"] == "missense",
            "#AlleleID"
        ].astype(str)
    )
)

checks["non_missense_ids_match_classification"] = (
    non_missense_ids == set(
        all_variants.loc[
            all_variants["missense_classification"] == "non_missense",
            "#AlleleID"
        ].astype(str)
    )
)

checks["no_overlap_between_output_cohorts"] = (
    missense_ids.isdisjoint(non_missense_ids)
)

checks["all_output_ids_exist_in_full_dataset"] = (
    missense_ids.union(non_missense_ids).issubset(all_ids)
)

# 4. Gene-level counts
gene_summary = {}

for gene in ["BRCA1", "BRCA2"]:
    subset = all_variants[all_variants["GeneSymbol"] == gene]
    gene_summary[gene] = {
        "total": int(len(subset)),
        "missense": int((subset["missense_classification"] == "missense").sum()),
        "non_missense": int((subset["missense_classification"] == "non_missense").sum()),
        "manual_review": int(
            subset["missense_classification"]
            .astype(str).str.startswith("manual_review").sum()
        ),
    }

checks["gene_counts_match_expected"] = gene_summary == {
    "BRCA1": {
        "total": 11714, "missense": 6374,
        "non_missense": 5340, "manual_review": 0
    },
    "BRCA2": {
        "total": 16588, "missense": 9819,
        "non_missense": 6768, "manual_review": 1
    },
}

# 5. Selected transcript gene consistency
selected = all_variants[
    all_variants["transcript_selection_method"] != "manual_review"
].copy()

checks["selected_transcripts_match_target_gene"] = (
    selected["GeneSymbol"].astype(str).str.upper()
    == selected["gene_symbol"].astype(str).str.upper()
).all()

# 6. Missense transcript biotype
checks["missense_transcripts_are_protein_coding"] = (
    missense["biotype"].astype(str) == "protein_coding"
).all()

# Report
report = {
    "total_variants": int(len(all_variants)),
    "missense_variants": int(len(missense)),
    "non_missense_variants": int(len(non_missense)),
    "classification_distribution": {
        str(k): int(v) for k, v in distribution.items()
    },
    "gene_summary": gene_summary,
    "checks": {k: bool(v) for k, v in checks.items()},
    "all_checks_passed": all(checks.values()),
}

report_path = ROOT / "data/metadata/primary_missense_integrity_audit.json"
report_path.parent.mkdir(parents=True, exist_ok=True)
report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

print("\nMISSENSE COHORT INTEGRITY AUDIT")
print("=" * 55)

for name, passed in checks.items():
    print(f"{'PASS' if passed else 'FAIL'} | {name}")

print("\nClassification distribution:")
for label, count in distribution.items():
    print(f"  {label}: {count}")

print("\nGene summary:")
for gene, values in gene_summary.items():
    print(f"  {gene}: {values}")

print(f"\nAudit report: {report_path}")

if not all(checks.values()):
    raise SystemExit("\nAUDIT FAILED — investigate failed checks.")

print("\nALL MISSENSE COHORT INTEGRITY CHECKS PASSED.")
