
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

COHORT_FILE = (
    ROOT / "data/processed/"
    "brca1_brca2_primary_missense_cohort.csv"
)

SOURCE_FILE = (
    ROOT / "data/raw/alphamissense/"
    "AlphaMissense_hg38.tsv.gz"
)

OUTPUT_DIR = ROOT / "data/processed"
METADATA_DIR = ROOT / "data/metadata"

MATCHED_FILE = (
    OUTPUT_DIR / "alphamissense_cohort_source_matches.csv"
)

UNMATCHED_FILE = (
    OUTPUT_DIR / "alphamissense_unmatched_cohort_variants.csv"
)

DUPLICATES_FILE = (
    OUTPUT_DIR / "alphamissense_duplicate_genomic_keys.csv"
)

REPORT_FILE = (
    METADATA_DIR / "alphamissense_extraction_report.json"
)

CHUNK_SIZE = 500_000

SOURCE_COLUMNS = [
    "#CHROM",
    "POS",
    "REF",
    "ALT",
    "genome",
    "uniprot_id",
    "transcript_id",
    "protein_variant",
    "am_pathogenicity",
    "am_class",
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

logger = logging.getLogger(__name__)


def normalize_chromosome(value):
    chromosome = str(value).strip()

    if chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]

    return chromosome.upper()


def build_source_key(df):
    chromosome = df["#CHROM"].map(normalize_chromosome)

    return (
        "GRCh38:"
        + chromosome
        + ":"
        + df["POS"].astype(str)
        + ":"
        + df["REF"].astype(str).str.upper()
        + ">"
        + df["ALT"].astype(str).str.upper()
    )


def main():

    if not COHORT_FILE.exists():
        raise FileNotFoundError(COHORT_FILE)

    if not SOURCE_FILE.exists():
        raise FileNotFoundError(SOURCE_FILE)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    METADATA_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("Loading BRCA1/BRCA2 missense cohort.")

    cohort = pd.read_csv(
        COHORT_FILE,
        low_memory=False,
        dtype={"variant_key": "string"},
    )

    if "variant_key" not in cohort.columns:
        raise ValueError("Cohort has no variant_key column.")

    if cohort["variant_key"].isna().any():
        raise ValueError("Cohort contains missing variant keys.")

    if cohort["variant_key"].duplicated().any():
        raise ValueError("Cohort variant keys are not unique.")

    target_keys = set(cohort["variant_key"].astype(str))

    logger.info("Cohort variants: %s", len(target_keys))

    matched_chunks = []
    matched_key_set = set()

    total_source_rows = 0
    matching_source_rows = 0
    matching_non_hg38_rows = 0

    logger.info("Scanning AlphaMissense source in chunks.")

    reader = pd.read_csv(
        SOURCE_FILE,
        sep="\t",
        skiprows=3,
        usecols=SOURCE_COLUMNS,
        chunksize=CHUNK_SIZE,
        compression="gzip",
        low_memory=False,
        dtype={
            "#CHROM": "string",
            "REF": "string",
            "ALT": "string",
            "genome": "string",
            "uniprot_id": "string",
            "transcript_id": "string",
            "protein_variant": "string",
            "am_class": "string",
        },
    )

    for chunk_number, chunk in enumerate(reader, start=1):

        total_source_rows += len(chunk)

        chunk["source_variant_key"] = build_source_key(chunk)

        matched = chunk[
            chunk["source_variant_key"].isin(target_keys)
        ].copy()

        if matched.empty:
            if chunk_number % 10 == 0:
                logger.info(
                    "Scanned %s source rows.",
                    f"{total_source_rows:,}",
                )
            continue

        matching_source_rows += len(matched)

        non_hg38 = (
            matched["genome"].str.lower() != "hg38"
        )

        matching_non_hg38_rows += int(non_hg38.sum())

        matched_key_set.update(
            matched["source_variant_key"].astype(str)
        )

        matched_chunks.append(matched)

        logger.info(
            "Chunk %s: found %s cohort records.",
            chunk_number,
            len(matched),
        )

    if matched_chunks:
        matches = pd.concat(
            matched_chunks,
            ignore_index=True,
        )
    else:
        matches = pd.DataFrame(columns=SOURCE_COLUMNS + [
            "source_variant_key"
        ])

    unmatched_keys = sorted(target_keys - matched_key_set)

    unmatched = cohort[
        cohort["variant_key"].astype(str).isin(unmatched_keys)
    ].copy()

    # Validate scores and categories without silently modifying them.
    if not matches.empty:

        matches["am_pathogenicity"] = pd.to_numeric(
            matches["am_pathogenicity"],
            errors="coerce",
        )

        invalid_scores = (
            matches["am_pathogenicity"].isna()
            | ~np.isfinite(matches["am_pathogenicity"])
            | ~matches["am_pathogenicity"].between(0, 1)
        )

        unexpected_classes = sorted(
            set(matches["am_class"].dropna().unique())
            - {
                "likely_benign",
                "ambiguous",
                "likely_pathogenic",
            }
        )

        duplicate_counts = (
            matches.groupby("source_variant_key")
            .size()
            .rename("source_record_count")
            .reset_index()
        )

        duplicate_keys = duplicate_counts[
            duplicate_counts["source_record_count"] > 1
        ].copy()

        class_counts = (
            matches["am_class"]
            .value_counts(dropna=False)
            .to_dict()
        )

        genome_counts = (
            matches["genome"]
            .value_counts(dropna=False)
            .to_dict()
        )

        score_summary = (
            matches["am_pathogenicity"]
            .describe()
            .to_dict()
        )

    else:

        invalid_scores = pd.Series(dtype=bool)
        unexpected_classes = []
        duplicate_keys = pd.DataFrame(
            columns=["source_variant_key", "source_record_count"]
        )
        class_counts = {}
        genome_counts = {}
        score_summary = {}

    matches.to_csv(MATCHED_FILE, index=False)
    unmatched.to_csv(UNMATCHED_FILE, index=False)
    duplicate_keys.to_csv(DUPLICATES_FILE, index=False)

    unique_matched_keys = len(matched_key_set)

    report = {
        "source": "AlphaMissense",
        "source_file": SOURCE_FILE.name,
        "source_release": "Zenodo record 10813168, version 3",
        "genome_assembly": "GRCh38",
        "retrieved_at_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "matching_method": (
            "Exact genomic identity: "
            "GRCh38:chromosome:position:REF>ALT"
        ),
        "chunk_size": CHUNK_SIZE,
        "cohort_variant_count": len(target_keys),
        "source_rows_scanned": total_source_rows,
        "matching_source_records": matching_source_rows,
        "unique_matched_genomic_keys": unique_matched_keys,
        "unmatched_cohort_variants": len(unmatched_keys),
        "source_records_with_non_hg38_genome": (
            matching_non_hg38_rows
        ),
        "duplicate_genomic_keys": len(duplicate_keys),
        "duplicate_key_records": (
            int(duplicate_keys["source_record_count"].sum())
            if not duplicate_keys.empty else 0
        ),
        "genome_distribution": {
            str(k): int(v) for k, v in genome_counts.items()
        },
        "class_distribution": {
            str(k): int(v) for k, v in class_counts.items()
        },
        "score_summary": {
            str(k): float(v) for k, v in score_summary.items()
            if pd.notna(v)
        },
        "invalid_score_records": int(invalid_scores.sum()),
        "unexpected_classes": unexpected_classes,
        "outputs": {
            "matched_records": str(MATCHED_FILE.relative_to(ROOT)),
            "unmatched_variants": str(UNMATCHED_FILE.relative_to(ROOT)),
            "duplicate_keys": str(DUPLICATES_FILE.relative_to(ROOT)),
        },
    }

    REPORT_FILE.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )

    logger.info("EXTRACTION COMPLETED")
    logger.info("Source rows scanned: %s", f"{total_source_rows:,}")
    logger.info("Matched source records: %s", f"{matching_source_rows:,}")
    logger.info("Unique matched variants: %s", f"{unique_matched_keys:,}")
    logger.info("Unmatched cohort variants: %s", f"{len(unmatched_keys):,}")
    logger.info("Duplicate genomic keys: %s", len(duplicate_keys))
    logger.info("Invalid scores: %s", int(invalid_scores.sum()))
    logger.info("Unexpected classes: %s", unexpected_classes)
    logger.info("Report: %s", REPORT_FILE)

    if matching_non_hg38_rows:
        raise RuntimeError(
            "Matching records contain non-hg38 genome annotations."
        )

    if invalid_scores.any():
        raise RuntimeError(
            "Invalid AlphaMissense scores detected."
        )

    if unexpected_classes:
        raise RuntimeError(
            "Unexpected AlphaMissense classes detected."
        )


if __name__ == "__main__":
    main()
