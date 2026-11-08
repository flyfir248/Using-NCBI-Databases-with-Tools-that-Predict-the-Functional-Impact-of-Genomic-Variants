
"""
ClinVar data acquisition module.

Research:
Integrative Assessment of Genomic Variant Functional Impact
Using NCBI Resources and Computational Prediction Algorithms.

Scope:
- Human BRCA1 and BRCA2 variants
- GRCh38 assembly
- Single nucleotide variants
- Clinical evidence retained without premature classification filtering

Source:
NCBI ClinVar variant_summary.txt.gz
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm


PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DIR = PROJECT_ROOT / "data" / "raw"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
METADATA_DIR = PROJECT_ROOT / "data" / "metadata"

CLINVAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/"
    "tab_delimited/variant_summary.txt.gz"
)

RAW_FILE = RAW_DIR / "clinvar_variant_summary.txt.gz"

OUTPUT_FILE = (
    PROCESSED_DIR /
    "brca1_brca2_grch38_candidate_snvs.csv"
)

REPORT_FILE = (
    METADATA_DIR /
    "clinvar_extraction_report.json"
)

MANIFEST_FILE = (
    METADATA_DIR /
    "clinvar_download_manifest.json"
)

TARGET_GENES = ["BRCA1", "BRCA2"]

CHUNK_SIZE = 250_000


REQUIRED_COLUMNS = [
    "GeneSymbol",
    "Assembly",
    "Type",
    "VariationID",
    "ClinicalSignificance",
    "ReviewStatus",
    "Chromosome",
    "ChromosomeAccession",
    "Start",
    "Stop",
    "ReferenceAllele",
    "AlternateAllele",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
]



def setup_directories() -> None:
    """Create project directories."""

    for directory in [
        RAW_DIR,
        PROCESSED_DIR,
        METADATA_DIR,
    ]:
        directory.mkdir(parents=True, exist_ok=True)


def configure_logging() -> None:
    """Configure logging."""

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def calculate_sha256(filepath: Path) -> str:
    """Calculate SHA256 checksum."""

    digest = hashlib.sha256()

    with filepath.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def download_clinvar(force: bool = False) -> dict:
    """
    Download the official ClinVar variant summary.

    Existing downloads are reused unless force=True.
    """

    if RAW_FILE.exists() and not force:
        logging.info("Existing ClinVar archive found. Reusing it.")

    else:
        logging.info("Downloading ClinVar archive...")

        temporary_file = RAW_FILE.with_suffix(
            RAW_FILE.suffix + ".part"
        )

        try:
            with requests.get(
                CLINVAR_URL,
                stream=True,
                timeout=(30, 300),
                headers={
                    "User-Agent": "NCBI-Functional-Variant-Research/1.0"
                },
            ) as response:

                response.raise_for_status()

                total_size = int(
                    response.headers.get("Content-Length", 0)
                )

                with temporary_file.open("wb") as file:

                    with tqdm(
                        total=total_size or None,
                        unit="B",
                        unit_scale=True,
                        desc="ClinVar download",
                    ) as progress:

                        for chunk in response.iter_content(
                            chunk_size=1024 * 1024
                        ):
                            if chunk:
                                file.write(chunk)
                                progress.update(len(chunk))

            temporary_file.replace(RAW_FILE)

        except Exception:
            temporary_file.unlink(missing_ok=True)
            raise

    checksum = calculate_sha256(RAW_FILE)

    manifest = {
        "source": "NCBI ClinVar",
        "source_url": CLINVAR_URL,
        "retrieved_at_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "filename": RAW_FILE.name,
        "file_size_bytes": RAW_FILE.stat().st_size,
        "sha256": checksum,
    }

    with MANIFEST_FILE.open("w", encoding="utf-8") as file:
        json.dump(manifest, file, indent=4)

    logging.info("ClinVar archive ready.")
    logging.info("SHA256: %s", checksum)

    return manifest


def is_single_nucleotide_variant(
    dataframe: pd.DataFrame,
) -> pd.Series:
    """Identify simple single-base substitutions."""

    ref = dataframe["ReferenceAllele"].str.strip()
    alt = dataframe["AlternateAllele"].str.strip()

    valid_bases = {"A", "C", "G", "T"}

    return (
        dataframe["Type"].str.lower()
        .eq("single nucleotide variant")
        & ref.isin(valid_bases)
        & alt.isin(valid_bases)
        & ref.ne(alt)
        & dataframe["Start"].eq(dataframe["Stop"])
    )



def extract_variants() -> tuple[pd.DataFrame, dict]:
    """
    Extract canonical BRCA1/BRCA2 GRCh38 candidate SNVs.

    ClinVar VCF-specific fields are used for genomic alleles.
    Original HGVS-style fields are preserved for provenance.
    """

    logging.info("Reading ClinVar variant summary...")

    selected_chunks = []

    total_rows = 0
    gene_rows = 0
    assembly_rows = 0
    snv_rows = 0
    canonical_snv_rows = 0
    ambiguous_allele_rows = 0

    for chunk in pd.read_csv(
        RAW_FILE,
        sep="\t",
        compression="gzip",
        dtype=str,
        keep_default_na=False,
        chunksize=CHUNK_SIZE,
        low_memory=False,
    ):

        total_rows += len(chunk)

        missing = set(REQUIRED_COLUMNS) - set(chunk.columns)

        if missing:
            raise ValueError(
                f"Missing ClinVar columns: {missing}"
            )

        # Stage 1: Exact target gene symbols.
        gene_mask = chunk["GeneSymbol"].isin(TARGET_GENES)
        gene_chunk = chunk.loc[gene_mask].copy()

        gene_rows += len(gene_chunk)

        # Stage 2: GRCh38.
        assembly_mask = gene_chunk["Assembly"].eq("GRCh38")
        assembly_chunk = gene_chunk.loc[assembly_mask].copy()

        assembly_rows += len(assembly_chunk)

        # Stage 3: Single nucleotide variants.
        snv_mask = assembly_chunk["Type"].eq(
            "single nucleotide variant"
        )

        snv_chunk = assembly_chunk.loc[snv_mask].copy()

        snv_rows += len(snv_chunk)

        if snv_chunk.empty:
            continue

        # Stage 4: Use VCF-specific allele fields.
        ref = snv_chunk["ReferenceAlleleVCF"].str.strip().str.upper()
        alt = snv_chunk["AlternateAlleleVCF"].str.strip().str.upper()

        valid_bases = {"A", "C", "G", "T"}

        valid_alleles = (
            ref.isin(valid_bases)
            & alt.isin(valid_bases)
            & ref.ne(alt)
        )

        ambiguous_allele_rows += int((~valid_alleles).sum())

        canonical_chunk = snv_chunk.loc[valid_alleles].copy()

        canonical_chunk["ReferenceAlleleVCF"] = ref[valid_alleles]
        canonical_chunk["AlternateAlleleVCF"] = alt[valid_alleles]

        # Stage 5: Validate VCF position.
        canonical_chunk["PositionVCF"] = pd.to_numeric(
            canonical_chunk["PositionVCF"],
            errors="coerce",
        )

        canonical_chunk = canonical_chunk.dropna(
            subset=["PositionVCF"]
        )

        canonical_chunk = canonical_chunk.loc[
            canonical_chunk["PositionVCF"] > 0
        ].copy()

        canonical_snv_rows += len(canonical_chunk)

        if not canonical_chunk.empty:
            selected_chunks.append(canonical_chunk)

    if not selected_chunks:
        raise RuntimeError(
            "No canonical BRCA1/BRCA2 GRCh38 SNVs found."
        )

    variants = pd.concat(
        selected_chunks,
        ignore_index=True,
    )

    variants["PositionVCF"] = variants[
        "PositionVCF"
    ].astype("int64")

    # Build a genomic allele key.
    variants["variant_key"] = (
        variants["Assembly"].astype(str)
        + ":"
        + variants["Chromosome"].astype(str)
        + ":"
        + variants["PositionVCF"].astype(str)
        + ":"
        + variants["ReferenceAlleleVCF"].astype(str)
        + ">"
        + variants["AlternateAlleleVCF"].astype(str)
    )

    # Retain distinct ClinVar records while removing exact duplicates.
    deduplication_columns = [
        "VariationID",
        "GeneSymbol",
        "Assembly",
        "ChromosomeAccession",
        "PositionVCF",
        "ReferenceAlleleVCF",
        "AlternateAlleleVCF",
    ]

    before_deduplication = len(variants)

    variants = variants.drop_duplicates(
        subset=deduplication_columns
    ).copy()

    duplicates_removed = (
        before_deduplication - len(variants)
    )

    # Candidate labels only; do not treat these as final benchmark labels.
    significance = (
        variants["ClinicalSignificance"]
        .str.strip()
        .str.lower()
    )

    variants["binary_clinvar_label_candidate"] = pd.NA

    positive = significance.isin([
        "pathogenic",
        "likely pathogenic",
    ])

    negative = significance.isin([
        "benign",
        "likely benign",
    ])

    variants.loc[
        positive,
        "binary_clinvar_label_candidate"
    ] = 1

    variants.loc[
        negative,
        "binary_clinvar_label_candidate"
    ] = 0

    report = {
        "source": "NCBI ClinVar",
        "retrieved_at_utc": datetime.now(
            timezone.utc
        ).isoformat(),

        "assembly": "GRCh38",
        "target_genes": TARGET_GENES,

        "total_clinvar_rows_scanned": total_rows,
        "target_gene_rows": gene_rows,
        "grch38_rows": assembly_rows,
        "candidate_snv_rows": snv_rows,

        "canonical_snv_rows_before_deduplication":
            canonical_snv_rows,

        "ambiguous_or_invalid_allele_rows":
            ambiguous_allele_rows,

        "duplicates_removed": duplicates_removed,

        "final_unique_candidate_variants": len(variants),

        "gene_distribution": (
            variants["GeneSymbol"]
            .value_counts()
            .to_dict()
        ),

        "clinical_significance_distribution": (
            variants["ClinicalSignificance"]
            .value_counts()
            .to_dict()
        ),

        "review_status_distribution": (
            variants["ReviewStatus"]
            .value_counts()
            .to_dict()
        ),

        "label_candidate_distribution": (
            variants["binary_clinvar_label_candidate"]
            .value_counts(dropna=False)
            .astype(str)
            .to_dict()
        ),

        "coordinate_representation": (
            "PositionVCF + ReferenceAlleleVCF + "
            "AlternateAlleleVCF"
        ),

        "important_note": (
            "These are canonical candidate SNVs, not confirmed "
            "missense variants. Transcript-level annotation is "
            "required. Clinical classifications are provisional "
            "and must be evaluated at the appropriate evidence level."
        ),
    }

    return variants, report



def save_results(
    variants: pd.DataFrame,
    report: dict,
) -> None:
    """Save processed variants and extraction report."""

    variants.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    with REPORT_FILE.open(
        "w",
        encoding="utf-8",
    ) as file:

        json.dump(
            report,
            file,
            indent=4,
            default=str,
        )

    logging.info("Processed variants saved: %s", OUTPUT_FILE)
    logging.info("Extraction report saved: %s", REPORT_FILE)


def main() -> None:

    parser = argparse.ArgumentParser(
        description="Download and extract ClinVar BRCA variants."
    )

    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Download ClinVar again even if archive exists.",
    )

    args = parser.parse_args()

    setup_directories()
    configure_logging()

    download_clinvar(force=args.force_download)

    variants, report = extract_variants()

    save_results(variants, report)

    logging.info("=" * 55)
    logging.info("CLINVAR EXTRACTION COMPLETED")
    logging.info("=" * 55)

    logging.info(
        "Final candidate variants: %s",
        len(variants),
    )

    logging.info(
        "BRCA1: %s",
        report["gene_distribution"].get("BRCA1", 0),
    )

    logging.info(
        "BRCA2: %s",
        report["gene_distribution"].get("BRCA2", 0),
    )


if __name__ == "__main__":
    main()
