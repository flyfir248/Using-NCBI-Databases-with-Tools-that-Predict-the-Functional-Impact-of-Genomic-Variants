import hashlib
import json
import logging
import math
import shutil
import time
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path

import pandas as pd
import requests


# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

INPUT_CSV = (
    ROOT / "data/processed/"
    "brca1_brca2_primary_missense_cohort.csv"
)

RUN_DIR = ROOT / "data/raw/vep/predictor_full_run"
BATCH_DIR = RUN_DIR / "batches"

PROCESSED_DIR = ROOT / "data/processed"
METADATA_DIR = ROOT / "data/metadata"

API_URL = "https://rest.ensembl.org/vep/homo_sapiens/region"

ASSEMBLY = "GRCh38"
BATCH_SIZE = 200
MAX_ATTEMPTS = 5
REQUEST_TIMEOUT = (15, 180)
PAUSE_BETWEEN_BATCHES = 1.0

PARAMETERS = {
    "mane": 1,
    "hgvs": 1,
    "canonical": 1,
    "protein": 1,
    "numbers": 1,
    "tsl": 1,
    "CADD": 1,
    "REVEL": 1,
}

REQUIRED_COLUMNS = [
    "#AlleleID",
    "GeneSymbol",
    "Chromosome",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
    "variant_key",
]

RETRY_STATUS_CODES = {408, 429, 500, 502, 503, 504}

OUTPUT_CSV = (
    PROCESSED_DIR / "vep_full_predictor_annotations.csv"
)

FAILURE_CSV = (
    PROCESSED_DIR / "vep_full_predictor_failures.csv"
)

REPORT_JSON = (
    METADATA_DIR / "vep_full_predictor_report.json"
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)


# ============================================================
# GENERAL UTILITIES
# ============================================================

def utc_now():
    return datetime.now(timezone.utc).isoformat()


def write_json_atomic(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)

    temporary_path = path.with_suffix(path.suffix + ".tmp")

    with open(temporary_path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, allow_nan=False)

    temporary_path.replace(path)


def sha256_file(path):
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def normalize_chromosome(value):
    chromosome = str(value).strip()

    if chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]

    if chromosome.endswith(".0"):
        chromosome = chromosome[:-2]

    return chromosome


def make_identity_key(chromosome, position, ref, alt):
    return (
        f"{normalize_chromosome(chromosome)}:"
        f"{int(position)}:"
        f"{str(ref).strip().upper()}:"
        f"{str(alt).strip().upper()}"
    )


def make_vep_input(row):
    chromosome = normalize_chromosome(row["Chromosome"])
    position = int(row["PositionVCF"])

    ref = str(row["ReferenceAlleleVCF"]).strip().upper()
    alt = str(row["AlternateAlleleVCF"]).strip().upper()

    if position <= 0:
        raise ValueError(f"Invalid genomic position: {position}")

    if ref not in {"A", "C", "G", "T"}:
        raise ValueError(f"Invalid reference allele: {ref}")

    if alt not in {"A", "C", "G", "T"}:
        raise ValueError(f"Invalid alternate allele: {alt}")

    if ref == alt:
        raise ValueError(f"Reference equals alternate allele: {ref}")

    # VEP region format: chromosome start ID REF ALT QUAL FILTER INFO
    return f"{chromosome} {position} . {ref} {alt} . . ."


def identity_key_from_vep_input(value):
    """
    Support both the full VEP region input representation and
    the compact representation returned by some VEP responses.
    """
    fields = str(value).split()

    if len(fields) < 5:
        raise ValueError(
            f"Unexpected VEP input representation: {value}"
        )

    chromosome = fields[0]
    position = int(fields[1])

    if "/" in fields[3]:
        ref, alt = fields[3].split("/", 1)
    else:
        ref, alt = fields[3], fields[4]

    return make_identity_key(
        chromosome,
        position,
        ref,
        alt,
    )


def identity_key_from_result(result):
    value = result.get("input")

    if not value:
        raise ValueError("VEP result has no input field.")

    return identity_key_from_vep_input(value)


def prepare_directories():
    for directory in [
        RUN_DIR,
        BATCH_DIR,
        PROCESSED_DIR,
        METADATA_DIR,
    ]:
        directory.mkdir(parents=True, exist_ok=True)


def safe_numeric_score(value, field_name, identity):
    if value is None:
        return None

    try:
        score = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Non-numeric {field_name} for {identity}: {value}"
        ) from exc

    if not math.isfinite(score):
        raise ValueError(
            f"Non-finite {field_name} for {identity}: {value}"
        )

    if field_name == "revel" and not 0.0 <= score <= 1.0:
        raise ValueError(
            f"REVEL outside [0,1] for {identity}: {score}"
        )

    return score


# ============================================================
# INPUT DATA AND RUN CONFIGURATION
# ============================================================

def load_candidates():
    logging.info("Loading primary missense cohort...")

    df = pd.read_csv(INPUT_CSV, low_memory=False)

    missing = [
        column for column in REQUIRED_COLUMNS
        if column not in df.columns
    ]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if df.empty:
        raise ValueError("Input cohort is empty.")

    if df[REQUIRED_COLUMNS].isna().any().any():
        raise ValueError(
            "Required source fields contain missing values."
        )

    # Recalculate the source row index from the current input order.
    # The input cohort may already contain a source_row_index column.
    if "source_row_index" in df.columns:
        df = df.drop(columns=["source_row_index"])

    df.insert(0, "source_row_index", range(len(df)))

    df["vep_identity_key"] = df.apply(
        lambda row: make_identity_key(
            row["Chromosome"],
            row["PositionVCF"],
            row["ReferenceAlleleVCF"],
            row["AlternateAlleleVCF"],
        ),
        axis=1,
    )

    df["vep_input"] = df.apply(make_vep_input, axis=1)

    # Reconstruct the cohort identity format:
    # GRCh38:chromosome:position:REF>ALT
    reconstructed_keys = df.apply(
        lambda row: (
            "GRCh38:"
            + normalize_chromosome(row["Chromosome"])
            + ":"
            + str(int(row["PositionVCF"]))
            + ":"
            + str(row["ReferenceAlleleVCF"]).upper()
            + ">"
            + str(row["AlternateAlleleVCF"]).upper()
        ),
        axis=1,
    )

    if not (
        reconstructed_keys.astype(str).values
        == df["variant_key"].astype(str).values
    ).all():
        mismatches = df.loc[
            reconstructed_keys.astype(str)
            != df["variant_key"].astype(str),
            [
                "#AlleleID",
                "variant_key",
                "vep_identity_key",
            ],
        ]

        raise ValueError(
            "Genomic identity does not match variant_key:\n"
            + mismatches.head(10).to_string(index=False)
        )

    if df["vep_identity_key"].duplicated().any():
        raise ValueError(
            "Duplicate genomic identities in the input cohort."
        )

    if df["#AlleleID"].duplicated().any():
        raise ValueError(
            "Duplicate ClinVar AlleleIDs in the input cohort."
        )

    if set(df["GeneSymbol"].unique()) != {"BRCA1", "BRCA2"}:
        raise ValueError(
            "Unexpected gene composition in input cohort."
        )

    logging.info("Loaded %s unique missense variants.", len(df))
    logging.info(
        "Gene distribution: %s",
        df["GeneSymbol"].value_counts().to_dict(),
    )

    return df


def build_run_configuration(df):
    return {
        "source_file": str(INPUT_CSV.relative_to(ROOT)),
        "source_sha256": sha256_file(INPUT_CSV),
        "input_variant_count": int(len(df)),
        "assembly": ASSEMBLY,
        "api_url": API_URL,
        "parameters": PARAMETERS,
        "batch_size": BATCH_SIZE,
        "transcript_selection": (
            "Exactly one target-gene MANE Select transcript "
            "per genomic variant."
        ),
        "predictor_fields": [
            "cadd_raw",
            "cadd_phred",
            "revel",
        ],
    }


def initialize_manifest(configuration):
    manifest_path = RUN_DIR / "run_manifest.json"

    if manifest_path.exists():
        with open(manifest_path, encoding="utf-8") as handle:
            existing = json.load(handle)

        for key, value in configuration.items():
            if existing.get(key) != value:
                raise RuntimeError(
                    f"Existing run configuration differs for '{key}'. "
                    "Do not resume with changed inputs or settings."
                )

        logging.info("Existing run configuration verified.")
        return existing

    manifest = {
        **configuration,
        "created_at_utc": utc_now(),
        "status": "initialized",
    }

    write_json_atomic(manifest_path, manifest)
    return manifest


def load_checkpoint():
    checkpoint_path = RUN_DIR / "checkpoint.json"

    if checkpoint_path.exists():
        with open(checkpoint_path, encoding="utf-8") as handle:
            checkpoint = json.load(handle)

        if not isinstance(checkpoint.get("batches"), dict):
            raise RuntimeError("Invalid checkpoint structure.")

        return checkpoint

    return {
        "created_at_utc": utc_now(),
        "updated_at_utc": utc_now(),
        "batches": {},
    }


def save_checkpoint(checkpoint):
    checkpoint["updated_at_utc"] = utc_now()

    write_json_atomic(
        RUN_DIR / "checkpoint.json",
        checkpoint,
    )


# ============================================================
# VEP REQUEST AND RETRY HANDLING
# ============================================================

def retry_delay(response, attempt):
    retry_after = response.headers.get("Retry-After")

    if retry_after:
        try:
            return min(float(retry_after), 120.0)
        except ValueError:
            try:
                retry_time = parsedate_to_datetime(retry_after)
                seconds = (
                    retry_time - datetime.now(timezone.utc)
                ).total_seconds()

                return min(max(seconds, 1.0), 120.0)
            except (TypeError, ValueError):
                pass

    return min(2 ** attempt, 60)


def submit_to_vep(session, variants):
    payload = {"variants": variants}

    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "BRCA-Functional-Variant-Analysis/2.0",
    }

    for attempt in range(1, MAX_ATTEMPTS + 1):
        logging.info(
            "VEP request attempt %s/%s for %s variants.",
            attempt,
            MAX_ATTEMPTS,
            len(variants),
        )

        try:
            response = session.post(
                API_URL,
                params=PARAMETERS,
                headers=headers,
                json=payload,
                timeout=REQUEST_TIMEOUT,
            )

        except (requests.Timeout, requests.ConnectionError) as exc:
            if attempt == MAX_ATTEMPTS:
                raise RuntimeError(
                    f"Network failure after {MAX_ATTEMPTS} attempts: {exc}"
                ) from exc

            wait = min(2 ** attempt, 60)
            logging.warning(
                "Network error: %s. Retrying in %s seconds.",
                exc,
                wait,
            )
            time.sleep(wait)
            continue

        if response.status_code in RETRY_STATUS_CODES:
            if attempt == MAX_ATTEMPTS:
                response.raise_for_status()

            wait = retry_delay(response, attempt)
            logging.warning(
                "HTTP %s. Retrying in %.1f seconds.",
                response.status_code,
                wait,
            )
            time.sleep(wait)
            continue

        response.raise_for_status()

        try:
            results = response.json()
        except ValueError as exc:
            raise RuntimeError(
                "VEP returned invalid JSON."
            ) from exc

        if not isinstance(results, list):
            raise RuntimeError(
                "Unexpected VEP response: expected a JSON list."
            )

        return results

    raise RuntimeError("VEP request failed after all attempts.")


# ============================================================
# STRICT RESPONSE VALIDATION
# ============================================================

def validate_response(expected_rows, results):
    if not isinstance(results, list):
        raise ValueError("VEP response is not a list.")

    expected_map = {
        row["vep_identity_key"]: row
        for _, row in expected_rows.iterrows()
    }

    returned_map = {}

    for result in results:
        if not isinstance(result, dict):
            raise ValueError(
                "VEP response contains a non-object result."
            )

        key = identity_key_from_result(result)

        if key in returned_map:
            raise ValueError(f"Duplicate VEP result: {key}")

        if key not in expected_map:
            raise ValueError(f"Unexpected VEP result: {key}")

        returned_map[key] = result

    missing = set(expected_map) - set(returned_map)

    if missing:
        raise ValueError(
            f"{len(missing)} submitted variants missing from response. "
            f"Examples: {list(missing)[:5]}"
        )

    if len(results) != len(expected_rows):
        raise ValueError(
            f"Expected {len(expected_rows)} results; "
            f"received {len(results)}."
        )

    for key, result in returned_map.items():
        source = expected_map[key]

        if result.get("assembly_name") != ASSEMBLY:
            raise ValueError(
                f"Assembly mismatch for {key}: "
                f"{result.get('assembly_name')}"
            )

        chromosome = normalize_chromosome(
            result.get("seq_region_name")
        )

        expected_chromosome = normalize_chromosome(
            source["Chromosome"]
        )

        if chromosome != expected_chromosome:
            raise ValueError(f"Chromosome mismatch for {key}")

        if int(result.get("start")) != int(source["PositionVCF"]):
            raise ValueError(f"Position mismatch for {key}")

        expected_alleles = (
            str(source["ReferenceAlleleVCF"]).upper()
            + "/"
            + str(source["AlternateAlleleVCF"]).upper()
        )

        returned_alleles = str(
            result.get("allele_string", "")
        ).upper()

        if returned_alleles != expected_alleles:
            raise ValueError(
                f"Allele mismatch for {key}: "
                f"expected {expected_alleles}; "
                f"received {returned_alleles}"
            )

    return returned_map


# ============================================================
# BATCH CHECKPOINT AND RAW DATA MANAGEMENT
# ============================================================

def batch_paths(batch_number):
    directory = BATCH_DIR / f"batch_{batch_number:04d}"

    return {
        "directory": directory,
        "request": directory / "request.json",
        "response": directory / "response.json",
        "error": directory / "error.json",
    }


def archive_invalid_response(response_path):
    timestamp = datetime.now(timezone.utc).strftime(
        "%Y%m%dT%H%M%S"
    )

    destination = response_path.parent / (
        f"invalid_response_{timestamp}.json"
    )

    shutil.move(str(response_path), str(destination))

    logging.warning(
        "Archived invalid response as %s",
        destination,
    )


def process_batch(batch_number, batch_df, session, checkpoint):
    paths = batch_paths(batch_number)
    paths["directory"].mkdir(parents=True, exist_ok=True)

    expected_inputs = batch_df["vep_input"].tolist()

    request_data = {
        "batch_number": batch_number,
        "created_at_utc": utc_now(),
        "api_url": API_URL,
        "assembly": ASSEMBLY,
        "parameters": PARAMETERS,
        "input_count": len(batch_df),
        "variants": expected_inputs,
        "variant_keys": batch_df["variant_key"].tolist(),
    }

    write_json_atomic(paths["request"], request_data)

    # Reuse an existing response only after full identity validation.
    if paths["response"].exists():
        try:
            with open(paths["response"], encoding="utf-8") as handle:
                previous_results = json.load(handle)

            validate_response(batch_df, previous_results)

            logging.info(
                "Batch %04d has a valid response; skipping API call.",
                batch_number,
            )

            checkpoint["batches"][str(batch_number)] = {
                "status": "completed",
                "input_count": len(batch_df),
                "result_count": len(previous_results),
                "validated_at_utc": utc_now(),
            }

            save_checkpoint(checkpoint)
            return True

        except Exception as exc:
            logging.warning(
                "Existing batch %04d response invalid: %s",
                batch_number,
                exc,
            )
            archive_invalid_response(paths["response"])

    try:
        results = submit_to_vep(session, expected_inputs)

        # Preserve the unmodified API response before validation.
        write_json_atomic(paths["response"], results)

        validate_response(batch_df, results)

        checkpoint["batches"][str(batch_number)] = {
            "status": "completed",
            "input_count": len(batch_df),
            "result_count": len(results),
            "completed_at_utc": utc_now(),
        }

        save_checkpoint(checkpoint)

        logging.info(
            "Batch %04d validated: %s/%s variants.",
            batch_number,
            len(results),
            len(batch_df),
        )

        return True

    except Exception as exc:
        logging.exception("Batch %04d failed.", batch_number)

        if paths["response"].exists():
            archive_invalid_response(paths["response"])

        write_json_atomic(
            paths["error"],
            {
                "batch_number": batch_number,
                "failed_at_utc": utc_now(),
                "error": str(exc),
                "input_count": len(batch_df),
                "inputs": expected_inputs,
            },
        )

        checkpoint["batches"][str(batch_number)] = {
            "status": "failed",
            "input_count": len(batch_df),
            "error": str(exc),
            "failed_at_utc": utc_now(),
        }

        save_checkpoint(checkpoint)
        return False


# ============================================================
# TARGET-GENE MANE SELECT TRANSCRIPT EXTRACTION
# ============================================================

def extract_target_mane_record(source, result):
    consequences = result.get("transcript_consequences") or []

    target_records = [
        item
        for item in consequences
        if item.get("gene_symbol") == source["GeneSymbol"]
        and item.get("mane_select")
    ]

    if len(target_records) != 1:
        return None, (
            "target_gene_mane_select_count_"
            + str(len(target_records))
        )

    consequence = target_records[0]

    identity = source["vep_identity_key"]

    cadd_raw = safe_numeric_score(
        consequence.get("cadd_raw"),
        "cadd_raw",
        identity,
    )

    cadd_phred = safe_numeric_score(
        consequence.get("cadd_phred"),
        "cadd_phred",
        identity,
    )

    revel = safe_numeric_score(
        consequence.get("revel"),
        "revel",
        identity,
    )

    terms = consequence.get("consequence_terms") or []

    if isinstance(terms, list):
        consequence_terms = "|".join(str(term) for term in terms)
    else:
        consequence_terms = str(terms)

    record = {
        "#AlleleID": source["#AlleleID"],
        "GeneSymbol": source["GeneSymbol"],
        "variant_key": source["variant_key"],
        "vep_identity_key": identity,
        "Chromosome": source["Chromosome"],
        "PositionVCF": source["PositionVCF"],
        "ReferenceAlleleVCF": source["ReferenceAlleleVCF"],
        "AlternateAlleleVCF": source["AlternateAlleleVCF"],
        "transcript_id": consequence.get("transcript_id"),
        "mane_select": consequence.get("mane_select"),
        "gene_id": consequence.get("gene_id"),
        "biotype": consequence.get("biotype"),
        "consequence_terms": consequence_terms,
        "hgvsc": consequence.get("hgvsc"),
        "hgvsp": consequence.get("hgvsp"),
        "protein_id": consequence.get("protein_id"),
        "amino_acids": consequence.get("amino_acids"),
        "codons": consequence.get("codons"),
        "cadd_raw": cadd_raw,
        "cadd_phred": cadd_phred,
        "revel": revel,
        "cadd_raw_available": cadd_raw is not None,
        "cadd_phred_available": cadd_phred is not None,
        "revel_available": revel is not None,
        "mane_target_gene_validated": True,
        "assembly_name": result.get("assembly_name"),
        "vep_input": result.get("input"),
    }

    return record, None


# ============================================================
# FINAL DATASET ASSEMBLY AND QUALITY CONTROL
# ============================================================

def save_final_outputs(
    candidates,
    batch_list,
    successful_batches,
    failed_batches,
):
    annotation_records = []
    failure_records = []

    for batch_number, batch_df in batch_list:
        if batch_number not in successful_batches:
            for _, row in batch_df.iterrows():
                failure_records.append({
                    "#AlleleID": row["#AlleleID"],
                    "GeneSymbol": row["GeneSymbol"],
                    "variant_key": row["variant_key"],
                    "vep_identity_key": row["vep_identity_key"],
                    "batch_number": batch_number,
                    "failure_type": "batch_annotation_failed",
                })
            continue

        paths = batch_paths(batch_number)

        with open(paths["response"], encoding="utf-8") as handle:
            results = json.load(handle)

        validated = validate_response(batch_df, results)

        source_map = {
            row["vep_identity_key"]: row
            for _, row in batch_df.iterrows()
        }

        for identity, result in validated.items():
            source = source_map[identity]

            record, error = extract_target_mane_record(
                source,
                result,
            )

            if error:
                failure_records.append({
                    "#AlleleID": source["#AlleleID"],
                    "GeneSymbol": source["GeneSymbol"],
                    "variant_key": source["variant_key"],
                    "vep_identity_key": identity,
                    "batch_number": batch_number,
                    "failure_type": error,
                })
                continue

            annotation_records.append(record)

    annotation_df = pd.DataFrame(annotation_records)
    failure_df = pd.DataFrame(failure_records)

    failure_columns = [
        "#AlleleID",
        "GeneSymbol",
        "variant_key",
        "vep_identity_key",
        "batch_number",
        "failure_type",
    ]

    if failure_df.empty:
        failure_df = pd.DataFrame(columns=failure_columns)

    if not annotation_df.empty:
        annotation_df = annotation_df.merge(
            candidates[
                [
                    "#AlleleID",
                    "source_row_index",
                    "variant_key",
                ]
            ],
            on=["#AlleleID", "variant_key"],
            how="left",
            validate="one_to_one",
        )

        if annotation_df["source_row_index"].isna().any():
            raise RuntimeError(
                "Output contains records not mapped to source cohort."
            )

        annotation_df = (
            annotation_df
            .sort_values("source_row_index")
            .drop(columns=["source_row_index"])
            .reset_index(drop=True)
        )

    annotation_df.to_csv(OUTPUT_CSV, index=False)
    failure_df.to_csv(FAILURE_CSV, index=False)

    # --------------------------------------------------------
    # Final identity accounting
    # --------------------------------------------------------

    expected_keys = set(candidates["variant_key"].astype(str))

    returned_keys = (
        set(annotation_df["variant_key"].astype(str))
        if not annotation_df.empty
        else set()
    )

    missing_keys = expected_keys - returned_keys
    unexpected_keys = returned_keys - expected_keys

    duplicate_keys = (
        int(annotation_df["variant_key"].duplicated().sum())
        if not annotation_df.empty
        else 0
    )

    duplicate_allele_ids = (
        int(annotation_df["#AlleleID"].duplicated().sum())
        if not annotation_df.empty
        else 0
    )

    # --------------------------------------------------------
    # Predictor coverage and distributions
    # --------------------------------------------------------

    predictor_summary = {}

    for column in ["cadd_raw", "cadd_phred", "revel"]:
        if column in annotation_df.columns:
            values = pd.to_numeric(
                annotation_df[column],
                errors="coerce",
            )

            available = values.dropna()

            predictor_summary[column] = {
                "available_count": int(available.count()),
                "missing_count": int(values.isna().sum()),
                "coverage_fraction": (
                    float(available.count() / len(candidates))
                    if len(candidates)
                    else None
                ),
                "minimum": (
                    float(available.min())
                    if not available.empty
                    else None
                ),
                "median": (
                    float(available.median())
                    if not available.empty
                    else None
                ),
                "mean": (
                    float(available.mean())
                    if not available.empty
                    else None
                ),
                "maximum": (
                    float(available.max())
                    if not available.empty
                    else None
                ),
            }
        else:
            predictor_summary[column] = {
                "available_count": 0,
                "missing_count": int(len(candidates)),
                "coverage_fraction": 0.0,
                "minimum": None,
                "median": None,
                "mean": None,
                "maximum": None,
            }

    gene_coverage = {}

    for gene, group in candidates.groupby("GeneSymbol"):
        gene_keys = set(group["variant_key"].astype(str))

        matched = (
            annotation_df[
                annotation_df["variant_key"].astype(str).isin(gene_keys)
            ]
            if not annotation_df.empty
            else pd.DataFrame()
        )

        gene_coverage[str(gene)] = {
            "input_variants": int(len(group)),
            "selected_mane_records": int(len(matched)),
            "cadd_raw_available": (
                int(matched["cadd_raw"].notna().sum())
                if not matched.empty else 0
            ),
            "cadd_phred_available": (
                int(matched["cadd_phred"].notna().sum())
                if not matched.empty else 0
            ),
            "revel_available": (
                int(matched["revel"].notna().sum())
                if not matched.empty else 0
            ),
        }

    # --------------------------------------------------------
    # Report and integrity status
    # --------------------------------------------------------

    all_variants_have_mane = (
        len(annotation_df) == len(candidates)
        and len(missing_keys) == 0
        and len(unexpected_keys) == 0
        and duplicate_keys == 0
        and duplicate_allele_ids == 0
        and len(failure_df) == 0
    )

    report = {
        "source": "Ensembl Variant Effect Predictor REST API",
        "created_at_utc": utc_now(),
        "assembly": ASSEMBLY,
        "api_url": API_URL,
        "parameters": PARAMETERS,
        "batch_size": BATCH_SIZE,
        "input_variant_count": int(len(candidates)),
        "total_batches": int(len(batch_list)),
        "successful_batches": int(len(successful_batches)),
        "failed_batches": int(len(failed_batches)),
        "selected_target_gene_mane_count": int(len(annotation_df)),
        "failed_variant_count": int(len(failure_df)),
        "missing_variant_keys": int(len(missing_keys)),
        "unexpected_variant_keys": int(len(unexpected_keys)),
        "duplicate_variant_keys": duplicate_keys,
        "duplicate_allele_ids": duplicate_allele_ids,
        "all_variants_have_exactly_one_target_gene_mane": (
            bool(all_variants_have_mane)
        ),
        "predictor_summary": predictor_summary,
        "gene_coverage": gene_coverage,
        "failure_type_distribution": (
            {
                str(key): int(value)
                for key, value in (
                    failure_df["failure_type"].value_counts().items()
                )
            }
            if not failure_df.empty else {}
        ),
        "outputs": {
            "annotations": str(OUTPUT_CSV.relative_to(ROOT)),
            "failures": str(FAILURE_CSV.relative_to(ROOT)),
        },
        "sha256": {
            "input": sha256_file(INPUT_CSV),
            "annotations": (
                sha256_file(OUTPUT_CSV)
                if OUTPUT_CSV.exists() else None
            ),
            "failures": (
                sha256_file(FAILURE_CSV)
                if FAILURE_CSV.exists() else None
            ),
        },
    }

    write_json_atomic(REPORT_JSON, report)

    logging.info("=" * 65)
    logging.info("FULL VEP PREDICTOR ANNOTATION SUMMARY")
    logging.info("=" * 65)
    logging.info("Input variants: %s", len(candidates))
    logging.info("Successful API batches: %s", len(successful_batches))
    logging.info("Failed API batches: %s", len(failed_batches))
    logging.info("Selected target-gene MANE records: %s", len(annotation_df))
    logging.info("Failed variants: %s", len(failure_df))
    logging.info("Missing identities: %s", len(missing_keys))

    for predictor, summary in predictor_summary.items():
        logging.info(
            "%s coverage: %s/%s",
            predictor,
            summary["available_count"],
            len(candidates),
        )

    logging.info("Report: %s", REPORT_JSON)

    return bool(all_variants_have_mane)


# ============================================================
# MAIN
# ============================================================

def main():
    prepare_directories()

    candidates = load_candidates()

    configuration = build_run_configuration(candidates)
    initialize_manifest(configuration)

    checkpoint = load_checkpoint()

    batch_list = [
        (
            index + 1,
            candidates.iloc[start:start + BATCH_SIZE].copy(),
        )
        for index, start in enumerate(
            range(0, len(candidates), BATCH_SIZE)
        )
    ]

    logging.info("Total variants: %s", len(candidates))
    logging.info("Batch size: %s", BATCH_SIZE)
    logging.info("Total batches: %s", len(batch_list))

    session = requests.Session()

    successful_batches = set()
    failed_batches = set()

    try:
        for batch_number, batch_df in batch_list:
            logging.info("-" * 60)
            logging.info(
                "Processing batch %04d/%04d (%s variants)",
                batch_number,
                len(batch_list),
                len(batch_df),
            )

            success = process_batch(
                batch_number,
                batch_df,
                session,
                checkpoint,
            )

            if success:
                successful_batches.add(batch_number)
                time.sleep(PAUSE_BETWEEN_BATCHES)
            else:
                failed_batches.add(batch_number)
                time.sleep(5)

    finally:
        session.close()

    all_successful = save_final_outputs(
        candidates,
        batch_list,
        successful_batches,
        failed_batches,
    )

    manifest_path = RUN_DIR / "run_manifest.json"

    with open(manifest_path, encoding="utf-8") as handle:
        manifest = json.load(handle)

    manifest["finished_at_utc"] = utc_now()
    manifest["status"] = (
        "completed" if all_successful else "incomplete"
    )

    write_json_atomic(manifest_path, manifest)

    if not all_successful:
        raise SystemExit(
            "VEP predictor annotation is incomplete. "
            "Review the failure CSV and batch error files, "
            "then rerun this script to resume."
        )

    logging.info("=" * 65)
    logging.info(
        "FULL VEP PREDICTOR ANNOTATION COMPLETED SUCCESSFULLY"
    )
    logging.info("=" * 65)


if __name__ == "__main__":
    main()
