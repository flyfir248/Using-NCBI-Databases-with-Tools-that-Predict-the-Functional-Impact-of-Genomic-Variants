
import hashlib
import json
import logging
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
    "brca1_brca2_grch38_candidate_snvs.csv"
)

RUN_DIR = ROOT / "data/raw/vep/full_run"
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
}

REQUIRED_COLUMNS = [
    "#AlleleID",
    "GeneSymbol",
    "Chromosome",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
]

RETRY_STATUS_CODES = {408, 429, 500, 502, 503, 504}

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
        json.dump(data, handle, indent=2)

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
        raise ValueError(f"Invalid position: {position}")

    if ref not in {"A", "C", "G", "T"}:
        raise ValueError(f"Invalid reference allele: {ref}")

    if alt not in {"A", "C", "G", "T"}:
        raise ValueError(f"Invalid alternate allele: {alt}")

    if ref == alt:
        raise ValueError(f"Reference and alternate alleles are identical: {ref}")

    # VEP region input:
    # chromosome start ID ref alt QUAL FILTER INFO
    return f"{chromosome} {position} . {ref} {alt} . . ."


def identity_key_from_vep_input(value):
    fields = str(value).split()

    if len(fields) < 5:
        raise ValueError(f"Unexpected VEP input representation: {value}")

    chromosome = fields[0]
    position = int(fields[1])
    ref = fields[3]
    alt = fields[4]

    return make_identity_key(chromosome, position, ref, alt)


def identity_key_from_result(result):
    if not result.get("input"):
        raise ValueError("VEP result has no input field.")

    return identity_key_from_vep_input(result["input"])


def prepare_directories():
    for directory in [
        RUN_DIR,
        BATCH_DIR,
        PROCESSED_DIR,
        METADATA_DIR,
    ]:
        directory.mkdir(parents=True, exist_ok=True)


# ============================================================
# INPUT DATA AND RUN CONFIGURATION
# ============================================================

def load_candidates():
    logging.info("Loading ClinVar candidates...")

    df = pd.read_csv(INPUT_CSV, low_memory=False)

    missing = [
        column
        for column in REQUIRED_COLUMNS
        if column not in df.columns
    ]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if df.empty:
        raise ValueError("Candidate dataset is empty.")

    if df[REQUIRED_COLUMNS].isna().any().any():
        raise ValueError(
            "Required source fields contain missing values. "
            "Resolve these before annotation."
        )

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

    if df["vep_identity_key"].duplicated().any():
        duplicates = df.loc[
            df["vep_identity_key"].duplicated(keep=False),
            ["#AlleleID", "GeneSymbol", "vep_identity_key"],
        ]

        raise ValueError(
            "Duplicate genomic identities detected:\n"
            + duplicates.head(20).to_string(index=False)
        )

    if df["#AlleleID"].duplicated().any():
        raise ValueError("Duplicate ClinVar AlleleIDs detected.")

    logging.info("Loaded %s unique candidates.", len(df))

    return df


def build_run_configuration(df):
    return {
        "source_file": str(INPUT_CSV.relative_to(ROOT)),
        "source_sha256": sha256_file(INPUT_CSV),
        "input_candidate_count": int(len(df)),
        "assembly": ASSEMBLY,
        "api_url": API_URL,
        "parameters": PARAMETERS,
        "batch_size": BATCH_SIZE,
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
                    "Do not resume this run with changed inputs or settings."
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
            return json.load(handle)

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

        # Non-retryable HTTP errors stop this batch.
        response.raise_for_status()

        try:
            results = response.json()
        except ValueError as exc:
            raise RuntimeError(
                "VEP returned a response that is not valid JSON."
            ) from exc

        if not isinstance(results, list):
            raise RuntimeError(
                "Unexpected VEP response structure: expected a JSON list."
            )

        return results

    raise RuntimeError("VEP request failed after all retry attempts.")


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
            raise ValueError("VEP response contains a non-object result.")

        key = identity_key_from_result(result)

        if key in returned_map:
            raise ValueError(f"Duplicate VEP result: {key}")

        if key not in expected_map:
            raise ValueError(f"Unexpected VEP result: {key}")

        returned_map[key] = result

    missing = set(expected_map) - set(returned_map)
    unexpected = set(returned_map) - set(expected_map)

    if missing:
        raise ValueError(
            f"{len(missing)} submitted variants are missing from VEP response. "
            f"Examples: {list(missing)[:5]}"
        )

    if unexpected:
        raise ValueError(
            f"Unexpected results: {list(unexpected)[:5]}"
        )

    if len(results) != len(expected_rows):
        raise ValueError(
            f"Expected {len(expected_rows)} results, "
            f"received {len(results)}."
        )

    for key, result in returned_map.items():
        source = expected_map[key]

        if result.get("assembly_name") != ASSEMBLY:
            raise ValueError(
                f"Assembly mismatch for {key}: "
                f"{result.get('assembly_name')}"
            )

        returned_chromosome = normalize_chromosome(
            result.get("seq_region_name")
        )

        expected_chromosome = normalize_chromosome(
            source["Chromosome"]
        )

        if returned_chromosome != expected_chromosome:
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
                f"expected {expected_alleles}, "
                f"received {returned_alleles}"
            )

    return returned_map


# ============================================================
# BATCH CHECKPOINT AND RAW DATA MANAGEMENT
# ============================================================

def batch_paths(batch_number):
    batch_name = f"batch_{batch_number:04d}"
    directory = BATCH_DIR / batch_name

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

    # Save request provenance before making the API call.
    request_data = {
        "batch_number": batch_number,
        "created_at_utc": utc_now(),
        "api_url": API_URL,
        "assembly": ASSEMBLY,
        "parameters": PARAMETERS,
        "input_count": len(batch_df),
        "variants": expected_inputs,
    }

    write_json_atomic(paths["request"], request_data)

    # A previous response can be reused only after full validation.
    if paths["response"].exists():
        try:
            with open(paths["response"], encoding="utf-8") as handle:
                previous_results = json.load(handle)

            validate_response(batch_df, previous_results)

            logging.info(
                "Batch %04d already has a valid response. Skipping API call.",
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
                "Existing batch %04d response failed validation: %s",
                batch_number,
                exc,
            )

            archive_invalid_response(paths["response"])

    try:
        results = submit_to_vep(session, expected_inputs)

        # Preserve the raw response before validation.
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
            "Batch %04d completed and validated: %s/%s",
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
# RESULT FLATTENING
# ============================================================

def flatten_batch(batch_df, results):
    validated = validate_response(batch_df, results)

    variant_records = []
    transcript_records = []

    for key, result in validated.items():
        source = batch_df.loc[
            batch_df["vep_identity_key"] == key
        ].iloc[0].to_dict()

        consequences = result.get("transcript_consequences") or []

        variant_record = {
            **source,
            "vep_input": result.get("input"),
            "assembly_name": result.get("assembly_name"),
            "seq_region_name": result.get("seq_region_name"),
            "start": result.get("start"),
            "end": result.get("end"),
            "allele_string": result.get("allele_string"),
            "most_severe_consequence": result.get(
                "most_severe_consequence"
            ),
            "transcript_consequence_count": len(consequences),
        }

        variant_records.append(variant_record)

        for consequence in consequences:
            transcript_record = {
                **source,
                "vep_input": result.get("input"),
                "assembly_name": result.get("assembly_name"),
                "most_severe_consequence": result.get(
                    "most_severe_consequence"
                ),
                **consequence,
            }

            terms = consequence.get("consequence_terms", [])

            transcript_record["consequence_terms"] = "|".join(
                terms if isinstance(terms, list) else [str(terms)]
            )

            transcript_records.append(transcript_record)

    return variant_records, transcript_records


# ============================================================
# FINAL DATASET ASSEMBLY AND QC
# ============================================================

def save_final_outputs(
    candidates,
    batch_list,
    successful_batches,
    failed_batches,
):
    variant_records = []
    transcript_records = []
    failure_records = []

    for batch_number, batch_df in batch_list:
        paths = batch_paths(batch_number)

        if batch_number in successful_batches:
            with open(paths["response"], encoding="utf-8") as handle:
                results = json.load(handle)

            variants, transcripts = flatten_batch(batch_df, results)

            variant_records.extend(variants)
            transcript_records.extend(transcripts)

        else:
            for _, row in batch_df.iterrows():
                failure_records.append({
                    "#AlleleID": row["#AlleleID"],
                    "GeneSymbol": row["GeneSymbol"],
                    "vep_identity_key": row["vep_identity_key"],
                    "vep_input": row["vep_input"],
                    "batch_number": batch_number,
                    "failure_type": "batch_annotation_failed",
                })

    variant_df = pd.DataFrame(variant_records)
    transcript_df = pd.DataFrame(transcript_records)
    failure_df = pd.DataFrame(failure_records)

    # Preserve a predictable schema even when no failures occur.
    if failure_df.empty:
        failure_df = pd.DataFrame(
            columns=[
                "#AlleleID",
                "GeneSymbol",
                "vep_identity_key",
                "vep_input",
                "batch_number",
                "failure_type",
            ]
        )

    if not variant_df.empty:
        variant_df = variant_df.sort_values(
            "source_row_index"
        ).reset_index(drop=True)

    if not transcript_df.empty:
        transcript_df = transcript_df.sort_values(
            "source_row_index"
        ).reset_index(drop=True)

    variant_path = PROCESSED_DIR / "vep_full_variants.csv"
    transcript_path = (
        PROCESSED_DIR / "vep_full_transcript_consequences.csv"
    )
    failure_path = PROCESSED_DIR / "vep_full_failures.csv"

    variant_df.to_csv(variant_path, index=False)
    transcript_df.to_csv(transcript_path, index=False)
    failure_df.to_csv(failure_path, index=False)

    # --------------------------------------------------------
    # MANE and canonical transcript coverage
    # --------------------------------------------------------

    if not transcript_df.empty:
        transcript_df["is_mane_select"] = (
            transcript_df["mane"]
            .fillna("")
            .astype(str)
            .str.contains("MANE_Select")
        )

        transcript_df["is_canonical"] = (
            pd.to_numeric(
                transcript_df["canonical"],
                errors="coerce",
            ) == 1
        )

        coverage = (
            transcript_df.groupby("#AlleleID")
            .agg(
                transcript_count=("transcript_id", "nunique"),
                mane_select_count=("is_mane_select", "sum"),
                canonical_count=("is_canonical", "sum"),
            )
            .reset_index()
        )

    else:
        coverage = pd.DataFrame(
            columns=[
                "#AlleleID",
                "transcript_count",
                "mane_select_count",
                "canonical_count",
            ]
        )

    coverage = candidates[
        ["#AlleleID", "GeneSymbol"]
    ].merge(
        coverage,
        on="#AlleleID",
        how="left",
    )

    for column in [
        "transcript_count",
        "mane_select_count",
        "canonical_count",
    ]:
        coverage[column] = coverage[column].fillna(0).astype(int)

    coverage["has_mane_select"] = (
        coverage["mane_select_count"] > 0
    )

    coverage["has_canonical"] = (
        coverage["canonical_count"] > 0
    )

    coverage_path = (
        METADATA_DIR / "vep_full_transcript_coverage.csv"
    )

    coverage.to_csv(coverage_path, index=False)

    # --------------------------------------------------------
    # Final integrity checks
    # --------------------------------------------------------

    expected_ids = set(candidates["#AlleleID"].astype(str))

    returned_ids = (
        set(variant_df["#AlleleID"].astype(str))
        if not variant_df.empty
        else set()
    )

    missing_ids = expected_ids - returned_ids
    unexpected_ids = returned_ids - expected_ids

    duplicate_ids = (
        int(variant_df["#AlleleID"].duplicated().sum())
        if not variant_df.empty
        else 0
    )

    if len(returned_ids) != len(variant_df):
        logging.error("Duplicate ClinVar IDs in final variant output.")

    if unexpected_ids:
        logging.error("Unexpected ClinVar IDs in final output.")

    if len(variant_df) + len(failure_df) != len(candidates):
        logging.error(
            "Final accounting does not equal the input candidate count."
        )

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------

    consequence_distribution = {}

    if not variant_df.empty:
        consequence_distribution = {
            str(key): int(value)
            for key, value in (
                variant_df["most_severe_consequence"]
                .fillna("not_returned")
                .value_counts()
                .items()
            )
        }

    gene_distribution = {
        str(key): int(value)
        for key, value in (
            candidates["GeneSymbol"].value_counts().items()
        )
    }

    report = {
        "source": "Ensembl Variant Effect Predictor REST API",
        "completed_at_utc": utc_now(),
        "assembly": ASSEMBLY,
        "api_url": API_URL,
        "parameters": PARAMETERS,
        "batch_size": BATCH_SIZE,
        "input_candidate_count": int(len(candidates)),
        "total_batches": len(batch_list),
        "successful_batches": len(successful_batches),
        "failed_batches": len(failed_batches),
        "returned_variant_count": int(len(variant_df)),
        "transcript_consequence_rows": int(len(transcript_df)),
        "failed_variant_count": int(len(failure_df)),
        "missing_clinvar_ids": len(missing_ids),
        "unexpected_clinvar_ids": len(unexpected_ids),
        "duplicate_clinvar_ids": duplicate_ids,
        "gene_distribution": gene_distribution,
        "consequence_distribution": consequence_distribution,
        "transcript_coverage": {
            "variants_with_mane_select": int(
                coverage["has_mane_select"].sum()
            ),
            "variants_without_mane_select": int(
                (~coverage["has_mane_select"]).sum()
            ),
            "variants_with_canonical": int(
                coverage["has_canonical"].sum()
            ),
            "variants_without_canonical": int(
                (~coverage["has_canonical"]).sum()
            ),
        },
        "outputs": {
            "variants": str(variant_path.relative_to(ROOT)),
            "transcript_consequences": str(
                transcript_path.relative_to(ROOT)
            ),
            "failures": str(failure_path.relative_to(ROOT)),
            "transcript_coverage": str(coverage_path.relative_to(ROOT)),
        },
    }

    report_path = METADATA_DIR / "vep_full_report.json"
    write_json_atomic(report_path, report)

    logging.info("=" * 65)
    logging.info("FULL VEP ANNOTATION SUMMARY")
    logging.info("=" * 65)
    logging.info("Input candidates: %s", len(candidates))
    logging.info("Successful batches: %s", len(successful_batches))
    logging.info("Failed batches: %s", len(failed_batches))
    logging.info("Returned variants: %s", len(variant_df))
    logging.info("Transcript rows: %s", len(transcript_df))
    logging.info("Failed variants: %s", len(failure_df))
    logging.info(
        "Variants with MANE Select: %s",
        coverage["has_mane_select"].sum(),
    )
    logging.info(
        "Variants with canonical transcript: %s",
        coverage["has_canonical"].sum(),
    )
    logging.info("Report: %s", report_path)

    return len(failed_batches) == 0 and (
        len(variant_df) == len(candidates)
        and len(missing_ids) == 0
        and len(unexpected_ids) == 0
        and duplicate_ids == 0
        and len(failure_df) == 0
    )


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
        (index + 1, candidates.iloc[start:start + BATCH_SIZE].copy())
        for index, start in enumerate(
            range(0, len(candidates), BATCH_SIZE)
        )
    ]

    logging.info("Total candidates: %s", len(candidates))
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
            "VEP annotation is incomplete. "
            "Review the failure CSV and batch error files, "
            "then rerun this script to resume."
        )

    logging.info("=" * 65)
    logging.info("FULL VEP ANNOTATION COMPLETED SUCCESSFULLY")
    logging.info("=" * 65)


if __name__ == "__main__":
    main()
