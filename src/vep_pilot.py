import json
import logging
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]

INPUT_CSV = ROOT / "data/processed/brca1_brca2_grch38_candidate_snvs.csv"

RAW_DIR = ROOT / "data/raw/vep"
PROCESSED_DIR = ROOT / "data/processed"
METADATA_DIR = ROOT / "data/metadata"

API_URL = "https://rest.ensembl.org/vep/homo_sapiens/region"

PILOT_SIZE_PER_GENE = 10
RANDOM_STATE = 42

PARAMETERS = {
    "mane": 1,
    "hgvs": 1,
    "canonical": 1,
    "protein": 1,
    "numbers": 1,
    "tsl": 1,
}

REQUIRED_COLUMNS = [
    "GeneSymbol",
    "Chromosome",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)


def prepare_directories():
    for directory in [RAW_DIR, PROCESSED_DIR, METADATA_DIR]:
        directory.mkdir(parents=True, exist_ok=True)


def load_candidates():
    logging.info("Loading ClinVar candidate variants...")

    df = pd.read_csv(INPUT_CSV, low_memory=False)

    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    logging.info("Loaded %s candidate variants.", len(df))

    return df


def select_pilot(df):
    logging.info("Selecting a deterministic pilot sample...")

    selected = []

    for gene in ["BRCA1", "BRCA2"]:
        gene_df = df[df["GeneSymbol"] == gene]

        if len(gene_df) < PILOT_SIZE_PER_GENE:
            raise ValueError(f"Not enough variants for {gene}")

        selected.append(
            gene_df.sample(
                n=PILOT_SIZE_PER_GENE,
                random_state=RANDOM_STATE,
            )
        )

    pilot = pd.concat(selected).reset_index(drop=True)

    logging.info(
        "Selected %s variants: %s",
        len(pilot),
        pilot["GeneSymbol"].value_counts().to_dict(),
    )

    return pilot


def make_vep_input(row):
    chromosome = str(row["Chromosome"]).strip()

    if chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]

    position = int(row["PositionVCF"])
    ref = str(row["ReferenceAlleleVCF"]).upper()
    alt = str(row["AlternateAlleleVCF"]).upper()

    if not ref or not alt or ref == "N" or alt == "N":
        raise ValueError(f"Invalid allele: {ref}>{alt}")

    # VEP region input: chromosome start ID ref alt QUAL FILTER INFO
    return f"{chromosome} {position} . {ref} {alt} . . ."


def submit_to_vep(variants):
    payload = {"variants": variants}

    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "BRCA-Functional-Variant-Analysis/1.0",
    }

    for attempt in range(1, 5):
        logging.info("Submitting VEP request (attempt %s/4)...", attempt)

        try:
            response = requests.post(
                API_URL,
                params=PARAMETERS,
                headers=headers,
                json=payload,
                timeout=120,
            )

            if response.status_code == 429 or response.status_code >= 500:
                if attempt == 4:
                    response.raise_for_status()

                retry_after = response.headers.get("Retry-After")

                try:
                    wait_seconds = int(retry_after) if retry_after else 2 ** attempt
                except ValueError:
                    wait_seconds = 2 ** attempt

                logging.warning(
                    "VEP returned HTTP %s. Retrying in %s seconds.",
                    response.status_code,
                    wait_seconds,
                )

                time.sleep(wait_seconds)
                continue

            response.raise_for_status()

            return response.json()

        except requests.RequestException:
            if attempt == 4:
                raise

            time.sleep(2 ** attempt)

    raise RuntimeError("VEP request failed after all retry attempts.")


def normalize_input(value):
    return " ".join(str(value).split())


def process_response(pilot, results):
    input_lookup = {
        normalize_input(row["vep_input"]): row.to_dict()
        for _, row in pilot.iterrows()
    }

    variant_records = []
    transcript_records = []

    for result in results:
        result_input = normalize_input(result.get("input", ""))

        source = input_lookup.get(result_input)

        if source is None:
            logging.warning("Could not match VEP result: %s", result_input)
            source = {}

        transcript_consequences = result.get(
            "transcript_consequences", []
        )

        variant_record = {
            **source,
            "vep_input": result_input,
            "assembly_name": result.get("assembly_name"),
            "seq_region_name": result.get("seq_region_name"),
            "start": result.get("start"),
            "end": result.get("end"),
            "allele_string": result.get("allele_string"),
            "most_severe_consequence": result.get(
                "most_severe_consequence"
            ),
            "transcript_consequence_count": len(transcript_consequences),
        }

        variant_records.append(variant_record)

        for consequence in transcript_consequences:
            transcript_records.append({
                **source,
                "vep_input": result_input,
                "assembly_name": result.get("assembly_name"),
                "most_severe_consequence": result.get(
                    "most_severe_consequence"
                ),
                **consequence,
                "consequence_terms": "|".join(
                    consequence.get("consequence_terms", [])
                ),
            })

    return (
        pd.DataFrame(variant_records),
        pd.DataFrame(transcript_records),
    )


def main():
    prepare_directories()

    df = load_candidates()
    pilot = select_pilot(df)

    pilot["vep_input"] = pilot.apply(make_vep_input, axis=1)

    request_payload = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "api_url": API_URL,
        "parameters": PARAMETERS,
        "variants": pilot["vep_input"].tolist(),
    }

    request_path = RAW_DIR / "vep_pilot_request.json"

    with open(request_path, "w", encoding="utf-8") as handle:
        json.dump(request_payload, handle, indent=2)

    results = submit_to_vep(pilot["vep_input"].tolist())

    response_path = RAW_DIR / "vep_pilot_response.json"

    with open(response_path, "w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)

    variant_df, transcript_df = process_response(pilot, results)

    variant_path = PROCESSED_DIR / "vep_pilot_variants.csv"
    transcript_path = (
        PROCESSED_DIR / "vep_pilot_transcript_consequences.csv"
    )

    variant_df.to_csv(variant_path, index=False)
    transcript_df.to_csv(transcript_path, index=False)

    report = {
        "source": "Ensembl Variant Effect Predictor REST API",
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "assembly": "GRCh38",
        "input_candidate_count": len(pilot),
        "vep_result_count": len(results),
        "transcript_consequence_rows": len(transcript_df),
        "gene_distribution": pilot["GeneSymbol"].value_counts().to_dict(),
        "consequence_distribution": (
            variant_df["most_severe_consequence"]
            .fillna("not_returned")
            .value_counts()
            .to_dict()
        ),
        "parameters": PARAMETERS,
        "api_url": API_URL,
    }

    report_path = METADATA_DIR / "vep_pilot_report.json"

    with open(report_path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)

    logging.info("=" * 55)
    logging.info("VEP PILOT COMPLETED")
    logging.info("=" * 55)
    logging.info("Input variants: %s", len(pilot))
    logging.info("VEP results: %s", len(results))
    logging.info("Transcript consequence rows: %s", len(transcript_df))
    logging.info("Variant output: %s", variant_path)
    logging.info("Transcript output: %s", transcript_path)
    logging.info("Report: %s", report_path)


if __name__ == "__main__":
    main()