from pathlib import Path
import json
import logging
import time
from datetime import datetime, timezone

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    ROOT / "data/processed/"
    "brca1_brca2_primary_missense_cohort.csv"
)

OUTPUT_DIR = ROOT / "data/raw/vep/predictor_pilot"

ENDPOINT = (
    "https://rest.ensembl.org/vep/human/region"
)

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

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    cohort = pd.read_csv(INPUT_FILE, low_memory=False)

    # Ten variants per gene, reproducible selection.
    sample = (
        cohort.groupby("GeneSymbol", group_keys=False)
        .sample(n=10, random_state=42)
        .sort_values(["GeneSymbol", "variant_key"])
    )

    variants = []

    for key in sample["variant_key"]:
        # Expected: GRCh38:chromosome:position:REF>ALT
        build, chrom, position, substitution = key.split(":")

        if build != "GRCh38":
            raise ValueError(f"Unexpected genome build: {build}")

        ref, alt = substitution.split(">")

        variants.append(
            f"{chrom} {position} {position} {ref}/{alt} 1"
        )

    request_payload = {"variants": variants}

    request_record = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "endpoint": ENDPOINT,
        "parameters": PARAMETERS,
        "variants": variants,
        "variant_keys": sample["variant_key"].tolist(),
    }

    with open(
        OUTPUT_DIR / "request.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(request_record, handle, indent=4)

    logging.info("Submitting %s pilot variants.", len(variants))
    logging.info("Predictor parameters: %s", PARAMETERS)

    response = None

    for attempt in range(1, 6):
        try:
            response = requests.post(
                ENDPOINT,
                params=PARAMETERS,
                headers={
                    "Content-Type": "application/json",
                    "Accept": "application/json",
                    "User-Agent": "BRCA-functional-variant-research",
                },
                json=request_payload,
                timeout=180,
            )

            if response.status_code == 429 or response.status_code >= 500:
                wait = min(15 * attempt, 90)
                logging.warning(
                    "HTTP %s. Retrying in %s seconds.",
                    response.status_code,
                    wait,
                )
                time.sleep(wait)
                continue

            response.raise_for_status()
            break

        except requests.RequestException:
            if attempt == 5:
                raise

            logging.exception("VEP request failed on attempt %s.", attempt)
            time.sleep(10 * attempt)

    if response is None or not response.ok:
        raise RuntimeError("VEP predictor pilot failed.")

    results = response.json()

    if not isinstance(results, list):
        raise TypeError("Unexpected VEP response format.")

    with open(
        OUTPUT_DIR / "response.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(results, handle, indent=2)

    if len(results) != len(variants):
        raise RuntimeError(
            f"Expected {len(variants)} results, received {len(results)}."
        )

    print("\nVEP PREDICTOR PILOT")
    print("Submitted variants:", len(variants))
    print("Returned variants:", len(results))
    print("HTTP status:", response.status_code)

    print("\nTOP-LEVEL RESPONSE FIELDS")

    for result in results[:3]:
        print("\nVariant:", result.get("input"))
        print("Fields:", sorted(result.keys()))

        consequences = result.get("transcript_consequences", [])

        for consequence in consequences[:2]:
            print(
                "Transcript:",
                consequence.get("transcript_id"),
            )
            print(
                "Consequence fields:",
                sorted(consequence.keys()),
            )

    print("\nPREDICTOR FIELD INSPECTION")

    for predictor in ["CADD", "REVEL"]:
        variant_hits = 0
        transcript_hits = 0

        for result in results:
            if predictor in result:
                variant_hits += 1

            for consequence in result.get(
                "transcript_consequences", []
            ):
                if predictor in consequence:
                    transcript_hits += 1

        print(
            f"{predictor}: "
            f"variant-level fields={variant_hits}; "
            f"transcript-level fields={transcript_hits}"
        )

    print("\nPilot request and response saved to:")
    print(OUTPUT_DIR)
    print("\nVEP PREDICTOR PILOT COMPLETED.")

if __name__ == "__main__":
    main()

