import hashlib
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]

OUTPUT_DIR = ROOT / "data/raw/alphamissense"
OUTPUT_FILE = OUTPUT_DIR / "AlphaMissense_hg38.tsv.gz"

MANIFEST_FILE = (
    ROOT / "data/metadata/"
    "alphamissense_download_manifest.json"
)

URL = (
    "https://zenodo.org/records/8208688/files/"
    "AlphaMissense_hg38.tsv.gz?download=1"
)

CHUNK_SIZE = 1024 * 1024

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

logger = logging.getLogger(__name__)


def calculate_sha256(path):

    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(CHUNK_SIZE), b""
        ):
            digest.update(chunk)

    return digest.hexdigest()


def main():

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if OUTPUT_FILE.exists():
        raise FileExistsError(
            f"File already exists: {OUTPUT_FILE}. "
            "Do not overwrite an existing download."
        )

    temporary_file = OUTPUT_FILE.with_suffix(
        OUTPUT_FILE.suffix + ".part"
    )

    logger.info("Downloading official AlphaMissense GRCh38 data.")

    try:

        digest = hashlib.sha256()
        downloaded = 0

        with requests.get(
            URL,
            stream=True,
            timeout=(30, 180),
        ) as response:

            response.raise_for_status()

            total = int(
                response.headers.get("Content-Length", 0)
            )

            with temporary_file.open("wb") as handle:

                for chunk in response.iter_content(
                    chunk_size=CHUNK_SIZE
                ):

                    if not chunk:
                        continue

                    handle.write(chunk)
                    digest.update(chunk)
                    downloaded += len(chunk)

                    if downloaded % (100 * CHUNK_SIZE) < CHUNK_SIZE:
                        logger.info(
                            "Downloaded %.1f MB",
                            downloaded / (1024 * 1024),
                        )

        if total and downloaded != total:
            raise RuntimeError(
                f"Download size mismatch: "
                f"{downloaded} != {total}"
            )

        # Validate gzip integrity before accepting the file.
        import gzip

        with gzip.open(temporary_file, "rb") as handle:
            header = handle.read(1024)

        if b"#CHROM" not in header:
            raise ValueError(
                "Expected AlphaMissense header not found."
            )

        temporary_file.replace(OUTPUT_FILE)

        manifest = {
            "source": "AlphaMissense",
            "source_url": URL,
            "filename": OUTPUT_FILE.name,
            "genome_assembly": "GRCh38.p13",
            "transcript_annotation": "GENCODE V32 canonical transcripts",
            "retrieved_at_utc": datetime.now(
                timezone.utc
            ).isoformat(),
            "file_size_bytes": downloaded,
            "sha256": digest.hexdigest(),
            "gzip_integrity": "passed",
            "header_validation": "passed",
            "license": "CC BY 4.0",
            "citation": (
                "Cheng et al. Accurate proteome-wide missense "
                "variant effect prediction with AlphaMissense. "
                "Science (2023)."
            ),
        }

        MANIFEST_FILE.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        MANIFEST_FILE.write_text(
            json.dumps(manifest, indent=2),
            encoding="utf-8",
        )

        logger.info("DOWNLOAD AND VALIDATION COMPLETED")
        logger.info("File: %s", OUTPUT_FILE)
        logger.info("Size: %.2f MB", downloaded / (1024 * 1024))
        logger.info("SHA256: %s", digest.hexdigest())
        logger.info("Manifest: %s", MANIFEST_FILE)

    except Exception:

        if temporary_file.exists():
            temporary_file.unlink()

        raise


if __name__ == "__main__":
    main()
