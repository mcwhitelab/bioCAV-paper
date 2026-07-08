#!/usr/bin/env python3

import argparse
from pathlib import Path
from urllib.parse import quote

import requests


def read_accessions(path: Path) -> list[str]:
    return [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


def batched(items: list[str], batch_size: int):
    for i in range(0, len(items), batch_size):
        yield items[i:i + batch_size]


def fetch_fasta(accessions: list[str], timeout: int) -> str:
    query = " OR ".join(f"accession:{acc}" for acc in accessions)
    url = (
        "https://rest.uniprot.org/uniprotkb/stream"
        f"?query={quote(query)}"
        "&format=fasta"
    )

    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    return response.text


def main():
    parser = argparse.ArgumentParser(
        description="Download FASTA sequences from UniProt for accessions in a text file."
    )
    parser.add_argument(
        "ids_file",
        type=Path,
        help="Text file containing one UniProt accession per line.",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("uniprot_sequences.fasta"),
        help="Output FASTA file. Default: uniprot_sequences.fasta",
    )
    parser.add_argument(
        "-b", "--batch-size",
        type=int,
        default=100,
        help="Number of accessions per UniProt request. Default: 100",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="Request timeout in seconds. Default: 60",
    )

    args = parser.parse_args()

    accessions = read_accessions(args.ids_file)

    if not accessions:
        raise SystemExit(f"No accessions found in {args.ids_file}")

    with args.output.open("w") as out:
        for i, batch in enumerate(batched(accessions, args.batch_size), start=1):
            print(f"Downloading batch {i} ({len(batch)} accessions)...")
            fasta = fetch_fasta(batch, args.timeout)
            out.write(fasta)
            if not fasta.endswith("\n"):
                out.write("\n")

    print(f"Done. Wrote {args.output}")


if __name__ == "__main__":
    main()
