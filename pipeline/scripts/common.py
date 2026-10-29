from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path


def config(path: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        key, value = line.split(":", 1)
        out[key.strip()] = value.split("#", 1)[0].strip().strip("'\"")
    return out


def ensure_parent(path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def dump_json(path: str | Path, value: object) -> None:
    p = ensure_parent(path)
    p.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def open_text(path: str | Path, mode: str = "rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode, encoding="utf-8")


def fastq(path: str | Path):
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline().strip()
            handle.readline()
            qual = handle.readline().strip()
            if not header.startswith("@") or len(seq) != len(qual):
                raise ValueError(f"Malformed FASTQ record in {path}")
            yield header[1:].split()[0], seq, qual


def write_tsv(path: str | Path, columns: list[str], rows) -> None:
    p = ensure_parent(path)
    with gzip.open(p, "wt", encoding="utf-8", newline="") if str(p).endswith(".gz") else open(p, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
