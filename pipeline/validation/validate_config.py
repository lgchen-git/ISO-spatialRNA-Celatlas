from __future__ import annotations
import argparse
import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from common import config


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = config(args.config)
    required = ("raw_fastq", "whitelist_csv", "bbox_csv", "barcode_to_pos_h5", "genome_fa", "gtf")
    missing = [key for key in required if not Path(cfg.get(key, "")).is_file()]
    sample = Path(cfg.get("sample_sheet", ""))
    if not sample.is_file():
        missing.append("sample_sheet")
    if missing:
        print("missing: " + ", ".join(missing))
        return 1
    with sample.open(encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="	"), {})
    if not row.get("sample_id"):
        print("sample sheet has no sample_id")
        return 1
    print("config: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
