from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--barcode", required=True); ap.add_argument("--structure", required=True)
    ap.add_argument("--output", required=True); ap.add_argument("--summary", required=True)
    a = ap.parse_args(); out = Path(a.output); out.parent.mkdir(parents=True, exist_ok=True)
    kept = total = 0
    with gzip.open(a.barcode, "rt", encoding="utf-8", errors="replace") as bsrc, gzip.open(a.structure, "rt", encoding="utf-8", errors="replace") as ssrc, gzip.open(out, "wt", encoding="utf-8", newline="") as dst:
        b = csv.DictReader(bsrc, delimiter="\t"); s = csv.DictReader(ssrc, delimiter="\t")
        writer = csv.DictWriter(dst, fieldnames=["read_id", "barcode", "umi", "x", "y"], delimiter="\t"); writer.writeheader()
        for brow, srow in zip(b, s):
            total += 1
            if brow.get("read_id") != srow.get("read_id"):
                raise SystemExit(f"Input order mismatch at read {total}: {brow.get('read_id')} != {srow.get('read_id')}")
            if brow.get("assign_status") not in {"exact", "rescued"} or brow.get("barcode_corr", ".") in {"", "."} or srow.get("umi_seq", ".") in {"", "."}:
                continue
            writer.writerow({"read_id": brow["read_id"], "barcode": brow["barcode_corr"], "umi": srow["umi_seq"], "x": brow.get("x", "-1"), "y": brow.get("y", "-1")})
            kept += 1
    Path(a.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(a.summary).write_text(json.dumps({"total_reads": total, "kept_barcode_umi_reads": kept, "source": "validated_reference_EXAMPLE_SAMPLE_barcode_and_structure_tables"}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
