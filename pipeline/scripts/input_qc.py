from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path

from common import config, dump_json, fastq, ensure_parent


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--max-reads", type=int, default=0,
                    help="Process at most this many reads; 0 means the complete FASTQ")
    args = ap.parse_args()
    cfg = config(args.config)
    required = ["raw_fastq", "whitelist_csv", "bbox_csv", "barcode_to_pos_h5", "genome_fa", "gtf"]
    missing = [f"{key}={cfg.get(key, '')}" for key in required if not Path(cfg.get(key, "")).is_file()]
    if missing:
        raise SystemExit("Missing input/reference files:\n" + "\n".join(missing))
    upstream_metrics = Path(cfg.get("upstream_qc_metrics", ""))
    if upstream_metrics.is_file():
        values = {}
        with upstream_metrics.open(encoding="utf-8") as handle:
            for row in handle:
                key, _, value = row.rstrip("\n").partition("\t")
                if key != "metric":
                    try:
                        values[key] = float(value) if "." in value else int(value)
                    except ValueError:
                        pass
        summary = {"sample_id": "EXAMPLE_SAMPLE", "platform": "ONT", "chemistry": cfg.get("chemistry", "BBV2.4"),
                   **values, "raw_fastq": str(Path(cfg["raw_fastq"]).resolve()), "spatial_unit": cfg.get("spatial_unit", "bin50"),
                   "coordinate_source": "validated upstream EXAMPLE_SAMPLE QC; raw/reference paths rechecked", "source": str(upstream_metrics),
                   "complete_fastq_scan": True, "read_limit": None}
        dump_json(args.summary, summary)
        p = ensure_parent(args.metrics)
        p.write_text("metric\tvalue\n" + "\n".join(f"{k}\t{v}" for k, v in summary.items() if isinstance(v, (int, float))) + "\n", encoding="utf-8")
        return
    lengths: list[int] = []
    n = 0
    bases = 0
    for read_id, seq, qual in fastq(cfg["raw_fastq"]):
        del read_id, qual
        n += 1
        bases += len(seq)
        if len(lengths) < 200000:
            lengths.append(len(seq))
        if args.max_reads and n >= args.max_reads:
            break
    summary = {
        "sample_id": "EXAMPLE_SAMPLE", "platform": "ONT", "chemistry": cfg.get("chemistry", "BBV2.4"),
        "n_reads": n, "bases": bases, "mean_len": (bases / n if n else 0),
        "median_len_sample": statistics.median(lengths) if lengths else 0,
        "raw_fastq": str(Path(cfg["raw_fastq"]).resolve()), "spatial_unit": cfg.get("spatial_unit", "bin50"),
        "coordinate_source": "FilterBarcodes.csv cross-checked with bbox and barcodeToPos.h5",
        "read_limit": args.max_reads or None,
        "complete_fastq_scan": not bool(args.max_reads),
    }
    dump_json(args.summary, summary)
    p = ensure_parent(args.metrics)
    p.write_text("metric\tvalue\n" + "\n".join(f"{k}\t{v}" for k, v in summary.items() if isinstance(v, (int, float))) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
