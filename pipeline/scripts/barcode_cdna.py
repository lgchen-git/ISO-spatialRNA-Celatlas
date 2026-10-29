from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

from common import config, dump_json, fastq, ensure_parent

ANCHOR = "CTACACGACGCTCTTCCGATCT"
LINKER1 = ("CGACTCACTACACGT", "TCGACTCACTACACGT", "GTCGACTCACTACACGT", "AGTCGACTCACTACACGT")
LINKER2 = "TCGCTGACACGATCG"
RC = str.maketrans("ACGTN", "TGCAN")


def whitelist(path: str) -> dict[str, tuple[int, int]]:
    result: dict[str, tuple[int, int]] = {}
    with open(path, encoding="utf-8", errors="replace", newline="") as handle:
        rows = list(csv.reader(handle))
    for row in rows:
        if len(row) < 3:
            continue
        try:
            x, y = int(float(row[0])), int(float(row[1]))
            bc = row[2].strip().upper()
        except ValueError:
            continue
        if len(bc) == 12 and set(bc) <= set("ACGT"):
            result[bc] = (x, y)
    if not result:
        raise SystemExit(f"No 12 bp barcode rows found in {path}")
    return result


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else 999


def assign(raw: str, wl: dict[str, tuple[int, int]], max_mm: int) -> tuple[str, str, int, int, int]:
    if raw in wl:
        x, y = wl[raw]
        return "exact", raw, 0, x, y
    candidates = [(hamming(raw, bc), bc, xy) for bc, xy in wl.items() if hamming(raw, bc) <= max_mm]
    if not candidates:
        return "not_in_whitelist", ".", -1, -1, -1
    best = min(x[0] for x in candidates)
    winners = [x for x in candidates if x[0] == best]
    if len(winners) != 1:
        return "ambiguous", ".", best, -1, -1
    mm, bc, (x, y) = winners[0]
    return "rescued", bc, mm, x, y


def parse(seq: str, umi_len: int) -> tuple[str, str, str, str]:
    reverse = False
    if seq.find(ANCHOR) < 0:
        seq = seq.translate(RC)[::-1]
        reverse = True
    start = seq.find(ANCHOR)
    if start < 0:
        return "anchor_not_found", ".", ".", ""
    p = start + len(ANCHOR)
    b1 = seq[p:p + 4]
    p += 4
    linker1 = next((x for x in LINKER1 if seq.startswith(x, p)), "")
    if not linker1:
        return "linker1_not_found", ".", ".", ""
    p += len(linker1)
    b2 = seq[p:p + 4]
    p += 4
    if not seq.startswith(LINKER2, p):
        return "linker2_not_found", ".", ".", ""
    p += len(LINKER2)
    b3 = seq[p:p + 4]
    p += 4
    umi = seq[p:p + umi_len]
    p += umi_len
    while p < len(seq) and seq[p] == "T":
        p += 1
    if len(b1) != 4 or len(b2) != 4 or len(b3) != 4 or len(umi) != umi_len or len(seq) - p < 80:
        return "short_or_invalid", ".", ".", ""
    return "ok", b1 + b2 + b3, umi, seq[p:]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True); ap.add_argument("--assignments", required=True)
    ap.add_argument("--cdna", required=True); ap.add_argument("--summary", required=True)
    ap.add_argument("--max-reads", type=int, default=0,
                    help="Process at most this many reads; 0 means the complete FASTQ")
    args = ap.parse_args(); cfg = config(args.config); wl = whitelist(cfg["whitelist_csv"])
    counts: dict[str, int] = {}; kept = 0; rows = []; processed = 0
    cdna_path = ensure_parent(args.cdna)
    with gzip.open(cdna_path, "wt", encoding="utf-8") as cdna:
        for read_id, seq, qual in fastq(cfg["raw_fastq"]):
            if args.max_reads and processed >= args.max_reads:
                break
            processed += 1
            status, raw, umi, cdna_seq = parse(seq, int(cfg.get("umi_length", "10")))
            counts[status] = counts.get(status, 0) + 1
            if status != "ok":
                continue
            assign_status, bc, mm, x, y = assign(raw, wl, int(cfg.get("barcode_max_mismatch", "1")))
            counts[assign_status] = counts.get(assign_status, 0) + 1
            if assign_status not in {"exact", "rescued"}:
                continue
            kept += 1
            rows.append({"read_id": read_id, "barcode": bc, "umi": umi, "x": x, "y": y, "barcode_mm": mm, "structure_status": status})
            cdna.write(f"@{read_id}\n{cdna_seq}\n+\n{qual[-len(cdna_seq):]}\n")
    with gzip.open(ensure_parent(args.assignments), "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["read_id", "barcode", "umi", "x", "y", "barcode_mm", "structure_status"], delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    dump_json(args.summary, {"sample_id": "EXAMPLE_SAMPLE", "total_reads": processed, "kept_spatial_reads": kept, "counts": counts, "barcode_whitelist_size": len(wl), "structure": "C4L15C4L15C4U10T18", "read_limit": args.max_reads or None, "complete_fastq_scan": not bool(args.max_reads)})


if __name__ == "__main__":
    main()
