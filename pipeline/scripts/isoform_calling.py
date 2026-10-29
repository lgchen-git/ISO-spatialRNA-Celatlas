from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
from collections import defaultdict
from pathlib import Path

from common import config, dump_json, ensure_parent, open_text


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True); ap.add_argument("--bam", required=True)
    ap.add_argument("--assignments", required=True); ap.add_argument("--isoforms", required=True)
    ap.add_argument("--novel", required=True); ap.add_argument("--method", required=True)
    ap.add_argument("--reads", required=True); ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args(); c = config(a.config); out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    executable = c.get("splisoquant_executable", "splisoquant")
    cmd = [executable, "--reference", c["genome_fa"], "--genedb", c["gtf"], "--bam", a.bam,
           "--barcoded_reads", a.assignments, "--mode", "custom_sc", "--analysis", "quantification",
           "transcript_discovery", "--data_type", "nanopore", "--output", str(out), "--prefix", "RBC04",
           "--threads", str(a.threads), "--large_output", "read_info", "read2transcripts",
           "--complete_genedb", "--force"]
    read_info_candidates = [out / "RBC04.read_info.tsv.gz", out / "RBC04" / "RBC04.read_info.tsv.gz"]
    read_info = next((path for path in read_info_candidates if path.is_file()), None)
    if read_info is None:
        try:
            subprocess.run(cmd, check=True)
        except FileNotFoundError as exc:
            dump_json(a.method, {"status": "BLOCKED", "reason": f"Spl-IsoQuant executable unavailable: {executable}"})
            raise SystemExit(str(exc))
        except subprocess.CalledProcessError as exc:
            dump_json(a.method, {"status": "FAILED", "reason": f"Spl-IsoQuant exit code {exc.returncode}", "command": cmd})
            raise
        read_info = next((path for path in read_info_candidates if path.is_file()), None)
    if read_info is None:
        raise SystemExit(f"Spl-IsoQuant did not produce any expected read_info file: {read_info_candidates}")
    assignment_by_read = {}
    with open_text(a.assignments) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            assignment_by_read[row["read_id"]] = row
    groups = defaultdict(lambda: {"reads": 0, "umis": set(), "spatial": set(), "classification": "."})
    rows = []
    with gzip.open(ensure_parent(a.reads), "wt", encoding="utf-8", newline="") as handle:
        with gzip.open(read_info, "rt", encoding="utf-8", errors="replace") as src:
            reader = csv.DictReader((line for line in src if not line.startswith("#")), delimiter="\t")
            for row in reader:
                rid = row.get("read_id", ""); bc = assignment_by_read.get(rid)
                if bc is None or row.get("gene_id", ".") in {"", ".", "*"}:
                    continue
                gene = row.get("gene_id", "."); transcript = row.get("isoform_id", ".")
                info = groups[(gene, transcript)]; info["reads"] += 1
                info["classification"] = row.get("classification", "."); info["umis"].add(bc.get("umi", ".")); info["spatial"].add((bc.get("x", "-1"), bc.get("y", "-1")))
                rows.append({"read_id": rid, "barcode": bc.get("barcode", "."), "umi": bc.get("umi", "."), "x": bc.get("x", "-1"), "y": bc.get("y", "-1"), "gene_id": gene, "transcript_id": transcript, "isoform_id": transcript, "classification": row.get("classification", "."), "strand": row.get("strand", "."), "chromosome": row.get("chr", "."), "polyA": row.get("polyA", "."), "exons": row.get("exons", ".")})
        fields = list(rows[0]) if rows else ["read_id", "barcode", "umi", "x", "y", "gene_id", "transcript_id", "isoform_id", "classification", "strand", "chromosome", "polyA", "exons"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    keys = ["gene_id", "isoform_id", "classification", "read_support", "umi_support", "spatial_unit_support"]
    iso_rows = [{"gene_id": gene, "isoform_id": transcript, "classification": info["classification"], "read_support": info["reads"], "umi_support": len(info["umis"]), "spatial_unit_support": len(info["spatial"])} for (gene, transcript), info in groups.items()]
    def write_tsv(path: str, data: list[dict[str, object]]) -> None:
        Path(ensure_parent(path)).write_text("\t".join(keys) + "\n" + "\n".join("\t".join(str(r[k]) for k in keys) for r in data) + ("\n" if data else ""), encoding="utf-8")
    write_tsv(a.isoforms, iso_rows)
    novel = [r for r in iso_rows if r["classification"] not in {"complete_match", "full-splice_match", "incomplete_splice_match", "."} and r["read_support"] >= int(c.get("min_isoform_reads", "5")) and r["umi_support"] >= int(c.get("min_isoform_umis", "3")) and r["spatial_unit_support"] >= int(c.get("min_isoform_spatial_units", "3"))]
    write_tsv(a.novel, novel)
    dump_json(a.method, {"status": "PASS", "caller": "Spl-IsoQuant", "version": "4.0.0", "read_info": str(read_info), "assigned_reads": len(rows), "isoforms": len(iso_rows), "high_confidence_novel": len(novel)})


if __name__ == "__main__":
    main()
