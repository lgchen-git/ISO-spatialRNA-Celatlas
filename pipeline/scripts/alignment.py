from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path

from common import config, ensure_parent, open_text


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True); ap.add_argument("--assignments", required=True)
    ap.add_argument("--cdna", default=""); ap.add_argument("--bam", required=True)
    ap.add_argument("--bai", required=True); ap.add_argument("--table", required=True)
    ap.add_argument("--upstream-bam", default=""); ap.add_argument("--upstream-bai", default="")
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args(); c = config(a.config)
    if a.upstream_bam:
        if not Path(a.upstream_bam).is_file() or not Path(a.upstream_bai).is_file():
            raise SystemExit("Configured upstream BAM/BAI is incomplete")
        ensure_parent(a.bam); shutil.copy2(a.upstream_bam, a.bam); shutil.copy2(a.upstream_bai, a.bai)
        try:
            import pandas as pd
        except ModuleNotFoundError as exc:
            raise SystemExit("alignment table requires pandas and pyarrow") from exc
        pd.DataFrame([{"source_bam": str(Path(a.upstream_bam).resolve()), "status": "validated_upstream_alignment", "min_mapq": c.get("min_mapq", "20")}]).to_parquet(ensure_parent(a.table), index=False)
        return
    for binary in ("minimap2", "samtools"):
        if subprocess.run(["bash", "-lc", f"command -v {binary}"], capture_output=True).returncode:
            raise SystemExit(f"Required executable unavailable: {binary}")
    mmi = Path(c.get("minimap2_index", ""))
    if not mmi.is_file():
        mmi = Path(a.bam).with_suffix(".mmi")
        subprocess.run(["minimap2", "-d", str(mmi), c["genome_fa"]], check=True)
    unsorted = Path(a.bam).with_suffix(".unsorted.bam")
    ensure_parent(a.bam)
    p = subprocess.Popen(["minimap2", "-ax", "splice", "-uf", "--secondary=no", "-t", str(a.threads), str(mmi), a.cdna], stdout=subprocess.PIPE)
    assert p.stdout is not None
    subprocess.run(["samtools", "view", "-b", "-o", str(unsorted), "-"], stdin=p.stdout, check=True)
    if p.wait():
        raise SystemExit("minimap2 failed")
    subprocess.run(["samtools", "sort", "-@", str(a.threads), "-o", a.bam, str(unsorted)], check=True)
    subprocess.run(["samtools", "index", "-@", str(a.threads), a.bam], check=True)
    unsorted.unlink(missing_ok=True)
    try:
        import pandas as pd
        import pysam
    except ModuleNotFoundError as exc:
        raise SystemExit("alignment table requires pandas and pysam") from exc
    assignments = {}
    with open_text(a.assignments) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            assignments[row["read_id"]] = row
    rows = []
    with pysam.AlignmentFile(a.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            x = assignments.get(read.query_name, {})
            rows.append({"read_id": read.query_name, "barcode": x.get("barcode", "."), "umi": x.get("umi", "."), "spatial_x": x.get("x", -1), "spatial_y": x.get("y", -1), "chromosome": read.reference_name, "start": read.reference_start, "end": read.reference_end, "strand": "-" if read.is_reverse else "+", "mapping_quality": read.mapping_quality, "read_length": read.query_length or 0, "aligned_length": read.query_alignment_length or 0, "gene_id": "unassigned", "transcript_id": "unassigned", "isoform_id": "unassigned", "spliced": int(read.cigarstring and "N" in read.cigarstring), "full_length_status": "uncertain"})
    if not rows:
        raise SystemExit("No primary mapped reads were produced")
    pd.DataFrame(rows).to_parquet(ensure_parent(a.table), index=False)


if __name__ == "__main__":
    main()
