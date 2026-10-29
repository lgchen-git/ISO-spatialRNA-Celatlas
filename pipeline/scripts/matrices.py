from __future__ import annotations

import argparse
import gzip
from pathlib import Path

from common import config, ensure_parent


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True); ap.add_argument("--reads", required=True); ap.add_argument("--caller-reads", required=True); ap.add_argument("--isoforms", required=True)
    ap.add_argument("--gene", required=True); ap.add_argument("--transcript", required=True); ap.add_argument("--isoform", required=True); ap.add_argument("--metadata", required=True)
    a = ap.parse_args(); c = config(a.config)
    try:
        import pandas as pd
    except ModuleNotFoundError as exc:
        raise SystemExit("matrix generation requires pandas and pyarrow") from exc
    # The alignment table supplies mapping/QC fields; the caller table supplies
    # validated gene/transcript/isoform assignments and spatial identity.
    reads = pd.read_csv(a.caller_reads, sep="\t", compression="gzip")
    if reads.empty:
        raise SystemExit("read-level table is empty")
    # The caller wrapper emits the canonical short coordinate names x/y.
    # Accept the older spatial_x/spatial_y names for compatibility with prior
    # intermediate tables.
    x_col = "x" if "x" in reads.columns else "spatial_x"
    y_col = "y" if "y" in reads.columns else "spatial_y"
    reads["spatial_unit_id"] = reads.apply(lambda r: f"bin_{int(float(r[x_col]))//int(c.get('bin_size','50'))}_{int(float(r[y_col]))//int(c.get('bin_size','50'))}", axis=1)
    reads["molecule_id"] = reads["barcode"].astype(str) + ":" + reads["umi"].astype(str) + ":" + reads["gene_id"].astype(str)
    reads = reads.drop_duplicates("molecule_id")
    valid = reads[reads["gene_id"].ne("unassigned")].copy()
    if valid.empty:
        raise SystemExit("No gene-assigned molecules; matrix generation stopped")
    def matrix(feature: str, path: str) -> None:
        tab = valid.groupby([feature, "spatial_unit_id"], as_index=False).size().rename(columns={"size": "molecules"})
        tab.to_csv(ensure_parent(path), sep="\t", index=False, compression="gzip")
    matrix("gene_id", a.gene); matrix("transcript_id", a.transcript); matrix("isoform_id", a.isoform)
    meta = valid[["spatial_unit_id", x_col, y_col]].drop_duplicates().copy()
    meta = meta.rename(columns={x_col: "x", y_col: "y"})
    meta["cluster"] = "unknown"; meta["cell_type"] = "NA"; meta["region"] = "NA"
    meta.to_csv(ensure_parent(a.metadata), sep="\t", index=False)


if __name__ == "__main__":
    main()
