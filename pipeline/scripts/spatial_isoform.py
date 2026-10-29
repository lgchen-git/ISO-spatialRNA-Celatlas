"""Run the supported All/bin-level Spl-IsoFind Moran's I route."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--isoform", required=True, type=Path)
    p.add_argument("--metadata", required=True, type=Path)
    p.add_argument("--isoforms", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--status", required=True, type=Path)
    p.add_argument("--stats", required=True, type=Path)
    return p.parse_args()


def write_status(path: Path, status: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def empty_stats(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(columns=["isoform_id", "morans_I", "p_value", "q_value", "gene_id", "total_molecules", "positive_bins"]).to_csv(path, sep="\t", index=False)


def main() -> None:
    a = parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    try:
        matrix = pd.read_csv(a.isoform, sep="\t", compression="gzip", dtype={"isoform_id": str, "spatial_unit_id": str})
        iso = pd.read_csv(a.isoforms, sep="\t", dtype={"isoform_id": str})
        meta = pd.read_csv(a.metadata, sep="\t", dtype={"spatial_unit_id": str})
    except Exception as exc:
        empty_stats(a.stats)
        write_status(a.status, {"status": "BLOCKED", "route": "Spl-IsoFind All/bin-level exploratory Moran's I", "cell_type_specific": "NOT_SUPPORTED", "reason": f"Unable to read inputs: {exc}"})
        return
    required = {"isoform_id", "spatial_unit_id", "molecules"}
    if not required.issubset(matrix.columns) or not {"isoform_id", "gene_id"}.issubset(iso.columns):
        empty_stats(a.stats)
        write_status(a.status, {"status": "BLOCKED", "route": "Spl-IsoFind All/bin-level exploratory Moran's I", "cell_type_specific": "NOT_SUPPORTED", "reason": "Required isoform matrix or annotation columns are missing"})
        return
    matrix["molecules"] = pd.to_numeric(matrix["molecules"], errors="coerce").fillna(0)
    annotated = matrix[matrix["molecules"].gt(0)].merge(iso[["isoform_id", "gene_id"]].drop_duplicates("isoform_id"), on="isoform_id", how="inner")
    meta = meta.drop_duplicates("spatial_unit_id").copy()
    parsed = meta["spatial_unit_id"].astype(str).str.extract(r"^bin_(\d+)_(\d+)$")
    meta["x"] = parsed[0].astype(float) * 50 + 25
    meta["y"] = parsed[1].astype(float) * 50 + 25
    meta = meta.dropna(subset=["x", "y"])
    summary = annotated.groupby(["gene_id", "isoform_id"], as_index=False).agg(total_molecules=("molecules", "sum"), positive_bins=("spatial_unit_id", "nunique"))
    counts = summary.groupby("gene_id")["isoform_id"].nunique().rename("gene_isoform_count")
    selected = summary.join(counts, on="gene_id")
    selected = selected[(selected["gene_isoform_count"] >= 2) & (selected["positive_bins"] >= 50)].nlargest(50, ["total_molecules", "positive_bins"])
    if selected.empty or meta.empty:
        empty_stats(a.stats)
        write_status(a.status, {"status": "BLOCKED", "route": "Spl-IsoFind All/bin-level exploratory Moran's I", "cell_type_specific": "NOT_SUPPORTED", "reason": "No eligible multi-isoform gene or valid 50 µm spatial bins"})
        return
    denominator = annotated.groupby(["gene_id", "spatial_unit_id"], as_index=False)["molecules"].sum().rename(columns={"molecules": "gene_molecules"})
    long = annotated[annotated["isoform_id"].isin(selected["isoform_id"])].merge(denominator, on=["gene_id", "spatial_unit_id"], how="left")
    long["psi"] = long["molecules"] / long["gene_molecules"].replace(0, np.nan)
    wide = long.pivot_table(index="spatial_unit_id", columns="isoform_id", values="psi", aggfunc="sum", fill_value=0).reindex(meta["spatial_unit_id"], fill_value=0)
    columns = ["isoform_id", "morans_I", "p_value", "q_value", "gene_id", "total_molecules", "positive_bins"]
    try:
        from SplIsoFind import sv
        coords = meta.set_index(pd.Index(range(len(meta))))[["x", "y"]]
        m_i, pval, qval = sv.moransI(wide, coords, nperm=999, k=10, mincells=50, mincellspergroup=20, imb=.05, celltypes=["All"], x_coord="x", y_coord="y", n_jobs=1, seed=0)
        stats = pd.DataFrame({"isoform_id": m_i.index.astype(str), "morans_I": m_i["All"].values, "p_value": pval["All"].values, "q_value": qval["All"].values})
        stats = stats.merge(selected[["isoform_id", "gene_id", "total_molecules", "positive_bins"]], on="isoform_id", how="left")
        status = {"status": "PASS", "route": "Spl-IsoFind All/bin-level exploratory Moran's I", "spatial_unit": "50 µm bin", "cell_type_specific": "NOT_SUPPORTED", "candidate_isoforms": int(wide.shape[1]), "moransI_nonnull": int(stats["morans_I"].notna().sum()), "nperm": 999, "k": 10}
    except Exception as exc:
        stats = pd.DataFrame(columns=columns)
        status = {"status": "BLOCKED", "route": "Spl-IsoFind All/bin-level exploratory Moran's I", "spatial_unit": "50 µm bin", "cell_type_specific": "NOT_SUPPORTED", "reason": f"Spl-IsoFind execution unavailable: {exc}"}
    a.stats.parent.mkdir(parents=True, exist_ok=True)
    stats.reindex(columns=columns).to_csv(a.stats, sep="\t", index=False)
    write_status(a.status, status)


if __name__ == "__main__":
    main()
