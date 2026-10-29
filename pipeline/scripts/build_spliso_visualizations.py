#!/usr/bin/env python3
"""Build isolated Spl-IsoQuant/Spl-IsoFind visualizations from EXAMPLE_SAMPLE results.

The input result tree is treated as read-only.  Spl-IsoQuant figures describe
the existing caller output (v4.0.0).  Spl-IsoFind is run only on gene-relative
PSI reconstructed at 50 um bins and only for the All/bin-level route; this is
an exploratory spatial autocorrelation analysis, not a cell-type-specific SVI.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


NPG = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"]
GENERAL = LinearSegmentedColormap.from_list("npg_general", ["#E5E7EB", "#4DBBD5", "#3C5488", "#E64B35"])
ACTIVITY = LinearSegmentedColormap.from_list("npg_activity", ["#E5E7EB", "#91D1C2", "#00A087", "#DC0000"])


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--source-root", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--max-isoforms", type=int, default=50)
    p.add_argument("--nperm", type=int, default=999)
    p.add_argument("--k", type=int, default=10)
    p.add_argument("--stats", type=Path, default=None, help="precomputed Spl-IsoFind statistics")
    return p.parse_args()


def style(ax, title: str, xlabel: str = "", ylabel: str = "") -> None:
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)


def save(fig: plt.Figure, stem: Path, dpi: int) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def add_record(records: list[dict], stem: Path, figure_id: str, section: str, title: str, sources: list[Path], limitation: str = "") -> None:
    records.append({
        "figure_id": figure_id,
        "section": section,
        "title": title,
        "png": str(stem.with_suffix(".png").resolve()),
        "pdf": str(stem.with_suffix(".pdf").resolve()),
        "source_paths": ";".join(str(x.resolve()) for x in sources if x.exists()),
        "generation_mode": "regenerated_from_validated_tables",
        "limitation": limitation,
    })


def spatial_frame(meta: pd.DataFrame) -> pd.DataFrame:
    required = {"spatial_unit_id", "x", "y"}
    missing = required.difference(meta.columns)
    if missing:
        raise ValueError(f"spatial metadata missing columns: {sorted(missing)}")
    meta = meta.drop_duplicates("spatial_unit_id").copy()
    parsed = meta["spatial_unit_id"].astype(str).str.extract(r"^bin_(\d+)_(\d+)$")
    parsed_ok = parsed.notna().all(axis=1)
    meta["x"] = pd.to_numeric(meta["x"], errors="coerce")
    meta["y"] = pd.to_numeric(meta["y"], errors="coerce")
    # spatial_metadata contains read-level rows whose x/y values vary within a
    # bin.  The unit ID is the authoritative 50 um grid coordinate.
    parsed_x = parsed[0].astype(float) * 50 + 25
    parsed_y = parsed[1].astype(float) * 50 + 25
    meta["x"] = meta["x"].where(~parsed_ok, parsed_x)
    meta["y"] = meta["y"].where(~parsed_ok, parsed_y)
    meta = meta.dropna(subset=["x", "y"])
    if meta.empty:
        raise ValueError("no finite spatial coordinates")
    return meta


def plot_spatial_patches(ax, frame: pd.DataFrame, value: str, title: str, cmap, vmin=None, vmax=None, bin_size: float = 50) -> ScalarMappable:
    vals = pd.to_numeric(frame[value], errors="coerce")
    ok = vals.notna()
    vals = vals.loc[ok]
    dat = frame.loc[ok]
    if vmin is None or vmax is None:
        vmin, vmax = np.percentile(vals, [2, 98]) if len(vals) else (0, 1)
        if vmax <= vmin:
            vmax = vmin + 1
    patches = [Rectangle((x - bin_size / 2, y - bin_size / 2), bin_size, bin_size) for x, y in zip(dat["x"], dat["y"])]
    collection = PatchCollection(patches, cmap=cmap, norm=Normalize(vmin, vmax), linewidth=0, rasterized=True)
    collection.set_array(vals.to_numpy())
    ax.add_collection(collection)
    ax.set_xlim(frame["x"].min() - bin_size / 2, frame["x"].max() + bin_size / 2)
    ax.set_ylim(frame["y"].min() - bin_size / 2, frame["y"].max() + bin_size / 2)
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="box")
    style(ax, title, "X (µm; right)", "Y (µm; down)")
    return collection


def gene_scale(frame: pd.DataFrame, gene_id: str, value: str, default=(0.0, 1.0)) -> tuple[float, float]:
    values = pd.to_numeric(frame.loc[frame["gene_id"].eq(gene_id), value], errors="coerce").dropna()
    if values.empty:
        return default
    low, high = np.percentile(values, [2, 98])
    if high <= low:
        high = low + 1
    return float(low), float(high)


def load_inputs(source: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    iso_path = source / "04_isoform_calling/all_isoforms.tsv"
    matrix_path = source / "05_matrices/isoform_x_spatial.tsv.gz"
    meta_path = source / "05_matrices/spatial_metadata.tsv"
    for path in (iso_path, matrix_path, meta_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    iso = pd.read_csv(iso_path, sep="\t")
    iso["isoform_id"] = iso["isoform_id"].astype(str)
    iso = iso[iso["isoform_id"].ne(".")].copy()
    iso["read_support"] = pd.to_numeric(iso["read_support"], errors="coerce").fillna(0)
    iso = iso.sort_values("read_support", ascending=False).drop_duplicates("isoform_id")
    matrix = pd.read_csv(matrix_path, sep="\t", compression="gzip", dtype={"isoform_id": str, "spatial_unit_id": str})
    matrix["molecules"] = pd.to_numeric(matrix["molecules"], errors="coerce").fillna(0)
    matrix = matrix[matrix["molecules"].gt(0)].copy()
    meta = spatial_frame(pd.read_csv(meta_path, sep="\t", dtype={"spatial_unit_id": str}))
    annotated = matrix.merge(iso[["isoform_id", "gene_id", "classification"]], on="isoform_id", how="inner")
    if annotated.empty:
        raise ValueError("no isoform matrix rows can be joined to all_isoforms.tsv")
    return iso, matrix, meta, annotated


def splisoquant_figures(records, figdir: Path, iso: pd.DataFrame, annotated: pd.DataFrame, meta: pd.DataFrame, dpi: int, sources: list[Path]) -> None:
    classes = iso["classification"].fillna("unknown").astype(str).value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.bar(classes.index, classes.values, color=NPG[: len(classes)])
    ax.tick_params(axis="x", rotation=45)
    style(ax, "Spl-IsoQuant isoform classification (existing output v4.0.0)", ylabel="Isoform records")
    save(fig, figdir / "splisoquant_isoform_classification", dpi)
    add_record(records, figdir / "splisoquant_isoform_classification", "splisoquant_isoform_classification", "Spl-IsoQuant", "Isoform classification", sources)

    support = iso[["read_support", "umi_support", "spatial_unit_support"]].apply(pd.to_numeric, errors="coerce").dropna()
    support = support[(support > 0).all(axis=1)]
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    sc = ax.scatter(support["umi_support"], support["read_support"], c=support["spatial_unit_support"], s=8, alpha=.4, linewidths=0, cmap=GENERAL, rasterized=True)
    ax.set_xscale("log"); ax.set_yscale("log")
    style(ax, "Spl-IsoQuant isoform support", "UMI support", "Read support")
    fig.colorbar(sc, ax=ax, label="Spatial bin support")
    save(fig, figdir / "splisoquant_support_vs_space", dpi)
    add_record(records, figdir / "splisoquant_support_vs_space", "splisoquant_support_vs_space", "Spl-IsoQuant", "Read, UMI and spatial support", sources)

    top = annotated.groupby(["gene_id", "isoform_id"], as_index=False)["molecules"].sum().nlargest(6, "molecules")
    values = annotated[annotated["isoform_id"].isin(top["isoform_id"])].groupby(["spatial_unit_id", "isoform_id"], as_index=False)["molecules"].sum()
    fig, axes = plt.subplots(2, 3, figsize=(11, 7.2), constrained_layout=True)
    for ax, row in zip(axes.ravel(), top.to_dict("records")):
        ident, gene = row["isoform_id"], row["gene_id"]
        panel = meta[["spatial_unit_id", "x", "y"]].merge(values[values["isoform_id"].eq(ident)][["spatial_unit_id", "molecules"]], on="spatial_unit_id", how="left")
        panel["molecules"] = panel["molecules"].fillna(0)
        low, high = gene_scale(annotated, gene, "molecules", (0, 1))
        sm = plot_spatial_patches(ax, panel, "molecules", f"{gene}\n{ident}", GENERAL, low, high)
        cb = fig.colorbar(sm, ax=ax, fraction=.035, pad=.015)
        cb.set_label("molecules/bin", fontsize=7)
        cb.ax.tick_params(labelsize=7)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Spl-IsoQuant isoform spatial maps (50 µm bin; coordinate-only)", x=.02, ha="left", fontsize=12, fontweight="bold")
    save(fig, figdir / "splisoquant_top_isoform_spatial", dpi)
    add_record(records, figdir / "splisoquant_top_isoform_spatial", "splisoquant_top_isoform_spatial", "Spl-IsoQuant", "Top isoform spatial maps", sources, "No tissue image; measured bin values only.")


def build_psi(iso: pd.DataFrame, annotated: pd.DataFrame, meta: pd.DataFrame, max_isoforms: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    iso_summary = annotated.groupby(["gene_id", "isoform_id"], as_index=False).agg(total_molecules=("molecules", "sum"), positive_bins=("spatial_unit_id", "nunique"))
    gene_counts = iso_summary.groupby("gene_id")["isoform_id"].nunique().rename("gene_isoform_count")
    candidates = iso_summary.join(gene_counts, on="gene_id")
    candidates = candidates[(candidates["gene_isoform_count"] >= 2) & (candidates["positive_bins"] >= 50)]
    selected = candidates.nlargest(max_isoforms, ["total_molecules", "positive_bins"]).copy()
    if selected.empty:
        raise ValueError("no gene with at least two isoforms and 50 positive spatial bins")
    denominator = annotated.groupby(["gene_id", "spatial_unit_id"], as_index=False)["molecules"].sum().rename(columns={"molecules": "gene_molecules"})
    long = annotated[annotated["isoform_id"].isin(selected["isoform_id"])].merge(denominator, on=["gene_id", "spatial_unit_id"], how="left")
    long["psi"] = long["molecules"] / long["gene_molecules"].replace(0, np.nan)
    wide = long.pivot_table(index="spatial_unit_id", columns="isoform_id", values="psi", aggfunc="sum", fill_value=0)
    units = meta.copy()
    wide = wide.reindex(units["spatial_unit_id"], fill_value=0)
    return selected, wide.reset_index(), units


def splisofind_figures(records, figdir: Path, iso: pd.DataFrame, annotated: pd.DataFrame, meta: pd.DataFrame, dpi: int, max_isoforms: int, nperm: int, k: int, sources: list[Path], stats_path: Path | None = None) -> dict:
    selected, wide, labels = build_psi(iso, annotated, meta, max_isoforms)
    feature_ids = [str(x) for x in wide.columns if x != "spatial_unit_id"]
    x = wide[feature_ids].astype(float)
    labels_for_find = labels.set_index(pd.Index(range(len(labels))))[["x", "y"]]
    if stats_path is not None and stats_path.is_file():
        stats = pd.read_csv(stats_path, sep="\t")
    else:
        try:
            from SplIsoFind import sv
        except Exception as exc:
            raise RuntimeError("Spl-IsoFind is required for All/bin-level Moran's I") from exc
        m_i, pval, qval = sv.moransI(x, labels_for_find, nperm=nperm, k=k, mincells=50, mincellspergroup=20, imb=.05, celltypes=["All"], x_coord="x", y_coord="y", n_jobs=1, seed=0)
        stats = pd.DataFrame({"isoform_id": m_i.index.astype(str), "morans_I": m_i["All"].values, "p_value": pval["All"].values, "q_value": qval["All"].values})
        stats = stats.merge(selected[["isoform_id", "gene_id", "total_molecules", "positive_bins"]], on="isoform_id", how="left")
        stats.to_csv(figdir.parent / "splisofind_all_bin_moransI.tsv", sep="\t", index=False)

    ranked = stats.dropna(subset=["morans_I"]).copy()
    ranked["neglog10_q"] = -np.log10(ranked["q_value"].clip(lower=1e-300))
    ranked = ranked.sort_values(["q_value", "morans_I"], ascending=[True, False])
    top_rank = ranked.head(30).sort_values("morans_I")
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.barh(top_rank["isoform_id"], top_rank["morans_I"], color="#3C5488")
    style(ax, "Spl-IsoFind All/bin-level spatial autocorrelation", "Moran's I", "Isoform")
    save(fig, figdir / "splisofind_moransI_rank", dpi)
    add_record(records, figdir / "splisofind_moransI_rank", "splisofind_moransI_rank", "Spl-IsoFind", "All/bin-level Moran's I ranking", sources, "Exploratory bin-level PSI; no cell-type-specific SVI.")

    top_ids = ranked.head(6)["isoform_id"].tolist()
    selected_gene = selected.set_index("isoform_id")["gene_id"].astype(str).to_dict()
    gene_scales = {}
    for gene in sorted(set(selected_gene.values())):
        gene_features = [ident for ident in feature_ids if selected_gene.get(ident) == gene]
        gene_values = x[gene_features].to_numpy().ravel()
        finite = gene_values[np.isfinite(gene_values)]
        if finite.size:
            low, high = np.percentile(finite, [2, 98])
            gene_scales[gene] = (float(low), float(high if high > low else low + 1))
        else:
            gene_scales[gene] = (0.0, 1.0)
    xplot = x.copy()
    nplots = max(1, len(top_ids))
    ncols = min(3, nplots)
    nrows = int(np.ceil(nplots / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.7 * nrows), squeeze=False, constrained_layout=True)
    axes = axes.ravel()
    cmap = ACTIVITY
    spatial_values = labels[["spatial_unit_id", "x", "y"]].copy()
    for ax, ident in zip(axes.ravel(), top_ids):
        gene = selected_gene[ident]
        panel = spatial_values.copy()
        panel["psi"] = x[ident].to_numpy()
        low, high = gene_scales[gene]
        sm = plot_spatial_patches(ax, panel, "psi", f"{gene}\n{ident}", cmap, low, high)
        cb = fig.colorbar(sm, ax=ax, fraction=.035, pad=.015)
        cb.set_label("gene-specific PSI", fontsize=7)
        cb.ax.tick_params(labelsize=7)
        ax.set_xticks([]); ax.set_yticks([])
    for ax in axes[len(top_ids):]:
        ax.set_axis_off()
    fig.suptitle("Spl-IsoFind top spatially variable isoforms (All/bin-level)", x=.02, ha="left", fontsize=12, fontweight="bold")
    save(fig, figdir / "splisofind_top_spatial_psi", dpi)
    add_record(records, figdir / "splisofind_top_spatial_psi", "splisofind_top_spatial_psi", "Spl-IsoFind", "Top spatial isoform PSI maps", sources, "Spl-IsoFind hexbin display of measured 50 µm bins; no H&E image.")

    burden = meta[["spatial_unit_id", "x", "y"]].copy().merge(annotated.groupby("spatial_unit_id")["isoform_id"].nunique().rename("detected_isoforms"), on="spatial_unit_id", how="left")
    burden["detected_isoforms"] = burden["detected_isoforms"].fillna(0)
    fig, ax = plt.subplots(figsize=(7.4, 6.1))
    sm = plot_spatial_patches(ax, burden, "detected_isoforms", "Observed isoform burden per 50 µm bin", ACTIVITY)
    fig.colorbar(sm, ax=ax, label="Distinct observed isoforms")
    save(fig, figdir / "splisofind_isoform_burden_spatial", dpi)
    add_record(records, figdir / "splisofind_isoform_burden_spatial", "splisofind_isoform_burden_spatial", "Spl-IsoFind", "Spatial isoform burden", sources, "Descriptive bin-level burden, not a cell-level SVI call.")
    return {"candidate_isoforms": len(feature_ids), "moransI_nonnull": int(ranked.shape[0]), "top_spatial_isoforms": top_ids}


def main() -> None:
    a = parse_args()
    out = a.outdir
    figdir = out / "01_figures"
    figdir.mkdir(parents=True, exist_ok=True)
    iso, matrix, meta, annotated = load_inputs(a.source_root)
    sources = [a.source_root / "04_isoform_calling/all_isoforms.tsv", a.source_root / "05_matrices/isoform_x_spatial.tsv.gz", a.source_root / "05_matrices/spatial_metadata.tsv", a.source_root / "04_isoform_calling/method_status.json"]
    records: list[dict] = []
    splisoquant_figures(records, figdir, iso, annotated, meta, a.dpi, sources)
    find_status = splisofind_figures(records, figdir, iso, annotated, meta, a.dpi, a.max_isoforms, a.nperm, a.k, sources, a.stats)
    manifest = out / "figure_manifest.tsv"
    columns = ["figure_id", "section", "title", "png", "pdf", "source_paths", "generation_mode", "limitation"]
    manifest.write_text("\t".join(columns) + "\n" + "\n".join("\t".join(str(row[c]) for c in columns) for row in records) + "\n", encoding="utf-8")
    status = {"status": "PASS", "figure_count": len(records), "spatial_unit": "50 µm bin", "coordinate_frame": "x right; y down", "h_and_e": False, "splisoquant_source_version": "4.0.0", "splisofind_route": "All/bin-level exploratory Moran's I", "splisofind_cell_type_specific": "NOT_SUPPORTED", "splisofind": find_status, "input_matrix_rows": int(matrix.shape[0]), "input_spatial_units": int(meta.shape[0])}
    (out / "visualization_status.json").write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "figure_count": len(records), "outdir": str(out.resolve()), **find_status}, ensure_ascii=False))


if __name__ == "__main__":
    main()
