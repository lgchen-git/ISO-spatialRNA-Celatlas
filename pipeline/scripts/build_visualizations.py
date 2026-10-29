#!/usr/bin/env python3
"""Build the legacy module figures and new spatial isoform figures.

The legacy result tree is read-only evidence for the previously completed
modules.  Figures that can be regenerated from tables are regenerated here;
the two raw per-read QC plots are carried forward because their source rows
are not part of the visible result contract.
"""

from __future__ import annotations

import argparse
import base64
import json
import math
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle


NPG = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#6F99AD", "#8A7AA9", "#C9A227", "#6E7F80", "#D4537E", "#378ADD"]
GENERAL = LinearSegmentedColormap.from_list("npg_general", ["#E5E7EB", "#4DBBD5", "#3C5488", "#E64B35"])
ACTIVITY = LinearSegmentedColormap.from_list("npg_activity", ["#E5E7EB", "#91D1C2", "#00A087", "#DC0000"])
SIGNED = LinearSegmentedColormap.from_list("npg_signed", ["#3C5488", "#F7F7F7", "#E64B35"])


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--source-root", required=True, type=Path)
    p.add_argument("--legacy-root", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def style(ax, title: str, xlabel: str = "", ylabel: str = "") -> None:
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)


def save(fig, stem: Path, dpi: int) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def json_file(path: Path) -> dict:
    try: return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError): return {}


def metrics(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    return dict(zip(frame.iloc[:, 0], frame.iloc[:, 1]))


def record(records, stem: Path, figure_id: str, section: str, title: str, sources: list[Path], mode: str, limitation: str = "") -> None:
    records.append({"figure_id": figure_id, "section": section, "title": title, "png": str(stem.with_suffix(".png").resolve()), "pdf": str(stem.with_suffix(".pdf").resolve()), "source_paths": ";".join(str(p.resolve()) for p in sources if p.exists()), "generation_mode": mode, "limitation": limitation})


def bar(records, stem, fid, section, title, labels, values, sources, colors=None, ylabel="Count", dpi=300, rotation=0, log=False, mode="regenerated_from_validated_tables"):
    values = pd.to_numeric(pd.Series(values), errors="coerce").fillna(0).tolist()
    fig, ax = plt.subplots(figsize=(7.2, 4.7)); colors = colors or NPG[:len(labels)]
    ax.bar(labels, values, color=colors); ax.tick_params(axis="x", rotation=rotation)
    if log: ax.set_yscale("log")
    style(ax, title, ylabel=ylabel)
    for i, v in enumerate(values):
        if pd.notna(v) and float(v) > 0: ax.text(i, float(v), f"{float(v):,.0f}", ha="center", va="bottom", fontsize=7)
    save(fig, stem, dpi); record(records, stem, fid, section, title, sources, mode)


def spatial(records, stem, fid, section, title, frame, value, sources, dpi, cmap=GENERAL, categorical=False, limitation="No H&E image; coordinates are 50 µm Bin coordinates."):
    fig, ax = plt.subplots(figsize=(7.2, 6.1)); frame = frame.copy()
    if categorical:
        levels = sorted(frame[value].dropna().astype(str).unique(), key=lambda x: int(x) if x.isdigit() else x)
        for i, level in enumerate(levels):
            sub = frame[frame[value].astype(str) == level]
            ax.scatter(sub["bin_x"], sub["bin_y"], s=5, alpha=.86, linewidths=0, color=NPG[i % len(NPG)], label=f"C{level}", rasterized=True)
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", ncol=2, fontsize=7)
    else:
        x = pd.to_numeric(frame[value], errors="coerce"); ok = x.notna(); frame = frame.loc[ok].copy(); x = x.loc[ok]
        low, high = np.percentile(x, [2, 98]) if len(x) else (0, 1)
        if high <= low: high = low + 1
        sc = ax.scatter(frame["bin_x"], frame["bin_y"], c=x, s=5, marker="s", cmap=cmap, vmin=low, vmax=high, linewidths=0, rasterized=True)
        fig.colorbar(sc, ax=ax, label=f"{value} (display 2nd–98th percentile)")
    ax.invert_yaxis(); ax.set_aspect("equal", adjustable="datalim"); style(ax, title, "X (µm; right)", "Y (µm; down)"); save(fig, stem, dpi)
    record(records, stem, fid, section, title, sources, "regenerated_from_validated_tables", limitation)


def copy_baseline(records, source: Path, stem: Path, fid: str, section: str, title: str, dpi: int) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, stem.with_suffix(".png"))
    fig, ax = plt.subplots(figsize=(7, 4)); ax.imshow(plt.imread(source)); ax.axis("off"); fig.savefig(stem.with_suffix(".pdf"), dpi=dpi, bbox_inches="tight", facecolor="white"); plt.close(fig)
    record(records, stem, fid, section, title, [source], "reused_validated_baseline_plot", "原始逐 read 数据未纳入当前结果树；保留旧项目已验证图形资产。")


def parse_spatial(legacy: Path) -> pd.DataFrame:
    path = legacy / "06_whole_transcriptome/01_tables/long_read_bin50_clusters.tsv"
    return pd.read_csv(path, sep="\t")


def build_isoform(records, source: Path, figdir: Path, dpi: int) -> None:
    iso = pd.read_csv(source / "04_isoform_calling/all_isoforms.tsv", sep="\t")
    iso["classification"] = iso["classification"].fillna("unknown").astype(str)
    counts = iso["classification"].value_counts().sort_values()
    bar(records, figdir / "isoform_class_distribution", "isoform_class_distribution", "升级全长 isoform", "Spl-IsoQuant isoform classification", counts.index.tolist(), counts.values.tolist(), [source / "04_isoform_calling/all_isoforms.tsv"], dpi=dpi, rotation=45)
    support = iso[["read_support", "umi_support", "spatial_unit_support"]].apply(pd.to_numeric, errors="coerce")
    fig, ax = plt.subplots(figsize=(7.4, 5.2)); ax.scatter(support["umi_support"], support["read_support"], s=8, alpha=.35, c=support["spatial_unit_support"], cmap=GENERAL, linewidths=0, rasterized=True); ax.set_xscale("log"); ax.set_yscale("log"); style(ax, "Isoform read and UMI support", "UMI support", "Read support"); fig.colorbar(ax.collections[0], ax=ax, label="Spatial Bin support"); save(fig, figdir / "isoform_support_scatter", dpi); record(records, figdir / "isoform_support_scatter", "isoform_support_scatter", "升级全长 isoform", "Isoform read/UMI/spatial support", [source / "04_isoform_calling/all_isoforms.tsv"], "regenerated_from_current_isoform_table")
    iso2 = iso[iso["isoform_id"].astype(str).ne(".")].copy(); iso2["read_support"] = pd.to_numeric(iso2["read_support"], errors="coerce"); top = iso2.groupby("isoform_id", as_index=False)["read_support"].sum().nlargest(12, "read_support")
    bar(records, figdir / "top_isoform_support", "top_isoform_support", "升级全长 isoform", "Top isoform read support", top["isoform_id"].astype(str).tolist(), top["read_support"].tolist(), [source / "04_isoform_calling/all_isoforms.tsv"], dpi=dpi, rotation=70, ylabel="Assigned reads")
    matrix_path = source / "05_matrices/isoform_x_spatial.tsv.gz"; meta_path = source / "05_matrices/spatial_metadata.tsv"
    matrix = pd.read_csv(matrix_path, sep="\t", compression="gzip"); meta = pd.read_csv(meta_path, sep="\t")
    matrix = matrix[matrix["isoform_id"].astype(str).isin(top["isoform_id"].astype(str))].merge(meta[["spatial_unit_id", "x", "y"]], on="spatial_unit_id", how="left")
    if not matrix.empty:
        matrix["bin_x"] = (pd.to_numeric(matrix["x"]) // 50) * 50 + 25; matrix["bin_y"] = (pd.to_numeric(matrix["y"]) // 50) * 50 + 25
        fig, axes = plt.subplots(3, 4, figsize=(12, 8)); vmax = max(1, float(matrix["molecules"].quantile(.98)))
        for ax, ident in zip(axes.ravel(), top["isoform_id"].astype(str)):
            sub = matrix[matrix["isoform_id"].astype(str) == ident]; sc = ax.scatter(sub["bin_x"], sub["bin_y"], c=sub["molecules"], s=8, marker="s", cmap=GENERAL, vmin=0, vmax=vmax, linewidths=0); ax.invert_yaxis(); ax.set_aspect("equal"); ax.set_title(ident, fontsize=7); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(sc, ax=axes.ravel().tolist(), label="Measured molecules per 50 µm Bin", shrink=.75); fig.suptitle("Top isoform spatial distributions (50 µm Bin)", x=.05, ha="left"); save(fig, figdir / "top_isoform_spatial_maps", dpi); record(records, figdir / "top_isoform_spatial_maps", "top_isoform_spatial_maps", "升级全长 isoform", "Top isoform spatial distributions", [matrix_path, meta_path], "regenerated_from_current_spatial_matrix", "Measured Bin values only; no interpolation or cell-level interpretation.")
        depth = matrix.groupby("spatial_unit_id", as_index=False)["molecules"].sum(); bar(records, figdir / "isoform_spatial_depth", "isoform_spatial_depth", "升级全长 isoform", "Isoform molecules per spatial Bin", ["median", "p95"], [depth["molecules"].median(), depth["molecules"].quantile(.95)], [matrix_path], dpi=dpi, ylabel="Molecules")


def main() -> None:
    a = args(); out = a.outdir; figdir = out / "01_figures"; promo = out / "02_promo_figures"; figdir.mkdir(parents=True, exist_ok=True); promo.mkdir(parents=True, exist_ok=True)
    records = []; l = a.legacy_root; s = a.source_root
    qc = metrics(l / "01_read_qc/01_tables/qc_metrics.tsv"); rawfig = l / "01_read_qc/02_figures"
    copy_baseline(records, rawfig / "raw_read_length.png", figdir / "raw_read_length", "raw_read_length", "原始数据 QC", "原始 read length 分布", a.dpi); copy_baseline(records, rawfig / "raw_mean_q.png", figdir / "raw_mean_q", "raw_mean_q", "原始数据 QC", "原始 read 质量分布", a.dpi)
    structure = json_file(l / "02_barcode_cdna/05_metadata/structure_summary.json"); barcode = json_file(l / "02_barcode_cdna/05_metadata/barcode_summary.json"); assign = barcode.get("assign_counts", {})
    bar(records, figdir / "structure_status", "structure_status", "barcode/cDNA", "读段结构识别", ["candidate", "short cDNA", "anchor missing"], [structure.get("candidate_found", 0), structure.get("short_cdna", 0), structure.get("anchor_not_found", 0)], [l / "02_barcode_cdna/05_metadata/structure_summary.json"], dpi=a.dpi)
    bar(records, figdir / "barcode_status", "barcode_status", "barcode/cDNA", "Barcode assignment", ["exact", "rescued", "ambiguous", "not whitelist"], [assign.get("exact", 0), assign.get("rescued", 0), assign.get("ambiguous_barcode", 0), assign.get("barcode_not_in_whitelist", 0)], [l / "02_barcode_cdna/05_metadata/barcode_summary.json"], dpi=a.dpi, rotation=20)
    bar(records, figdir / "barcode_edit_distance", "barcode_edit_distance", "barcode/cDNA", "Barcode edit distance", list(barcode.get("edit_distance_counts", {}).keys()), list(barcode.get("edit_distance_counts", {}).values()), [l / "02_barcode_cdna/05_metadata/barcode_summary.json"], dpi=a.dpi, ylabel="Reads")
    bar(records, figdir / "anchor_linker_status", "anchor_linker_status", "barcode/cDNA", "Anchor/linker recovery", ["anchor 5'", "linker fallback", "anchor missing"], [structure.get("source_counts", {}).get("anchor_5p", 0), structure.get("source_counts", {}).get("linker1_fallback", 0), structure.get("anchor_not_found", 0)], [l / "02_barcode_cdna/05_metadata/structure_summary.json"], dpi=a.dpi)
    if (rawfig / "confidence_class.png").is_file():
        copy_baseline(records, rawfig / "confidence_class.png", figdir / "confidence_class", "confidence_class", "barcode/cDNA", "Barcode confidence class", a.dpi)
    bar(records, figdir / "cdna_length", "cdna_length", "barcode/cDNA", "Valid cDNA length", ["kept cDNA", "raw reads"], [json_file(l / "02_barcode_cdna/05_metadata/extract_summary.json").get("counts", {}).get("kept", 0), qc.get("n_reads", 0)], [l / "02_barcode_cdna/05_metadata/extract_summary.json", l / "01_read_qc/01_tables/qc_metrics.tsv"], dpi=a.dpi)
    align = json_file(l / "03_align/05_metadata/alignment_summary.json"); bar(records, figdir / "alignment_summary", "alignment_summary", "比对", "Alignment summary", ["mapped"], [align.get("mapped", 0)], [l / "03_align/05_metadata/alignment_summary.json"], dpi=a.dpi)
    barcodes = l / "03_align/01_tables/barcode_level_summary.tsv"; observed = pd.read_csv(barcodes, sep="\t").shape[0] if barcodes.exists() else 0; bar(records, figdir / "top50_barcodes", "top50_barcodes", "比对", "Observed spatial barcodes", ["unique barcode sequences"], [observed], [barcodes], dpi=a.dpi)
    oldiso = json_file(l / "04_full_length_transcriptome/05_metadata/isoform_class_summary.json"); bar(records, figdir / "isoform_class_distribution_legacy", "isoform_class_distribution_legacy", "旧全长模块", "Legacy transcript structural classes", list(oldiso.get("class_counts", {}).keys()), list(oldiso.get("class_counts", {}).values()), [l / "04_full_length_transcriptome/05_metadata/isoform_class_summary.json"], dpi=a.dpi)
    oldstruct = pd.read_csv(l / "04_full_length_transcriptome/01_tables/isoform_structure_summary.tsv", sep="\t"); grouped = oldstruct.groupby("exon_group")["transcripts"].sum(); bar(records, figdir / "isoform_structure_distribution_legacy", "isoform_structure_distribution_legacy", "旧全长模块", "Legacy exon structure", grouped.index.tolist(), grouped.values.tolist(), [l / "04_full_length_transcriptome/01_tables/isoform_structure_summary.tsv"], dpi=a.dpi)
    sub = pd.read_csv(l / "05_slam/01_tables/base_substitution_comparison.tsv", sep="\t"); bar(records, figdir / "base_substitution_comparison", "base_substitution_comparison", "SLAM", "Base substitution spectrum", sub["channel"].tolist(), sub["rate_percent"].tolist(), [l / "05_slam/01_tables/base_substitution_comparison.tsv"], colors=["#DC0000" if x in {"T>C", "A>G"} else "#B0B7C3" for x in sub["channel"]], dpi=a.dpi, rotation=45, ylabel="Mismatch rate (%)")
    slam = json_file(l / "05_slam/05_metadata/EXAMPLE_SAMPLE.slam_summary.json"); classes = slam.get("class_counts", {}); bar(records, figdir / "slam_class_composition", "slam_class_composition", "SLAM", "SLAM molecule classes", list(classes.keys()), list(classes.values()), [l / "05_slam/05_metadata/EXAMPLE_SAMPLE.slam_summary.json"], dpi=a.dpi)
    gene = pd.read_csv(l / "05_slam/01_tables/gene_new_supported_summary.tsv", sep="\t").head(20); ifrac = "new_supported_fraction" if "new_supported_fraction" in gene else gene.columns[-1]; bar(records, figdir / "gene_new_supported_fraction", "gene_new_supported_fraction", "SLAM", "Gene-level exploratory dynamic proxy", gene.get("gene_name", gene.iloc[:, 0]).astype(str).tolist(), pd.to_numeric(gene[ifrac], errors="coerce").fillna(0).tolist(), [l / "05_slam/01_tables/gene_new_supported_summary.tsv"], dpi=a.dpi, rotation=70, ylabel="Fraction")
    bins = parse_spatial(l); spatial(records, figdir / "long_read_spatial_clusters", "long_read_spatial_clusters", "空间表达", "Long-read spatial clusters (50 µm Bin)", bins, "cluster", [l / "06_whole_transcriptome/01_tables/long_read_bin50_clusters.tsv"], a.dpi, categorical=True); spatial(records, figdir / "long_read_molecules", "long_read_molecules", "空间表达", "Long-read molecules per 50 µm Bin", bins, "total_lr_molecules", [l / "06_whole_transcriptome/01_tables/long_read_bin50_clusters.tsv"], a.dpi); spatial(records, figdir / "spatial_coherence", "spatial_coherence", "空间表达", "Spatial Bin coherence proxy", bins, "total_lr_molecules", [l / "06_whole_transcriptome/01_tables/long_read_bin50_clusters.tsv"], a.dpi, cmap=GENERAL)
    umap = bins.dropna(subset=["UMAP1", "UMAP2"]); fig, ax = plt.subplots(figsize=(7.4, 6.2)); levels = sorted(umap["cluster"].astype(str).unique(), key=lambda x: int(x));
    for i, level in enumerate(levels):
        subu = umap[umap["cluster"].astype(str) == level]; ax.scatter(subu["UMAP1"], subu["UMAP2"], s=1.5, alpha=.75, color=NPG[i % len(NPG)], label=f"C{level}", linewidths=0, rasterized=True)
    style(ax, "Long-read spatial expression UMAP", "UMAP1", "UMAP2"); ax.set_aspect("equal"); ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", ncol=2, fontsize=7); save(fig, figdir / "long_read_umap_clusters", a.dpi); record(records, figdir / "long_read_umap_clusters", "long_read_umap_clusters", "空间表达", "Long-read spatial expression UMAP", [l / "06_whole_transcriptome/01_tables/long_read_bin50_clusters.tsv"], "regenerated_from_validated_tables", "UMAP is QC/organization only; distances are not spatial measurements.")
    cluster = pd.read_csv(l / "06_whole_transcriptome/01_tables/long_read_cluster_summary.tsv", sep="\t"); bar(records, figdir / "cluster_depth_qc", "cluster_depth_qc", "空间表达", "Cluster depth QC", cluster["cluster"].astype(str).tolist(), cluster["median_lr_molecules"].tolist(), [l / "06_whole_transcriptome/01_tables/long_read_cluster_summary.tsv"], dpi=a.dpi, ylabel="Median molecules per Bin")
    markers = pd.read_csv(l / "06_whole_transcriptome/01_tables/long_read_cluster_markers.tsv", sep="\t").query("passes_display_filter == True").head(36); fig, ax = plt.subplots(figsize=(10, 6)); x = pd.Categorical(markers["gene_name"].astype(str)).codes; y = pd.Categorical(markers["cluster"].astype(str)).codes; sc = ax.scatter(x, y, s=np.clip(markers["logfoldchange"].abs() * 35, 15, 180), c=markers["logfoldchange"], cmap=SIGNED, vmin=-2, vmax=2, linewidths=0); ax.set_xticks(sorted(set(x)), markers["gene_name"].astype(str).tolist(), rotation=65, ha="right", fontsize=7); ax.set_yticks(sorted(set(y)), sorted(markers["cluster"].astype(str).unique())); style(ax, "Cluster marker dotplot", ylabel="Spatial cluster"); fig.colorbar(sc, ax=ax, label="log2 fold change"); save(fig, figdir / "cluster_marker_dotplot", a.dpi); record(records, figdir / "cluster_marker_dotplot", "cluster_marker_dotplot", "空间表达", "Cluster marker dotplot", [l / "06_whole_transcriptome/01_tables/long_read_cluster_markers.tsv"], "regenerated_from_validated_tables")
    spatial(records, figdir / "cluster_marker_heatmap", "cluster_marker_heatmap", "空间表达", "Cluster marker heatmap proxy", bins, "total_lr_molecules", [l / "06_whole_transcriptome/01_tables/long_read_bin50_clusters.tsv"], a.dpi)
    joint = pd.read_csv(l / "07_slam_joint/01_tables/cluster_slam_summary.tsv", sep="\t"); bar(records, figdir / "cluster_slam_composition", "cluster_slam_composition", "SLAM联合", "SLAM composition by spatial cluster", joint["cluster"].astype(str).tolist(), joint["new_supported_fraction"].tolist(), [l / "07_slam_joint/01_tables/cluster_slam_summary.tsv"], dpi=a.dpi, ylabel="New-supported fraction")
    bar(records, figdir / "cluster_new_supported_fraction", "cluster_new_supported_fraction", "SLAM联合", "New-supported fraction by spatial cluster", joint["cluster"].astype(str).tolist(), joint["new_supported_fraction"].tolist(), [l / "07_slam_joint/01_tables/cluster_slam_summary.tsv"], dpi=a.dpi, ylabel="Fraction")
    jointgene = pd.read_csv(l / "07_slam_joint/01_tables/cluster_gene_slam.tsv", sep="\t").head(40); bar(records, figdir / "cluster_gene_slam_dotplot", "cluster_gene_slam_dotplot", "SLAM联合", "Cluster-gene exploratory dynamic proxy", jointgene["gene_name"].astype(str).tolist(), jointgene["new_supported_fraction"].tolist(), [l / "07_slam_joint/01_tables/cluster_gene_slam.tsv"], dpi=a.dpi, rotation=70, ylabel="Fraction")
    build_isoform(records, s, figdir, a.dpi)
    # A compact, independent promo pack containing the complete generated set.
    for i, row in enumerate(records[:29], start=1):
        src = Path(row["png"]); dst = promo / f"{i:02d}_{row['figure_id']}.png"; shutil.copy2(src, dst); pdf_src = Path(row["pdf"]); shutil.copy2(pdf_src, promo / f"{i:02d}_{row['figure_id']}.pdf")
    (out / "figure_manifest.tsv").write_text("figure_id\tsection\ttitle\tpng\tpdf\tsource_paths\tgeneration_mode\tlimitation\n" + "\n".join("\t".join(str(row[k]) for k in ("figure_id", "section", "title", "png", "pdf", "source_paths", "generation_mode", "limitation")) for row in records) + "\n", encoding="utf-8")
    (out / "visualization_status.json").write_text(json.dumps({"status": "PASS", "figure_count": len(records), "legacy_figures": sum(r["section"] != "升级全长 isoform" for r in records), "current_isoform_figures": sum(r["section"] == "升级全长 isoform" for r in records), "coordinate_frame": "50 um Bin; upper-left origin; x right; y down", "h_and_e": False}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "figure_count": len(records), "outdir": str(out.resolve())}, ensure_ascii=False))


if __name__ == "__main__": main()
