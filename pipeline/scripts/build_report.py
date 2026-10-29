from __future__ import annotations

import argparse
import base64
import csv
import html
import json
from pathlib import Path


def rows(path: str, limit: int = 20):
    with open(path, encoding="utf-8", errors="replace") as handle:
        data = list(csv.DictReader(handle, delimiter="\t"))
    return data[:limit], len(data)


def data_uri(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        return ""
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def figure_cards(manifest: Path, section: str | None = None) -> tuple[str, list[dict[str, str]]]:
    if not manifest.is_file():
        return '<p class="empty">当前没有可展示图形。</p>', []
    records = list(csv.DictReader(manifest.open(encoding="utf-8"), delimiter="\t"))
    cards = []
    selected = [item for item in records if section is None or item.get("section") == section]
    for index, item in enumerate(selected, 1):
        uri = data_uri(Path(item.get("png", "")))
        if not uri:
            continue
        title = html.escape(item.get("title", item.get("figure_id", "figure")))
        limitation = html.escape(item.get("limitation", ""))
        cards.append(f'<article class="figure-card"><button class="expand figure-expand" type="button" data-title="{title}"><img src="{uri}" alt="{title}"></button><figcaption><strong>图 {index} · {title}</strong><span>{limitation}</span></figcaption></article>')
    return '<div class="figure-grid">' + "".join(cards) + "</div>", selected


def table(data, total):
    if not data:
        return '<p class="empty">无可展示记录。</p>'
    cols = list(data[0]); head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(r.get(c, '')))}</td>" for c in cols) + "</tr>" for r in data)
    return f'<p>共 {total:,} 条，当前展示前 {len(data):,} 条。</p><div class="table-scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def main() -> None:
    ap = argparse.ArgumentParser()
    for name in ("config", "outdir", "qc", "isoforms", "novel", "matrices", "events", "svi"):
        ap.add_argument("--" + name, required=True)
    ap.add_argument("--figures", action="append", required=False)
    ap.add_argument("--run-metadata", required=False)
    a = ap.parse_args(); out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    qc, qc_n = rows(a.qc); iso, iso_n = rows(a.isoforms); novel, novel_n = rows(a.novel)
    events = json.loads(Path(a.events).read_text(encoding="utf-8")); svi = json.loads(Path(a.svi).read_text(encoding="utf-8"))
    if a.figures:
        registry = []
        for manifest_path in a.figures:
            _, current = figure_cards(Path(manifest_path))
            registry.extend(current)
        def cards_for_sections(sections):
            chunks = []
            for manifest_path in a.figures:
                for section in sections:
                    chunks.append(figure_cards(Path(manifest_path), section)[0])
            return "".join(chunks)
        figures = cards_for_sections([None])
        legacy_figures = cards_for_sections(["原始数据 QC", "barcode/cDNA", "比对", "旧全长模块", "SLAM", "空间表达", "SLAM联合"])
        isoform_figures = cards_for_sections(["升级全长 isoform", "Spl-IsoQuant", "Spl-IsoFind"])
    else:
        figures = legacy_figures = isoform_figures = '<p class="empty">未提供可视化清单。</p>'
        registry = []
    def card(label, value): return f'<div class="card"><span>{html.escape(label)}</span><strong>{html.escape(str(value))}</strong></div>'
    sections = [
        ("overview", "1 项目概览", '<div class="cards">' + card("样本", "EXAMPLE_SAMPLE") + card("物种", "Mus musculus") + card("平台", "ONT / BBV2.4") + card("空间单位", "50 µm Bin") + '</div><div class="notice">本报告为单样本、50 µm Bin 分析；没有细胞分割、细胞类型注释或生物学重复。空间图显示实测 Bin 值，不代表 single-cell 分辨率。</div>'),
        ("qc", "2 Long-read QC / 原始数据与旧模块", table(qc, qc_n)),
        ("legacy", "3 已验证旧模块可视化", '<p>下列图形依据旧项目 EXAMPLE_SAMPLE 已验证结果表和元数据重绘；原始逐 read 图形作为已验证基线资产保留。</p>' + legacy_figures),
        ("isoform", "4 Transcript and Isoform / 本次升级的全长 isoform", f'<p>候选 isoform 记录：{iso_n:,}；高置信 novel 记录：{novel_n:,}。</p>{table(iso, iso_n)}<h3>高置信 novel isoform</h3>{table(novel, novel_n)}'),
        ("isoform_figures", "5 Isoform 专属可视化", isoform_figures),
        ("matrix", "6 Spatial matrices", f'<p>isoform × spatial-unit 矩阵：<code>{html.escape(a.matrices)}</code>。</p>'),
        ("events", "7 Splicing、APA、TSS、DTU 和 switch", f'<pre>{html.escape(json.dumps(events, indent=2, ensure_ascii=False))}</pre>'),
        ("svi", "8 Spatially Variable Isoforms / 空间可变 isoform", f'<pre>{html.escape(json.dumps(svi, indent=2, ensure_ascii=False))}</pre>'),
        ("methods", "9 Methods and Limitations / 方法、限制与可追溯性", '<p>barcode 解析采用 BBV2.4 三段 4 bp barcode、两个 linker 和 10 bp UMI；全长 isoform 主调用采用 Spl-IsoQuant 4.0.0；空间单位为 50 µm Bin。</p><p>当前无 H&E 或可验证的图像配准，因此空间图使用白底坐标图，不进行组织图像叠加。当前无细胞分割和细胞类型标签，不输出 cell-level 或 cell-type-specific 结论。单样本无重复，所有 SLAM/dynamic proxy 仅作描述性结果。</p>'),
    ]
    nav = "".join(f'<a href="#{sid}">{html.escape(title)}</a>' for sid, title, _ in sections)
    content = "".join(f'<section id="{sid}" class="report-section"><h2>{html.escape(title)}</h2>{body}</section>' for sid, title, body in sections)
    page = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>EXAMPLE_SAMPLE Spatial Full-Length Transcriptome</title><style>:root{{--navy:#253746;--cyan:#4DBBD5;--ink:#27323a;--muted:#667681;--line:#dce3e7;--soft:#f5f8f9}}*{{box-sizing:border-box}}body{{margin:0;background:#f2f5f6;color:var(--ink);font:13px Arial,"Noto Sans SC",sans-serif;line-height:1.55}}.sidebar{{position:fixed;inset:0 auto 0 0;width:230px;background:var(--navy);color:#fff;padding:24px 18px;overflow:auto}}.sidebar h1{{font-size:17px}}.sidebar a{{display:block;color:#dce7eb;text-decoration:none;padding:7px;border-left:2px solid transparent}}.sidebar a:hover{{border-left-color:var(--cyan);background:#344b5d}}main{{margin-left:230px;max-width:1400px;padding:26px 32px 60px}}.title,section{{background:#fff;border:1px solid var(--line);border-radius:7px;padding:21px 24px;margin-bottom:16px}}.title{{border-top:4px solid var(--navy)}}h1{{margin:0;color:var(--navy)}}h2{{color:var(--navy);border-bottom:1px solid var(--line);padding-bottom:7px}}h3{{color:#315468}}.cards{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:9px}}.card{{background:var(--soft);border-top:3px solid var(--cyan);padding:10px}}.card span,.card strong{{display:block}}.card span{{color:var(--muted);font-size:10px}}.card strong{{color:var(--navy);font-size:18px}}.notice{{padding:11px;background:#fff8e5;border-left:3px solid #e4aa38;margin:12px 0}}.figure-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;margin:14px 0}}.figure-card{{border:1px solid var(--line);margin:0;background:#fff}}.expand{{display:block;width:100%;border:0;background:#fff;padding:9px;cursor:zoom-in}}.expand img{{display:block;max-width:100%;max-height:610px;margin:auto}}figcaption{{display:flex;flex-direction:column;padding:8px 10px;border-top:1px solid var(--line);background:#fafbfb}}figcaption span{{color:var(--muted);font-size:10px;margin-top:3px}}.table-scroll{{overflow:auto;max-height:500px}}table{{border-collapse:collapse;width:100%;font-size:11px}}th,td{{padding:6px 8px;border-bottom:1px solid #e2eaed;text-align:left;white-space:nowrap}}th{{background:#edf4f5;position:sticky;top:0}}pre{{white-space:pre-wrap;background:var(--soft);padding:12px}}code{{overflow-wrap:anywhere}}.modal{{display:none;position:fixed;inset:0;background:#101820e8;z-index:5;padding:25px;align-items:center;justify-content:center}}.modal.open{{display:flex}}.modal-box{{position:relative;background:#fff;padding:38px 14px 14px;max-width:96vw;max-height:94vh}}.modal-box img{{max-width:92vw;max-height:82vh}}.modal-close{{position:absolute;right:8px;top:3px;border:0;background:none;font-size:26px}}.footer{{color:var(--muted);font-size:10px}}@media(max-width:900px){{.sidebar{{position:static;width:auto}}main{{margin:0;padding:10px}}.cards,.figure-grid{{grid-template-columns:1fr}}}}@media print{{.sidebar,.modal{{display:none!important}}main{{margin:0;padding:0;width:100%}}section{{break-inside:avoid;box-shadow:none}}.figure-card{{break-inside:avoid}}}} </style></head><body><aside class="sidebar"><h1>EXAMPLE_SAMPLE</h1><p>Spatial Full-Length Transcriptome</p><nav class="sidebar">{nav}</nav></aside><main><div class="title"><h1>EXAMPLE_SAMPLE 小鼠三代空间全长转录组分析</h1><p>ONT · BBV2.4 · 50 µm Bin · 离线自包含</p></div>{content}<div class="footer">Generated by pipeline/scripts/build_report.py; embedded figure records: {len(registry)}</div></main><div id="figure-modal" class="modal"><div class="modal-box"><button class="modal-close" type="button">×</button><div class="modal-title"></div><img alt=""></div></div><script>const modal=document.getElementById('figure-modal'),big=modal.querySelector('img'),title=modal.querySelector('.modal-title');document.querySelectorAll('.expand').forEach(b=>b.onclick=()=>{{const i=b.querySelector('img');big.src=i.src;big.alt=i.alt;title.textContent=b.dataset.title;modal.classList.add('open')}});function closeModal(){{modal.classList.remove('open');big.src=''}}modal.querySelector('.modal-close').onclick=closeModal;modal.onclick=e=>{{if(e.target===modal)closeModal()}};document.onkeydown=e=>{{if(e.key==='Escape')closeModal()}};</script></body></html>'''
    report = out / "EXAMPLE_SAMPLE_report.html"; report.write_text(page, encoding="utf-8")
    if a.run_metadata:
        metadata = {"status": "PASS", "sample_id": "EXAMPLE_SAMPLE", "report": str(report.resolve()), "config": str(Path(a.config).resolve()), "figure_manifests": [str(Path(x).resolve()) for x in (a.figures or [])]}
        Path(a.run_metadata).parent.mkdir(parents=True, exist_ok=True)
        Path(a.run_metadata).write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    manifest = out / "delivery_manifest.tsv"
    with manifest.open("w", encoding="utf-8") as handle:
        handle.write("description\tsource\ttarget\tstatus\n")
        if a.figures:
            for item in registry:
                target = Path(item["png"]); handle.write(f"figure {item['figure_id']}\t{target}\t{target}\tup_to_date\n")
        handle.write(f"EXAMPLE_SAMPLE final offline HTML report\t{report.resolve()}\t{report.resolve()}\tup_to_date\n")


if __name__ == "__main__": main()
