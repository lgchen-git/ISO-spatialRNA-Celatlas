from __future__ import annotations
import argparse, csv, re
from pathlib import Path

def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--report", required=True); ap.add_argument("--manifest", required=True); a = ap.parse_args()
    report = Path(a.report); errors = []; text = report.read_text(encoding="utf-8", errors="replace") if report.is_file() else ""
    if not text: errors.append("report missing or empty")
    if '<aside class="sidebar">' not in text or 'class="footer"' not in text: errors.append("navigation/footer missing")
    if re.search(r'<(?:img|script)[^>]+(?:src|href)=["\'](?!data:|#)', text, re.I): errors.append("external resource found")
    for section in ("项目概览", "Long-read QC", "Transcript and Isoform", "Spatially Variable Isoforms", "Methods and Limitations"):
        if section not in text: errors.append("missing section: " + section)
    manifest = Path(a.manifest); rows = list(csv.DictReader(manifest.open(encoding="utf-8"), delimiter="\t")) if manifest.is_file() else []
    matches = [r for r in rows if Path(r.get("target", "")).resolve() == report.resolve()]
    if len(matches) != 1: errors.append("manifest must contain exactly one final report")
    if not rows or rows[-1] not in matches: errors.append("final report row must be last")
    if errors: print("ok: false"); print("\n".join("ERROR: " + e for e in errors)); return 1
    print("ok: true"); return 0

if __name__ == "__main__": raise SystemExit(main())
