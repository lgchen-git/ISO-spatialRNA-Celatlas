from __future__ import annotations

import argparse
from pathlib import Path

from common import config, dump_json


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--config", required=True); ap.add_argument("--reads", required=True); ap.add_argument("--isoforms", required=True); ap.add_argument("--metadata", required=True); ap.add_argument("--outdir", required=True); ap.add_argument("--status", required=True)
    a = ap.parse_args(); c = config(a.config); out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    try:
        import pandas as pd
    except ModuleNotFoundError as exc:
        raise SystemExit("isoform event analysis requires pandas") from exc
    reads = pd.read_parquet(a.reads)
    if reads.empty: raise SystemExit("No reads for event analysis")
    status = {"AS": "BLOCKED", "APA": "EXPLORATORY_PENDING_INTERNAL_PRIMING_QC", "TSS": "EXPLORATORY_PENDING_5PRIME_COMPLETENESS", "DTU": "OBSERVATIONAL_BIN_LEVEL_ONLY", "isoform_switch": "OBSERVATIONAL_BIN_LEVEL_ONLY", "reason": "Event caller requires transcript assignments and coverage filters; no unsupported result is emitted"}
    dump_json(a.status, status)
    (out / "event_status.tsv").write_text("module\tstatus\n" + "\n".join(f"{k}\t{v}" for k, v in status.items() if k != "reason") + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
