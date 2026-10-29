from pathlib import Path
import subprocess, tempfile

def test_report():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        (root / "qc.tsv").write_text("metric\tvalue\nn_reads\t1\n")
        (root / "iso.tsv").write_text("gene\tisoform\nG\tI\n")
        (root / "novel.tsv").write_text("gene\tisoform\n")
        (root / "events.json").write_text('{"AS":"BLOCKED"}')
        (root / "svi.json").write_text('{"status":"BLOCKED"}')
        subprocess.run(["python3", "pipeline/scripts/build_report.py", "--config", "config/run.example.yaml", "--outdir", str(root / "out"), "--qc", str(root / "qc.tsv"), "--isoforms", str(root / "iso.tsv"), "--novel", str(root / "novel.tsv"), "--matrices", str(root / "iso.tsv"), "--events", str(root / "events.json"), "--svi", str(root / "svi.json")], check=True)
        subprocess.run(["python3", "pipeline/validation/validate_report.py", "--report", str(root / "out/EXAMPLE_SAMPLE_report.html"), "--manifest", str(root / "out/delivery_manifest.tsv")], check=True)

if __name__ == "__main__": test_report(); print("report test: OK")
