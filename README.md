# ISO-Celatlas Spatial Full-Length Transcriptome

Snakemake workflow for ONT BBV2.4 mouse spatial full-length transcriptomics.
The canonical spatial unit is a 50 µm bin. The pipeline preserves the
reference project's custom barcode decoder and FLAMES/gffcompare route, then
adds molecule normalization, sparse isoform matrices, isoform events, SVI
gating, all prior module visualizations, Spl-IsoQuant/Spl-IsoFind spatial
figures, and one offline Chinese report.

```bash
make setup-layout
make lint
make validate
make dry-run
make run CORES=1
```

Raw data and references remain external and read-only. The project `outputs/`
path points to a new versioned data-side result directory; the established
results are never overwritten by this workflow.
