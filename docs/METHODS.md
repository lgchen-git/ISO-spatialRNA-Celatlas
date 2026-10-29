# Methods

Input QC streams the raw ONT FASTQ. Barcode parsing uses the BBV2.4
`C4L15C4L15C4U10T18` structure, whitelist correction up to one mismatch, and
the supplied chip coordinate files. Alignment uses minimap2. FLAMES is the
primary isoform discovery caller and gffcompare classifies known/novel
structures. Matrices are deduplicated by barcode, UMI, and gene and stored as
long sparse tables. The canonical execution path is Snakemake; logs and
benchmarks are written to `log/`, transient runtime material to `tmp/`, and
deliverables to `outputs/the configured sample/<stage_dir>/`.

AS, APA, TSS, DTU, isoform switch, and SVI are only reportable after their
coverage, coordinate, and design gates pass. A single sample has no biological
replicate; bin-level contrasts are observational and are not biological
replication. No cell-level result is generated without validated segmentation.
