# Runtime environment

The workflow is designed for a conda environment created from
`pipeline/envs/iso-spatial-full-length.yml`:

```bash
micromamba create -n iso-spatial-full-length -f pipeline/envs/iso-spatial-full-length.yml
micromamba activate iso-spatial-full-length
```

The required command-line layer includes Snakemake, minimap2, samtools,
bedtools, gffread, and IsoQuant. Python dependencies include pysam, pandas,
NumPy/SciPy, Scanpy/anndata, pyarrow, and statsmodels. Spl-IsoFind is pinned
in the environment file. Spl-IsoQuant/FLAMES are method adapters and must be
installed or exposed through their validated source entrypoints before the
corresponding modules are enabled.

The canonical entrypoint is `pipeline/run_snakemake.sh`, which uses the
configured environment and places runtime cache/temp files under `tmp/`.
Spl-IsoQuant remains an explicit executable configured in the run YAML.
