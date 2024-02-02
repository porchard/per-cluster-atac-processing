# Per-cluster processing of snATAC data

This pipeline processed snATAC-data data on a per-cluster pseudobulk basis.

## Dependencies
[Singularity (v. 3)](https://docs.sylabs.io/guides/3.0/user-guide/) and [NextFlow](https://www.nextflow.io/) (>= v. 20.10.0). Containers with the software for each step are pulled from the [Sylabs cloud library](https://cloud.sylabs.io/library) or [Docker hub](https://hub.docker.com/).


## Running

To run the pipeline, you'll need to provide a config.json file like this:

```python
{
    "libraries": {
        "Sample_3172-CV-hg19": {
            "bam": "/path/to/library.bam", # path to pruned library ATAC bam file (from snATAC pipeline)
            "clusters": "/path/to/library.clusters.txt" # two-column TSV file (no header). First column is *RNA* barcode, second column is the cluster assignment for that barcode
        }
    }
}
```

You'll also need to update the `nextflow.config` file in this directory.

Then run the pipeline:

```bin
nextflow run -resume -params-file config.json --genome hg19 --atac_barcodes /path/to/atac-barcode-whitelist.txt.gz --rna_barcodes /path/to/rna-barcode-whitelist.txt.gz --markers Myh1,Myh2,Myh4 --results /path/to/results /path/to/per-cluster-atac-processing/main.nf
```

Where `--genome` is the name of the reference genome to use, and `--markers` is a comma-separated list of marker genes of interest (these genes must be included in the `gene_bed` file).

## Output
* `bam/per-library-pass-qc-nuclei`: Bam files subsetted to pass QC barcodes (per-library)
* `bam/per-library-per-cluster`: Per-library, per-cluster bam files
* `bam/per-cluster`: Per-cluster bam files
* `bam/aggregate`: Aggregate bam file (all clusters and all libraries)
* `peaks/broad`: MACS2 broad peak calling output
* `peaks/narrow`: MACS2 narrow peak calling output (including peak summits)
* `peaks/summit-extension`: Extended summits (default 150 bp either side; overlaps are removed, keeping the one with the highest score, as in paper: 'Epigenomic State Transitions Characterize Tumor Progression in Mouse Lung Adenocarcinoma')
* `bigwig`: Per-cluster bigwig files
* `plot-marker-gene-signal`: ATAC per-cluster marker gene plot