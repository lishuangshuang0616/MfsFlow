# MhsFlow

**MhsFlow** is a comprehensive and flexible pipeline designed for processing high-sensitivity full-length transcriptome data, automated to handle workflows from raw FASTQ data to gene expression quantification and optimized for MGI sequencing chemistry.

## Features

- **End-to-End Pipeline**: Handles filtering, mapping, counting, and statistical analysis.
- **Simple CLI Configuration**: Customer-facing parameters are passed on the command line; a run YAML is generated automatically for provenance.
- **Multi-Species Support**: Optimized for Human and Mouse genomes.
- **Efficient Processing**: Distinct stages for Filtering, Mapping, Counting, and Summarising.
- **Detailed Reporting**: Generates comprehensive statistics and analysis results.

## Prerequisites

Ensure you have the following installed on your system:

- **Python 3.8+**
- **Python packages**: Install via pip:
  ```bash
  pip install -r requirements.txt
  ```
- **STAR**: For genome alignment.
- **Samtools**: For BAM file manipulation.
- **Pigz**: For parallel gzip compression.
- **Seqkit**: For FASTQ manipulation.
- **FeatureCounts**: For read counting.

## Installation

1.  **Clone the Repository**
    ```bash
    git clone https://github.com/Yzh25/Mhsflt_toolkit.git
    cd Mhsflt_toolkit
    ```

2.  **Environment Setup**
    The toolkit relies on specific directory structures and environment settings. Run the release script to initialize the environment:
    ```bash
    sh release_env.sh
    ```

## Reference

1.  **Build Reference Index**
    Download the genome sequence (FASTA) and gene annotation (GTF) files from [GENCODE](https://www.gencodegenes.org). Build the genome index using STAR.

## Usage

Run the analysis pipeline with command-line parameters:

```bash
python3 run_analysis_pipeline.py \
  --fastqs /path/to/fastq_dir \
  --genomeDir /path/to/reference \
  --sample SAMPLE_NAME \
  --plate 1 \
  --outdir /path/to/output \
  --threads 20
```

By default, the pipeline expects STAR index at `<genomeDir>/star` and GTF at `<genomeDir>/genes/genes.gtf`.
The pipeline generates `XPRESS_PROCESSING/config/run_config.yaml` automatically for reproducibility.

Input modes:

- If R2 is 20 bp longer than R1, the last 20 bp of R2 are treated as the read barcode.
- If R1 and R2 have the same length, provide `--samplesheet`; barcode is assigned per FASTQ pair.
- Other read length structures are rejected.

Samplesheet format for equal-length R1/R2 data:

```csv
read1,read2,barcode
A_R1.fq.gz,A_R2.fq.gz,ACGTACGTACGTACGTACGT
B_R1.fq.gz,B_R2.fq.gz,TGCATGCATGCATGCATGCA
```

The barcode sequence is looked up in the expected barcode table generated from
`--plate`, `--manual`, or `--expectBarcode`. The pipeline derives the well ID
and whether the FASTQ pair is UMI or internal from that lookup.

Barcode modes:

```bash
--plate 1                  # automatic plate barcode list
--manual 20,21             # manual barcode IDs
--expectBarcode barcodes.tsv
--discoverBarcodes         # infer plate/manual barcode set from observed reads
```

For `--manual`, `--plate`, and `--expectBarcode`, barcode detection is strict:
if none of the expected barcodes are observed, the run stops instead of falling
back to top-count barcodes. `--discoverBarcodes` compares observed read barcodes
against the bundled manual/plate barcode lists with exact/Hamming-1 matching,
writes a discovery report, then continues with the inferred barcode set.

## Output Structure

Results are saved under your specified `out_dir`:

- `XPRESS_PROCESSING/config/`: Generated run configuration and expected barcode tables.
- `XPRESS_PROCESSING/logs/`: Pipeline logs.
- `XPRESS_PROCESSING/barcodes/`: Kept barcodes, barcode binning, and barcode discovery report.
- `XPRESS_PROCESSING/expression/`: MEX matrices and `<project>.h5ad`.
- `XPRESS_PROCESSING/stats/`: QC tables, Q30 statistics, and plots.
- `XPRESS_PROCESSING/<project>.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam`: Final UB-corrected sorted BAM.
- `XPRESS_PROCESSING/intermediate/`: Temporary BAM and split FASTQ working files.
- `outs/`: Final reports and customer-facing outputs.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
