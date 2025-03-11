# Metagenomic Taxonomy Pipeline

## Overview

This Nextflow pipeline is designed for metagenomic analysis, specifically for identifying species in a host sample using a combination of de novo assembly, BLAST searches, and taxonomic classification. The pipeline processes paired-end FASTQ files, performs quality control, removes host contamination (if specified), assembles reads into contigs, and identifies species using BLASTn and BLASTx (via DIAMOND). The pipeline outputs tally files, extracted FASTA files, and a heatmap of identified microorganisms.

## Requirements

- Nextflow (version 20.07.1 or later)
- Java 8 or later
- Required software tools (e.g., Trimmomatic, SPAdes, BLAST, DIAMOND, etc.) should be installed and available in your `PATH`.

## Usage
Create a new conda env with environment.yml to get all dependencies

```bash
conda env create -f environment.yml
```

To run the pipeline, use the following command:

```bash
nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR --local_blast_nt BLAST_DATABASE_DIR --local_diamond DIAMOND_DIR
```

### Required Parameters

- `--input INPUT_DIR`: Directory containing paired-end FASTQ files.
- `--output OUTPUT_DIR`: Directory to store pipeline output files. If the directory does not exist, it will be created.
- `--local_blast_nt BLAST_DATABASE_DIR`: Directory containing the local BLAST nt database.
- `--local_diamond DIAMOND_DIR`: Directory containing the local DIAMOND database.

### Optional Parameters

- `--host_fasta HOST_FASTA`: FASTA file of the host genome for host read removal.
- `--host_bt2_index INDEX_DIRECTORY`: Directory containing a pre-built Bowtie2 index for host read removal.
- `--threads INT`: Number of threads to use for parallel processing (default: 10).
- `--ref REFERENCE_FASTA`: Reference genome FASTA file for contig alignment.
- `--minLen INT`: Minimum read length to keep after trimming (default: 75 bp).
- `--minTrimQual INT`: Minimum average base quality for trimming (default: 20).
- `--ncbi_tax_dir`: Directory to store the NCBI taxonomy database for Taxonomizr.
- `--blast_tax_dir`: Directory containing the local NCBI taxonomy database for BLAST.

### Example Command

```bash
nextflow run main.nf --input /path/to/fastq_files --output /path/to/output --local_blast_nt /path/to/blast_nt_db --local_diamond /path/to/diamond_db --host_fasta /path/to/host.fasta --threads 16
```

## Pipeline Steps

1. **Quality Control (QC)**: Initial QC is performed using FASTQC on raw reads.
2. **Trimming**: Adapter and quality trimming are performed using Trimmomatic.
3. **PCR Duplicate Removal**: PCR duplicates are removed using prinseq.
4. **Host Read Removal**: If a host reference or Bowtie2 index is provided, host reads are removed.
5. **De Novo Assembly**: Reads are assembled into contigs using SPAdes.
6. **Contig Alignment**: If a reference genome is provided, contigs are aligned to it using minimap2.
7. **BLASTn Search**: Contigs are searched against the nt database using BLASTn.
8. **BLASTx Search**: Unmatched contigs are searched against the nr database using DIAMOND (BLASTx).
9. **Taxonomic Classification**: Results from BLASTn and BLASTx are tallied and classified using the NCBI taxonomy database.
10. **Heatmap Generation**: A heatmap of identified microorganisms is generated.

## Output Files

- **QC Reports**: FASTQC reports for raw, trimmed, and deduplicated reads.
- **Contigs and Scaffolds**: Assembled contigs and scaffolds in FASTA format.
- **BLAST Results**: Tally files and FASTA files for BLASTn and BLASTx hits.
- **Heatmap**: A heatmap visualizing the abundance of identified microorganisms.

## Notes

- Ensure that the required databases (BLAST nt, DIAMOND nr, and NCBI taxonomy) are correctly set up and accessible.
- If using a large host genome, providing a pre-built Bowtie2 index (`--host_bt2_index`) will significantly speed up the host read removal step.
- The pipeline supports both paired-end and single-end reads. Use the `--single_read` option for single-end reads.

## Help

To display the help message and available options, run:

```bash
nextflow run main.nf --help
```

## License

This pipeline is distributed under the MIT License. See the `LICENSE` file for more details.
