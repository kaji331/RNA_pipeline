# RNA pipeline

This pipeline contains three tool-chains: STAR-RSEM/HISAT2-featureCounts/kallisto. It can automatically get read counts from fastq files (paired-end) with qc report.

## Installation

```bash
git clone https://github.com/kaji331/RNA_pipeline.git
```

### Required

Before run this pipeline, you need to do something else:

- Create a conda environment named rseqc with RSeQC installed
- Create a conda environment named multiqc with MultiQC installed
- Install STAR/RSEM, HISAT2/featureCounts(Subread) or kallisto, add them to environment path

## Usage

Before running, you should modify the config_parameters.txt:

- PARALLEL_JOB is the number of simultaneously running samples
- TRIMMOMATIC is the path of trimmomatic jar
- STAR_REF is the path of STAR index made by genomic fasta
- RSEM_REF is the path of RSEM index made by genomic fasta
- HISAT2_REF is the path of HISAT2 index made by genomic fasta
- GTF is the path of GTF file for featureCounts
- KALLISTO_REF is the path of kallisto index made by transcript fasta
- BED is the path of bed12 file for RSeQC
- SHELL is the path of shell, if you have different shell or path (default:/bin/bash)

> Sometimes, you'd better use RSEM to read a gff file only with transcript coding and it will convert gff to gtf for featureCounts. You can convert gtf to bed file with any tools following BED12 rules.

Then, you should make a list of R1 reads (R2 reads list are optional). If you need to combine some fastq files, you also need a list file. See example/ folder.

1. python3 Step0_combine.py (optional)
2. python3 Step1_qc.py
3. Read fastq qc report
4. python3 Step2_trimming.py (optional)
5. python3 Step3_quant.py