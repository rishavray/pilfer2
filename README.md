# PILFER v2
An updated snakemake version of PILFER algorithm

This version uses ShortStack algorithm for alignment and BEDTOOLS for quantification followed by PILFER run.

Clone this repository, and add a fastq folder in the resources folder containing your fastq files.

Then update any other configuration parameters including the location of the bowtie index.

Run the pipeline

```
snakemake --cores 1
```


