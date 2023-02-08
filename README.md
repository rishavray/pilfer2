# PILFER v2
An updated snakemake version of PILFER algorithm

This version uses [ShortStack](https://github.com/MikeAxtell/ShortStack) algorithm for alignment and [BEDTOOLS](https://github.com/arq5x/bedtools2) for quantification followed by PILFER run. The piRNAs are matched against the [piRBase](http://bigdata.ibp.ac.cn/piRBase/index.php) curated gold standard piRNA set.

PILFER algorith is also now modularized and uses pandas for reading the bedfiles. It then calculates the CPM to call the clusters. Additionally the raw count table is aggregated.

Clone this repository, and add a fastq folder in the resources folder containing your fastq files.

Then update any other configuration parameters including the location of the bowtie index.

Run the pipeline

```
snakemake --cores 1
```


