# BRCA-PrimaryRecurrent


Scripts for somatic variant calling and allele-specific copy number profiling in matched tumor/normal whole exome sequencing data.

Required Software:

    Python 3.6+
    vcfpy
    Snakemake
    samtools
    bcftools
    GATK 4
    Strelka2
    VarDictJava
    VarScan2
    R
    Sequenza (R Package)
    Sequenza-utils

## Usage
```
snakemake -s mutect2.snake --configfile config.yaml
snakemake -s strelka2.snake --configfile config.yaml
snakemake -s vardictjava.snake --configfile config.yaml
snakemake -s sequenza.snake --configfile config.yaml
snakemake -s varscan2.snake --configfile config.yaml
```

## Author
Jennifer Shah
<jennshah@pennmedicine.upenn.edu>
Brad Wubbenhorst
<bwubb@pennmedicine.upenn.edu>
