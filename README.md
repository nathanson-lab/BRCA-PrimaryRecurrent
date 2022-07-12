# BRCA-PrimaryRecurrent


Scripts for somatic variant calling and allele-specific copy number profiling in matched tumor/normal whole exome sequencing data.

## Required Software:

* Python 3.7+  
   [https://www.python.org/downloads/]
   
* vcfpy  
   [https://vcfpy.readthedocs.io/en/stable/]
   
* Snakemake  
   [https://snakemake.readthedocs.io/en/stable/getting_started/installation.html]
   
* samtools  
   [https://github.com/samtools/samtools]
   
* bcftools  
   [https://github.com/samtools/bcftools]
   
* GATK 4  
   [https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2]
   
* Strelka2  
   [https://github.com/Illumina/strelka]
   
* VarDictJava  
   [https://github.com/AstraZeneca-NGS/VarDict]
   
* VarScan2  
   [https://github.com/dkoboldt/varscan]
   
* R  
   [https://www.r-project.org/]
   
* Sequenza (R Package)  
   [https://cran.r-project.org/web/packages/sequenza/index.html]
   
* Sequenza-utils  
   [https://bitbucket.org/sequenzatools/sequenza-utils/src/master/]

## Usage
```
snakemake -s mutect2.snake --configfile config.yaml
snakemake -s strelka2.snake --configfile config.yaml
snakemake -s vardictjava.snake --configfile config.yaml
snakemake -s sequenza.snake --configfile config.yaml
snakemake -s varscan2.snake --configfile config.yaml
```

## Authors
Jennifer Shah
<jennshah@pennmedicine.upenn.edu>

Brad Wubbenhorst
<bwubb@pennmedicine.upenn.edu>

## Paper Citation
