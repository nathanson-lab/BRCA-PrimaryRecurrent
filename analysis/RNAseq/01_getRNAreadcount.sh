#Get RNAseq read count from each STAR-aligned tumor bam file (optional, but will plot later to check for outliers)
#Based on STAR aligner usage, we know that all reads in these files should be mapped (unmapped reads were stored in *.fastx files elsewhere)
#So the total reads = mapped reads for all bams, but below we are still using a samtools command to just get mapped reads (just to be safe)
#Requires
	#samtools installed
	#no arguments or inputs needed for script itself
	#$RNA_seq points to directory containing tumor subdirectories and bams (see below)
#Required organization:
	#To be run out of $RNA_seq directory in which each tumor has its own subdirectory 
	#Inside each tumor subdirectory: STAR/input.bam
#Outfile is a tsv file containing tumor name (col 1) and number of mapped reads (col 2)

#!/bin/bash
#write header 
echo -e "Bam_file\tMapped_reads" > RNAseq.readCounts.tsv
for folder in $(ls $RNA_seq); do 
reads=$(samtools view -c -F 4 $RNA_seq/$folder/STAR/input.bam)
echo -e $folder$'\t'$reads >> RNAseq.readCounts.tsv;
done 
