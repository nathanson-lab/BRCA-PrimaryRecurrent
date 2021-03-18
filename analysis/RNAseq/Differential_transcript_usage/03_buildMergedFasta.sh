#Make a fasta file from the merged gtf file created in buildMergedGTF.sh
	#You will need this file for isoformSwitchAnalyzer
#Required 
	# merged gtf file generated in buildMergedGTF.sh (called merged.reference.gtf)
	# reference fasta for gencode v19 (GRCh37.p13.genome.fa; download from link below and unzip)
	# Samtools and gffread installed (see documentation below)
#Required organization
	# $Samtools should point to samtools install 
	# $Gffread should point to gffread install 
	# merged.reference.gtf and GRCh37.p13.genome.fa in working directory 
#Output will be a fasta file for the merged.reference.gtf (named merged.reference.fa)
#Relevant docs
	# Reference fasta for gencode v19 : pick "Genome sequence" from https://www.gencodegenes.org/human/release_19.html 
	# Gffread documentation https://github.com/gpertea/gffread
	# Samtools manual http://www.htslib.org/doc/#manual-pages

#!/bin/bash
#index the unzipped fasta file 
$Samtools faidx GRCh37.p13.genome.fa
#Create fasta for merged gtf file using the indexed reference fasta 
$Gffread merged.reference.gtf -g GRCh37.p13.genome.fa -w merged.reference.fa