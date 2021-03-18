#Run Stringtie for a single STAR-aligned RNA-seq bam file to generate Ballgown-formatted abundance files (see Stringtie documentation)
#Requires
	#Stringtie installed; use $Stringtie to point to location of install
	#Argument 1 is the name of the tumor subdirectory in which to place outfiles
	#Argument 2 is the path to the bam file from which to generate abundance files
	#Reference gtf file (see note below)
#Required organization
	#$RNA_seq points to directory in which each tumor has its own subdirectory
	#Reference gtf file (gencode.v19.annotation.gtf) also in $RNA_seq directory
#IMPORTANT: Annotation will be done with gencode_v19 transcriptome. 
	#This must be the same reference gtf/gff used for STAR alignment 
#Relevant docs
	#Stringtie manual: https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
	#GTF reference file: https://www.gencodegenes.org/human/release_19.html

#!/bin/bash
sample="$1"
bam="$2"

if test -z $sample; then
	echo "Error: no sample name given. Please enter sample name for output directory" #if no sample name entered as argument
	exit 1
elif test -z $bam; then
	echo "Error: no sample file given. Please enter bam input file" #if no bam file entered as argument
	exit 2
elif [[ $bam != *.bam ]]; then
	echo "Error: input file must be a bam file" #if file entered is not a bam
	exit 3
else
	$Stringtie $bam -eB -G ~/$RNA_seq/gencode.v19.annotation.gtf -o ~/$RNA_seq/$sample/stringtie.abundance.gtf -p 8
	exit 0
fi