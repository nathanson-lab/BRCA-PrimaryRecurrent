#Run Stringtie for a single RNAseq bam file for differential transcript usage analysis to generate abundance files in Ballgown format
#Required inputs
	#Argument 1 : name of tumor (and its subdirectory)
	#Argument 2 : path to tumor RNAseq bam
	#merged.reference.gtf, a merged gtf (generated in buildMergedGTF.sh) --> used as the reference gtf file this time
		#IMPORTANT: this merged gtf should be built using only and exactly the files intended for DTU analysis 
		#should sit in working directory
#Required organization
	#To be run out of a directory in which each tumor has a subdirectory 
	# $Stringtie should point to location of Stringtie installation
#Output
	#New directory (forDTU/) of Ballgown-formatted abundance files to use for DTU analysis only

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
	$Stringtie $bam -eB -G merged.reference.gtf -o $sample/forDTU/stringtie.abundance.gtf -p 8
	exit 0
fi