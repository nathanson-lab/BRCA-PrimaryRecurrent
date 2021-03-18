#Format Sequenza segments.txt files to one big cohort file of segments, for use with GISTIC 2.0
#No inputs or arguments
#outfile will be in GISTIC format and named SegmentsforGistic.tsv
#Required organization:
	#To be run out of the directory containing subdirectories for each tumor
	#Will compile one master file for all segments.txt files in the subdirectories 
	#If you want to limit the master file to a given list of tumors, use segments2gistic.byCohort.sh instead
#Seg.CN will be derived from depth ratio (see documentation)
	#awk doesn't have a log2 option, so here we use log(2*depth.ratio)/log(2)-1 to create the expression log2(depth.ratio)-1 
	#Num markers = N.BAF
#Relevant documentation
	#Converting Sequenza output to GISTIC input : http://crazyhottommy.blogspot.com/2017/11/run-gistic2-with-sequenza-segmentation.html
	#GISTIC2.0 : https://www.genepattern.org/modules/docs/GISTIC_2.0

#!/bin/bash
#open up the master file and write header
printf "Sample\tChromosome\tStart Position\tEnd Position\tNum markers\tSeg.CN\n" > SegmentsforGistic.tsv
#Num markers = N.BAF; Seg.CN = depth.ratio (but we will log transform this)
for file in $(ls -R */*_segments.txt); do dir=$(dirname "$file"); awk -v dir="$dir" 'FNR>1 {print dir"\t"$1"\t"$2"\t"$3"\t"$5"\t"(log(2*$7)/log(2)-1)}' "$file" >> SegmentsforGistic.tsv; done 