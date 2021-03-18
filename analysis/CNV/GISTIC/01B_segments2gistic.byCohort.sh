#Format Sequenza segments.txt files from a supplied list of tumors into one big cohort file, for use with GISTIC2.0
#Required input:
	#Only argument is a txt file of tumors of interest, one tumor per line
	#IMPORTANT: as written in this script, that file is "all_PR_WES_breast.txt." NEED TO EDIT FOR YOUR FILE NAME.
#outfile will be in GISTIC format
	#IMPORTANT: as written, outfile will be named "SegmentsforGistic.tsv". EDIT TO DESIRED OUTFILE NAME.
#Required organization:
	#To be run out of the directory containing subdirectories for each tumor
	#Will compile one master file for all segments.txt files corresponding to tumors in tumor list
	#If you want to generate a master file for all tumors in a directory, use segments2gistic.sh instead
#Seg.CN will be derived from depth ratio (see documentation)
	#awk doesn't have a log2 option, so here we use log(2*depth.ratio)/log(2)-1 to create the expression log2(depth.ratio)-1 
	#Num markers = N.BAF
#Relevant documentation
	#Converting Sequenza output to GISTIC input : http://crazyhottommy.blogspot.com/2017/11/run-gistic2-with-sequenza-segmentation.html
	#GISTIC2.0 : https://www.genepattern.org/modules/docs/GISTIC_2.0
#!/bin/bash
printf "Sample\tChromosome\tStart Position\tEnd Position\tNum markers\tSeg.CN\n" > SegmentsforGistic.breast.tsv
while read -r line; do dir=$(printf "%s\n" "$line"); file=$(ls "$dir" | grep segments.txt); awk -v dir="$dir" 'FNR>1 {print dir"\t"$1"\t"$2"\t"$3"\t"$5"\t"(log(2*$7)/log(2)-1)}' "$dir"/"$file"; done < all_PR_WES_breast.txt >> SegmentsforGistic.breast.tsv 