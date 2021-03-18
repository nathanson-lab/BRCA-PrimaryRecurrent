#Compare significantly altered regions from two bed files (including, but not limited to, GISTIC2.0 bed outfiles)
#Required inputs:
	#Arguments 1 and 2 : filepaths to bed files for comparison (order doesn't matter) 
		#Bed files should be 3-column bed-format (chromosome, start, end) for each peak, with no header. 
		#Chromosome input should be of type "chr1" rather than numeric ("1")
	#Argument 3 : name of outfile (.tsv)
#Output: a tsv file with Jaccard statistics for the intersection of bed files
	

#!/bin/bash
#Check number of arguments first
if [[ $# -ne 3 ]]; then
	echo "Too few arguments (need 3); exiting..."
	exit 1
else
	infile1=$1
	infile2=$2
	outfile=$3
	#clean up line endings as needed
	dos2unix $infile1
	dos2unix $infile2
	#sort the input bed files
	cat $infile1 > temp1.tsv; sort -k1,1 -k2,2n temp1.tsv > temp1.sorted.tsv
	cat $infile2 > temp2.tsv; sort -k1,1 -k2,2n temp2.tsv > temp2.sorted.tsv
	#run bedtools jaccard to get intersection statistics
	bedtools jaccard -a temp1.sorted.tsv -b temp2.sorted.tsv > $outfile
	rm temp*.tsv
	exit 0
fi
