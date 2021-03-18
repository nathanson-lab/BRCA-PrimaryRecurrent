#Bin all CNV segments from a given tumor into 7 bins based on total CN
#Required inputs:
	#Only argument is the file to bin, a tab-separated bed file of segments from Sequenza.
		#The 5th field of this file =totalCN*100 for each segment
#Outfiles: 7 binned *.bed files per tumor (CN0,CN1,CN0-1,CN2-3,CN4-5,CN4+,CN6+)
	#CN0 = copy number deletions
	#CN1 = copy number losses
	#CN2-3 = copy neutral;combines segments of CN=2 and segments of CN=3 without intersection (just append lines)
	#CN4-5 = copy number gains; combines segments of CN=4 and segments of CN=5 without intersection 
	#CN6+ = copy number amplifications; combines segments with total CN of 6 or greater (inclusive) without intersection 
	#CN0-1 = copy number losses + deletions; combines segments of CN=0 and segments of CN=1 without intersection
	#CN4+ = copy number gains + amplifications; combines segments with total CN of 4 or greater (inclusive) without intersection
#Required organization:
	#To be run out of a directory in which each tumor has a subdirectory 
	#Outfiles will be named in the convention of <infile_string>_CN*.bed; need to keep this convention if running runAnnotSV.py
	#Each binned segments outfile will get placed in the tumor's directory

#!/bin/bash
file="$1"
if test -z $file; then
	echo "Error: no file name given. Please enter segments.bed file name" #if no file name entered as argument
	exit 1
else
	#Initialize each segments output file with a header
	printf "chr\tCNV_start\tCNV_end\tdescription\ttotalCN\tstrand\n" > $(dirname $file)/$(basename $file .bed)_CN0.bed
	printf "chr\tCNV_start\tCNV_end\tdescription\ttotalCN\tstrand\n" > $(dirname $file)/$(basename $file .bed)_CN1.bed
	printf "chr\tCNV_start\tCNV_end\tdescription\ttotalCN\tstrand\n" > $(dirname $file)/$(basename $file .bed)_CN0-1.bed
	printf "chr\tCNV_start\tCNV_end\tdescription\ttotalCN\tstrand\n" > $(dirname $file)/$(basename $file .bed)_CN2-3.bed
	printf "chr\tCNV_start\tCNV_end\tdescription\ttotalCN\tstrand\n" > $(dirname $file)/$(basename $file .bed)_CN4-5.bed
	printf "chr\tCNV_start\tCNV_end\tdescription\ttotalCN\tstrand\n" > $(dirname $file)/$(basename $file .bed)_CN4+.bed
	printf "chr\tCNV_start\tCNV_end\tdescription\ttotalCN\tstrand\n" > $(dirname $file)/$(basename $file .bed)_CN6+.bed

	##Deletions, Losses
	#print all CN=0 segments to CN0 file
	awk -F "\t" '{if($5==0){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' $file >> $(dirname $file)/$(basename $file .bed)_CN0.bed
	#print all CN=1 segments to CN1 file
	awk -F "\t" '{if($5==100){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' $file >> $(dirname $file)/$(basename $file .bed)_CN1.bed
	#print all CN=0 and CN=1 segments to CN0-1 file (deletions + losses)
	awk -F "\t" '{if($5==0 || $5==100){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' $file >> $(dirname $file)/$(basename $file .bed)_CN0-1.bed

	##Neutral
	#print all CN=2 and CN=3 segments to the CN2-3 file
	awk -F "\t" '{if($5==200 || $5==300){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' $file >> $(dirname $file)/$(basename $file .bed)_CN2-3.bed

	##Gains, Amplifications
	#print all CN=4 and CN=5 segments to the CN4-5 file
	awk -F "\t" '{if($5==400 || $5==500){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' $file >> $(dirname $file)/$(basename $file .bed)_CN4-5.bed
	#print all CN>=4 segments to CN4+ file (gains + amplifications)
	awk -F "\t" '{if($5>=400){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' $file >> $(dirname $file)/$(basename $file .bed)_CN4+.bed
	#print all CN>=6 segments to CN6+ file
	awk -F "\t" '{if($5>=600){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' $file >> $(dirname $file)/$(basename $file .bed)_CN6+.bed
	exit 0
fi