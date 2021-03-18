#Take in a master list of gene-annotated CNVs of a given CNV type for a cohort and rank genes by frequency for GSEA Preranked
	#For example, "identify genes most commonly subjected to copy number deletions for breast tumors in the cohort"
#Required input
	#Only argument is a cohort file (*.CN*.tsv) of annotated CNVs of a given type (an outfile from compile_annotSV_results.py) 
		#Each line in this file should be annotated with patient (col 1) and tumor ID (col 9) as a result of runAnnotSV.py
#Output
	#Outfile is a <infile string>_allGenes.rnk file that can be used as input for GSEA Preranked 
	#Outfile has 2 columns: 
		#Col 1 : Genes
		#Col 2 : Number of patients in which the gene is subject to that CNV type 
			#A gene can only be counted once per patient, to limit skewing towards patients with multiple tumors

#!/bin/bash

#check for file input
file=$(basename $1 .tsv)
if test -z $file; then
	echo "Error: no master file name given" #if no file name entered as argument
	exit 1
else 
	# take the patient and gene columns from the input ONLY; skip header
	tail -n+2 $1 | \
	cut -f 1,9 $1 | \
	# sort for unique lines (gene/patient combinations) --> get unique occurences of 1 gene per patient
		#Here we are discarding multiple occurrences of the same gene in one patient (due to it appearing in CNVs from multiple tumors from that patient)
	sort -u -k 1,1 -k 2,2 -s | \
	#Cut out the patient column so we just have a list of genes; sort this list 
	cut -f2 | sort | \
	#Annotate unique entries with # of occurrences in the file; need to employ sed to format uniq -c output
	uniq -c | sed -e 's/ *//' -e 's/ /\t/' | \
	sort -nr -k 1 > "$file"_allGenes.tsv
	#switch the order of the columns in all_Genes.tsv to get in right format for rnk file
	paste "$file"_allGenes.tsv "$file"_allGenes.tsv | cut -f2,3 > "$file"_allGenes.rnk
	#Sort by occurrence of each gene; write the output to a file (no header)
	#Get rid of allGenes.tsv file (redundant)
	rm "$file"_allGenes.tsv
	exit 0
fi
