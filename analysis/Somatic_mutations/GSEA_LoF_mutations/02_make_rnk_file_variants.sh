#Take in a list of filtered variants for a cohort, rank genes by frequency of mutation, and produce a ranked gene list (*.rnk) for GSEA Preranked 
#Required inputs: 
	#1. Only argument: tsv file of filtered mutations for a cohort
		#For example, the output of compile_LOFvariants_byCohort.py (but could also be variants filtered a different way)
		#Col 8 of infile must be patient ID
		#Col 12 of infile must be RefSeq gene ID (Gene.refGene field from Annovar annotation)
	#2. One_gene_per_line.py script (called in line 23). Place this in working directory or edit line 23 to desired path
#Outfile is a 2-column list of genes annotated by # of patients in which they are mutated --> can be fed into GSEA Preranked as input 
	#Col 1: gene name
	#Col 2: number of occurrences of that gene in the supplied infile (each gene limited to once occurrence per patient) 

#!/bin/bash

#input should be a tsv file with somatic mutations across a cohort of patients. First use tail to remove the header so it doesn't get sorted 
infile=$1
outfile=$(echo -e $(basename $1 .tsv).allGenes.rnk)
tail -n +2 $infile | \
# take the patient and gene columns from the input ONLY
cut -f 8,12 | \
# sort for unique lines (gene/patient combinations) --> get unique occurences of 1 gene per patient
	#Here we are discarding multiple mutations of the same gene in the same patient (so a gene can only be counted as mutated once per patient)
sort -u -k 1,1 -k 2,2 -s | \
#Need to deal with lines that have multiple genes separated by ";" --> use python script to parse. Output has one gene per line.
python one_gene_per_line.py | \
#Cut out the patient column so we just have a list of genes; sort this list 
cut -f2 | sort | \
#Annotate unique entries with # of occurrences in the file; need to employ sed to format uniq -c output
uniq -c | sed -e 's/ *//' -e 's/ /\t/' | \
sort -nr -k 1 > temp.tsv
paste temp.tsv temp.tsv| cut -f2,3 > $outfile
rm temp.tsv 
#Sort by occurrence of each gene; write the output to a file (no header)