#Helper script to parse lines of variant call input that have multiple genes (i.e, Gene.refGene column contains something like "ABC;DEF")
	#Will create one new line per gene, for those lines that have multiple genes 
#To be called from make_rnk_file_variants.sh ; place in working directory for that script
#Input from stdin
#Output to stdout

import csv
import sys

lines=sys.stdin.readlines()
for line in lines:
	# parse the line to make a list (of 2 elements)
	line_parsed=line.split('\t')
	# check if the 2nd item in the list (the gene field) has semicolons (indicating multiple genes on one line) 
	if ";" in line_parsed[1]:
		#parse the gene list based on ";", then give each gene its own line in the output. The first field is still the patient ID.
		genes=line_parsed[1].split(';')
		for gene in genes:
			sys.stdout.write((line_parsed[0].strip())+'\t'+(gene.strip())+'\n')
	#if no parsing needed, don't do anything to the line 
	else :
		sys.stdout.write(line)