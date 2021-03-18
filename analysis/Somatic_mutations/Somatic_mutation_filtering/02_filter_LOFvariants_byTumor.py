#Filter an already-filtered list of variants one step further, for probable loss of function variants 
#Required inputs
	#1. Argument -i : a tsv file of Annovar-annotated somatic variants already filtered for alternative allele depth, population frequency, presence in exon, etc.
		#This file can be the outfile of filter_variants_byTumor.py 
#Outfile will be filtered for probable loss of function variants
	#Same format and fields as infile
	#Naming convention for outfile is <infile_string>.LOF.tsv
	#Will have >1 gene per line in some cases
		#this will be fixed if you tally up the most common mutations with make_rnk_file_variants.sh; see also one_gene_per_line.py

import csv
import os
import argparse

#Take in a list of variants (can be one tumors' variants or many tumors' variants concatenated into one file)
parser=argparse.ArgumentParser()
parser.add_argument('-i',action='store',dest='input',help='Enter the filepath of your already-filtered variants file')
args=parser.parse_args()

#Check input
if args.input==None:
	print("Error: please specify an input file of variants with -i")
	exit()
if not os.path.exists(args.input):
	print("Input file does not exist")
	exit()

	#read in variants line by line from input file, selecting for those likely to be loss of function
with open(args.input,'r') as infile, open((os.path.splitext(args.input)[0])+'.LOF.tsv','w') as outfile:
	reader=csv.DictReader(infile,delimiter='\t')
	writer=csv.DictWriter(outfile,delimiter='\t',fieldnames=reader.fieldnames)
	writer.writeheader()
	for row in reader:
		#select for exonic; leaving out any splice site mutations (pathogenicity can be hard to predict for these)
		if row["GenomicRegion.refGene"]=='exonic':
			#select for stopgains (nonsense) and frameshifts; filtering out non-frameshift insertions/deletions and stop loss variants
			if row["ExonicFunc.refGene"]=="frameshift_deletion" or row["ExonicFunc.refGene"]=="frameshift_insertion" or row["ExonicFunc.refGene"]=="stopgain":
				writer.writerow(row)
			#select for nonsynonymous SNVs predicted to be disease-causing via REVEL (filtering out those with no REVEL entry or score<=0.5)
			elif row["ExonicFunc.refGene"]=="nonsynonymous_SNV":
				if row["REVEL_score"]!=".":
					if float(row["REVEL_score"])>0.5:
						writer.writerow(row)
					else:
						()
				else:
					()
			else:
				()