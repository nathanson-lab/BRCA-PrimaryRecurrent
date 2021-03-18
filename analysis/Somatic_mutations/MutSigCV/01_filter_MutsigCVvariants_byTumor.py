#Filter a list of Annovar-annotated variants for MutSigCV (can be from one tumor or multiple, but must be one file)
#Required inputs:
	#1. Argument -i : path to tsv file of unfiltered variants
#Criteria for filtering (want to leave IN nonsynonymous variants or the algorithm will error out)
	#Popfreq_max<0.01
	#Tumor_ALT_AlleleDepth>=5
#Outfile will be filtered for mutations meeting criteria above
	#Same format and fields as infile
	#Naming convention for outfile is <infile_string>.filtered.MutSigCV.tsv

import csv
import os
import argparse

#Take in a list of variants, which can be one tumor's unfiltered file or many tumors' variants. Must be concatenated into one file.
parser=argparse.ArgumentParser()
parser.add_argument('-i',action='store',dest='input',help='Enter the filepath of your variants file')
args=parser.parse_args()

#Check input
if args.input==None:
	print("Error: please specify an input file of variants with -i")
	exit()
if not os.path.exists(args.input):
	print("Input file does not exist")
	exit()

#Read each line of input into the output based on the desired filtering criteria. Exclude rows that do not meet all filtering criteria.
with open(args.input,'r') as infile, open((os.path.splitext(args.input)[0])+'.filtered.MutSigCV.tsv','w') as outfile:
	reader=csv.DictReader(infile,delimiter='\t')
	writer=csv.DictWriter(outfile,delimiter='\t',fieldnames=reader.fieldnames)
	writer.writeheader()
	for row in reader:
		if int(row['Tumor_ALT_AlleleDepth'])>=5: #select for Alt AD >=5
			if row['Popfreq_max']=="." or float(row['Popfreq_max'])<0.01: #select for entries without popfreq_max entry (".") or <0.01
				writer.writerow(row) 
			else:
				()
		else:
			()