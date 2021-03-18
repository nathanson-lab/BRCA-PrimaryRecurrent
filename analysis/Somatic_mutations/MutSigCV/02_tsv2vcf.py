#Format a *.tsv file of Annovar-annotated variants into vcf (v4.3) file for use as input to vcf2maf.pl
#Required inputs:
	#1. Argument -i : tsv file of variants from a single tumor, containing the following fields for each variant: Chr, Start, Ref, Alt
		#For example, the outfile of filter_MutsigCVvariants_byTumor.py
		#here we are also including Tumor_ALT_AlleleFrac field, but this isn't required for vcf2maf to run
#Outfile will be same variants but in pared-down vcf v4.3 format
#Relevant documentation
	#vcf2maf https://github.com/mskcc/vcf2maf
	#vcf v4.3 http://samtools.github.io/hts-specs/VCFv4.3.pdf

import csv
import os
import argparse
import subprocess

#Take in a tsv file (need one file per tumor; don't combine into batches)
parser=argparse.ArgumentParser()
parser.add_argument('-i',action='store',dest='input',help='Enter the tsv filepath')
args=parser.parse_args()

#Check input
if args.input==None:
	print("Error: please specify an input vcf file of variants with -i")
elif not os.path.exists(args.input):
	print("Input file does not exist")
else:
	()

#Name output file named based on input file, and construct header 
output_path=(os.path.splitext(args.input)[0])+'.vcf'

#Construct a header for the output file
outfields=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']

#read in variants line by line from input file, and write line by line to output file 
with open(args.input,'r') as infile, open(output_path,'w') as outfile:
	reader=csv.DictReader(infile,delimiter='\t')
	next(reader)
	#Write the header documentation required for vcf format
	outfile.write('##fileformat=VCFv4.3\n')
	outfile.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
	#Switch to DictWriter for the variant rows
	writer=csv.DictWriter(outfile,delimiter='\t',fieldnames=outfields)
	#add header to outfile (3rd line down)
	writer.writeheader()
	for row in reader:
		info="AF="+row['Tumor_ALT_AlleleFrac']
		newDict={'#CHROM':row['Chr'],
				'POS':row['Start'],
				'ID':'.',
				'REF':row['Ref'],
				'ALT':row['Alt'],
				'QUAL':'.',
				'FILTER':'PASS',
				'INFO':info}
		writer.writerow(newDict)