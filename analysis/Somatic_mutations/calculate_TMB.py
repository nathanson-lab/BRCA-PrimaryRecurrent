#Calculate tumor mutational burden (TMB) for each tumor in a list, using unfiltered tsv files of variants for each tumor
#TMB = (# of nonsynonymous SNVs + indels)/(Mbp size of the capture used for whole exome sequencing)
#Required inputs:
	#1. Argument -i : a txt file of tumors and their capture, where each line looks like "3939-Brca1Br169,v5"
		#Captures must be one of: v5, v6+COSMIC, or v7. 
		#Captures above refer to the SureSelect All Exon captures from Agilent. Script is written for hg19 capture sizes.
	#2. Argument -o : path to outfile, a tsv file containing fields for tumor (col 1) and mutational burden (col 2)
#Required organization:
	#To be run out of a directory in which each tumor has a subdirectory
	#Each tumor in the input tumor list should have a file called all.somatic.variants.tsv inside its subdirectory
	#The all.somatic.variants.tsv files should be Annovar-annotated and unfiltered
#Script will only count variants for TMB if they meet the following criteria: 
	#exonic (must be in "well-covered" area for capture)
	#Alt AD>=5, to get rid of artifact
	#nonsynonymous SNV or indel

import csv
import argparse
import os

#parse and check the tumor list file 
parser=argparse.ArgumentParser()
parser.add_argument('-i',action='store',dest='infile',help='Enter the filepath of your txt file of tumors and their captures')
parser.add_argument('-o',action='store',dest='outfile',help='Enter the filepath of your tsv file of mutational burden output by tumor')
args=parser.parse_args()

if args.infile==None:
	print("Error: specify input filepath with -i")
	exit()
if args.outfile==None:
	print("Error: specify output filepath with -o")
	exit()

#Open a new file for output
with open(args.outfile,'w') as outfile:
	writer=csv.writer(outfile,delimiter='\t')
	writer.writerow(['Tumor','Tumor_mutational_burden'])

with open(args.infile,'r') as tumor_list, open(args.outfile,'a') as outfile:
	for line in tumor_list:
		parsed=line.split(",")
		#Find the name and capture for the tumor on a given line; specify capture size (all sizes for hg19)
		tumor=parsed[0].strip()
		capture=parsed[1].strip()
		if capture=="v5":
			capture_size=50.00
		elif capture=="v6+COSMIC":
			capture_size=66.00
		elif capture=="v7":
			capture_size=48.2
		else:
			print("Invalid capture specified for "+tumor+". Quitting...")
		
		#Count the number of lines that are exonic and not synonymous variants; make sure to take from the UNFILTERED list of variants
		if os.path.exists(tumor+'/all.somatic.variants.tsv'):
			with open(tumor+'/all.somatic.variants.tsv','r') as variants:
				reader=csv.DictReader(variants,delimiter='\t')
				count=0
				for row in reader:
					if row["GenomicRegion.refGene"]=='exonic' or row["GenomicRegion.refGene"]=='exonic;splicing' or row["GenomicRegion.refGene"]=='ncRNA_exonic':
						if row["ExonicFunc.refGene"]!='synonymous_SNV':
							if float(row["Tumor_ALT_AlleleDepth"])>=5:
								count+=1
			#Calculate the tumor mutational burden
			tmb=round(count/capture_size,2)
			#Write to outfile
			writer=csv.writer(outfile,delimiter='\t')
			writer.writerow([tumor,tmb])
		else:
			print("No variants file found for "+tumor)