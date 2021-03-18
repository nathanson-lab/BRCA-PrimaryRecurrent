#Make a master list of somatic loss of function variants for a given cohort of tumors
#Script takes in a list of tumors, finds the appropriate file of all.somatic.variants.filtered.LOF.tsv for each, and appends to a master list for the cohort
#Required inputs 
	#1. Argument -i : txt file of tumors. One tumor name per line, named in the convention of "patient ID"-"tumor name", for example: 3939-Brca1Br169
	#2. Argument -o : name of outfile 
#Required file conventions
	#1. To be run out of a directory in which each tumor has its own subdirectory.
	#2. Inside each tumor subdirectory, LoF variants for that tumor should exist in file called all.somatic.variants.filtered.LOF.tsv
		#this file is the outfile of filter_LOFvariants_byTumor.py

import csv
import argparse
import os

#Parse and check the arguments
parser=argparse.ArgumentParser()
parser.add_argument('-i',action='store',dest='infile',help='Enter the filepath of your txt file of tumors')
parser.add_argument('-o',action='store',dest='outfile',help='Enter the name of your output file of compiled LOF variants')
args=parser.parse_args()

if args.infile==None:
	print("Error: specify tumor list filepath with -i")
	exit()
if args.outfile==None:
	print("Error: specify filepath for output file with -o")
	exit()
else:
	()

#start a list of tumors that are listed in input file but are missing a variants file
missing_variants_files=[]

#Open an empty file for writing output, then add a header 
with open(args.outfile,'w') as master_file:
	writer=csv.writer(master_file,delimiter='\t')
	writer.writerow(['Chr','Start','End','Ref','Alt','Popfreq_max','REVEL_score','Patient','TumorID','GenomicRegion.refGene','ExonicFunc.refGene','Gene.refGene','Exon.refGene','NTChange.refGene','AAChange.refGene','REVEL_rankscore','Tumor_Zyg','Tumor_Total_Depth','Tumor_ALT_AlleleDepth','Tumor_ALT_AlleleFrac','CLNSIG','NumPASS','MUTECT2','STRELKA2','VARDICTJAVA','VARSCAN2'])

#For each line (tumor) in the input file, find the */all.somatic.variants.filtered.LOF.tsv file and write it to the master outfile 
with open(args.infile,'r') as input:
	for line in input:
		tumor=line.strip()
		if os.path.exists(tumor+'/all.somatic.variants.filtered.LOF.tsv'):
			with open(tumor+'/all.somatic.variants.filtered.LOF.tsv','r') as tumor_file, open(args.outfile,'a') as master_file:
				reader=csv.reader(tumor_file,delimiter='\t')
				writer=csv.writer(master_file,delimiter='\t')
				#skip the header in the tumor variants file
				next(reader)
				#add every following line of tumor variants to the master file 
				for row in reader:
					writer.writerow(row)
		else:
			print(tumor+" variants file does not exist")
			missing_variants_files.append(tumor)

#Write the tumors missing variants to a log file 
if not missing_variants_files:
	print("Compilation complete, all files found")
else:
	with open((os.path.splitext(args.outfile)[0])+'_log.txt','w') as missing_variants_log:
		writer=csv.writer(missing_variants_log, delimiter=",")
		missing_variants_log.write("These tumors were missing variants files: \n")
		writer.writerow(missing_variants_files)