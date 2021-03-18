#Compile all AnnotSV results (annotated CNVs) for a cohort
#Required inputs
	#1. Argument -i : txt file of tumors in the cohort (one tumor per line)
	#2. Argument -prefix : string for naming 7 outfiles; should be a description of the cohort (i.e., "primaryTumors")
#Required organization
	#To be run out of a directory in which each tumor has a subdirectory named for it
	#AnnotSV results (annotated CNVs) are in dated AnnotSV directory and named segments_CN*.genes.tsv 
		#Should be outfiles from runAnnotSV.py
	#Each tumor should have one set of AnnotSV results per CNV type (OR no results)
		#Can be from different dates and in different AnnotSV directories
		#But must have a single set of segments_CN0.genes.tsv results (for example) per tumor, or will throw error
		#If segments_CN*.genes.tsv is missing for a tumor, this is okay (tumor may not have had any amplifications, for example)
#Output
	#7 master files of annotated segments for the cohort (one for each CNV type; see binSegments.sh for details)
		#all named with convention of <user-entered prefix>.CN*.tsv
	#logfile reporting out which tumors were missing annotated segments files of each type

import csv
import sys 
import argparse
import subprocess
import os
import glob

#Parse and check the arguments
parser=argparse.ArgumentParser()
parser.add_argument('-i',action='store',dest='infile',help='Enter the filepath of your txt file of tumors')
parser.add_argument('-prefix',action='store',dest='prefix',help='Enter the string you want in each output file name (i.e., primaryTumors)')
args=parser.parse_args()

if args.infile==None:
	print("Error: specify input filepath with -i")
	exit()
elif args.prefix==None:
	print("Error: specify string for filenaming i.e., 'primaries' using -prefix")
	exit()
else:
	()

#name your output files based on the user-entered string
CN0_outfile=args.prefix+".CN0.tsv" #deletions
CN01_outfile=args.prefix+".CN0-1.tsv" #deletions+losses
CN1_outfile=args.prefix+".CN1.tsv" #losses
CN23_outfile=args.prefix+".CN2-3.tsv" #neutral
CN4up_outfile=args.prefix+".CN4+.tsv" #gains and amplifications
CN45_outfile=args.prefix+".CN4-5.tsv" #gains
CN6up_outfile=args.prefix+".CN6+.tsv" #amplifications

#open up the files to overwrite previous results 
open(CN0_outfile, 'w').close()
open(CN01_outfile, 'w').close()
open(CN1_outfile, 'w').close()
open(CN23_outfile, 'w').close()
open(CN4up_outfile, 'w').close()
open(CN45_outfile, 'w').close()
open(CN6up_outfile, 'w').close()

#initialize some lists to store info about missing files
missing_CN0=[]
missing_CN01=[]
missing_CN1=[]
missing_CN23=[]
missing_CN4up=[]
missing_CN45=[]
missing_CN6up=[]

#For each line (tumor) in the input file, find the respective segments_CN*.genes.tsv file and write its contents to the master file for that bin
with open(args.infile,'r') as input:
	for line in input:
		name=line.strip()
		#Deletions (CN=0)
		#find the file using glob; glob returns a list of length 1, so we will just save the first element in the list (the file name) as file variable
		file=(glob.glob(os.path.join(name,"*AnnotSV.results/segments_CN0.genes.tsv")))
		#glob only returns files that exist, so we will move on only if the list of the glob result is 1
		if not file: 
			print(name+" segments_CN0.genes.tsv file does not exist")
			missing_CN0.append(name)
		elif len(file)==1:
			CN0=open(file[0],'r')
			reader=csv.reader(CN0,delimiter='\t')
			with open(CN0_outfile,'a') as output:
				writer=csv.writer(output,delimiter='\t')
				for row in reader:
					writer.writerow(row)
		else:
			print("Multiple segments_CN0.genes.tsv files present for "+name+". Please delete or move duplicate files, then retry compilation. Quitting...")
			print(file)
			exit()
		
		#Deletions + Losses (CN=0-1)
		file=(glob.glob(os.path.join(name,"*AnnotSV.results/segments_CN0-1.genes.tsv")))
		if not file:
			print(name+" segments_CN0-1.genes.tsv file does not exist")
			missing_CN01.append(name)
		elif len(file)==1:
			CN01=open(file[0],'r')
			reader=csv.reader(CN01,delimiter='\t')
			with open(CN01_outfile,'a') as output:
				writer=csv.writer(output,delimiter='\t')
				for row in reader:
					writer.writerow(row)
		else:
			print("Multiple segments_CN0-1.genes.tsv files present for "+name+". Please delete or move duplicate files, then retry compilation. Quitting...")
			exit()
		
		#Losses (CN=1) 
		file=(glob.glob(os.path.join(name,"*AnnotSV.results/segments_CN1.genes.tsv")))
		if not file:
			print(name+" segments_CN1.genes.tsv file does not exist")
			missing_CN1.append(name)
		elif len(file)==1:
			CN1=open(file[0],'r')
			reader=csv.reader(CN1,delimiter='\t')
			with open(CN1_outfile,'a') as output:
				writer=csv.writer(output,delimiter='\t')
				for row in reader:
					writer.writerow(row)
		else:
			print("Multiple segments_CN1.genes.tsv files present for "+name+". Please delete or move duplicate files, then retry compilation. Quitting...")
			exit()
		
		#Neutral (CN=2-3)
		file=(glob.glob(os.path.join(name,"*AnnotSV.results/segments_CN2-3.genes.tsv")))
		if not file:
			print(name+" segments_CN2-3.genes.tsv file does not exist")
			missing_CN23.append(name)
		elif len(file)==1:
			CN23=open(file[0],'r')
			reader=csv.reader(CN23,delimiter='\t')
			with open(CN23_outfile,'a') as output:
				writer=csv.writer(output,delimiter='\t')
				for row in reader:
					writer.writerow(row)
		else:
			print("Multiple segments_CN2-3.genes.tsv files present for "+name+". Please delete or move duplicate files, then retry compilation. Quitting...")
			exit()
			
		#Gains + Amplifications(CN=4+)
		file=(glob.glob(os.path.join(name,"*AnnotSV.results/segments_CN4+.genes.tsv")))
		if not file:
			print(name+" segments_CN4+.genes.tsv file does not exist")
			missing_CN4up.append(name)
		elif len(file)==1:
			CN4up=open(file[0],'r')
			reader=csv.reader(CN4up,delimiter='\t')
			with open(CN4up_outfile,'a') as output:
				writer=csv.writer(output,delimiter='\t')
				for row in reader:
					writer.writerow(row)
		else:
			print("Multiple segments_CN4+.genes.tsv files present for "+name+". Please delete or move duplicate files, then retry compilation. Quitting...")
			exit()
			
		#Gains (CN=4-5)
		file=(glob.glob(os.path.join(name,"*AnnotSV.results/segments_CN4-5.genes.tsv")))
		if not file:
			print(name+" segments_CN4-5.genes.tsv file does not exist")
			missing_CN45.append(name)
		elif len(file)==1:
			CN45=open(file[0],'r')
			reader=csv.reader(CN45,delimiter='\t')
			with open(CN45_outfile,'a') as output:
				writer=csv.writer(output,delimiter='\t')
				for row in reader:
					writer.writerow(row)
		else:
			print("Multiple segments_CN4-5.genes.tsv files present for "+name+". Please delete or move duplicate files, then retry compilation. Quitting...")
			exit()

		#Amplifications (CN=6+)
		file=(glob.glob(os.path.join(name,"*AnnotSV.results/segments_CN6+.genes.tsv")))
		if not file:
			print(name+" segments_CN6+.genes.tsv file does not exist")
			missing_CN6up.append(name)
		elif len(file)==1:
			CN6up=open(file[0],'r')
			reader=csv.reader(CN6up,delimiter='\t')
			with open(CN6up_outfile,'a') as output:
				writer=csv.writer(output,delimiter='\t')
				for row in reader:
					writer.writerow(row)
		else:
			print("Multiple segments_CN6+.genes.tsv files present for "+name+". Please delete or move duplicate files, then retry compilation. Quitting...")
			exit()

#Get rid of multiple header lines in outfiles
subprocess.run("cat %s > temp.tsv; sort temp.tsv | uniq > %s" %(CN0_outfile,CN0_outfile), shell=True)
subprocess.run("cat %s > temp.tsv; sort temp.tsv | uniq > %s" %(CN01_outfile,CN01_outfile), shell=True)
subprocess.run("cat %s > temp.tsv; sort temp.tsv | uniq > %s" %(CN1_outfile,CN1_outfile), shell=True)
subprocess.run("cat %s > temp.tsv; sort temp.tsv | uniq > %s" %(CN23_outfile,CN23_outfile), shell=True)
subprocess.run("cat %s > temp.tsv; sort temp.tsv | uniq > %s" %(CN4up_outfile,CN4up_outfile), shell=True)
subprocess.run("cat %s > temp.tsv; sort temp.tsv | uniq > %s" %(CN45_outfile,CN45_outfile), shell=True)
subprocess.run("cat %s > temp.tsv; sort temp.tsv | uniq > %s" %(CN6up_outfile,CN6up_outfile), shell=True)
subprocess.run("rm temp.tsv", shell=True)

#Record which tumors had no genes reported for the various segments_CN*.bed files
with open(args.prefix+'.annotation.log.txt','w') as logfile:
	writer=csv.writer(logfile, delimiter=",")
	logfile.write("These tumors were missing segments_CN0.genes.tsv files: \n")
	writer.writerow(missing_CN0)
	logfile.write("\nThese tumors were missing segments_CN0-1.genes.tsv files: \n")
	writer.writerow(missing_CN01)
	logfile.write("\nThese tumors were missing segments_CN1.genes.tsv files: \n")
	writer.writerow(missing_CN1)
	logfile.write("\nThese tumors were missing segments_CN2-3.genes.tsv files: \n")
	writer.writerow(missing_CN23)
	logfile.write("\nThese tumors were missing segments_CN4+.genes.tsv files: \n")
	writer.writerow(missing_CN4up)
	logfile.write("\nThese tumors were missing segments_CN4-5.genes.tsv files: \n")
	writer.writerow(missing_CN45)
	logfile.write("\nThese tumors were missing segments_CN6+.genes.tsv files: \n")
	writer.writerow(missing_CN6up)
