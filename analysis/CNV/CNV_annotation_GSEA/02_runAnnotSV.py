#Annotate a binned file of Sequenza segments with genes using AnnotSV2.3
#Required input:
	#Only argument: any binned segments.bed file resulting from binSegments.sh
		#should have a number in the name, like *segments_CN0.bed
	#Edit bash .profile such that $ANNOTSV points to location of AnnotSV executable
#Arguments to AnnotSV (called by subprocess here; edit line 61 to change flags)
	#Tumors aligned to GRCh37
	#"split" annotation means one gene per line 
	#overlap argument is dictated by the infile name
		#For segments that are copy neutral, gains, or amplifications: overlap=100; genes only reported if 100% covered by CNV
		#For segments that are copy losses or deletions: overlap=50; genes reported if >=50% covered by CNV
		#This means (for example) that a gene is considered 
			# "deleted" if at least half the gene locus has CN=0 (as this would likely make the gene non-functional)
			# "amplified" only if entire gene locus has CN>=6 (as this would be required for increased gene dosage effect)
#Required organization
	#To be run out of a directory in which each tumor has a subdirectory, such that $dirname <file> would yield the tumor name
	#Tumors (subdirectories) should be named in the convention of <patient ID>-<tumor name>, for ex "3939-Brca1Br169"
#Output
	#new directory of AnnotSV results labeled with the date
	#outfile of annotated genes placed in that directory and named <infile string>.genes.tsv
		#each line annotated with patient and tumor ID, for easy appending/sorting downstream
	#a log file for the AnnotSV run
	#if runAnnotSV.py is run on multiple binned segments files from the same patient at once, all results go into same directory 


import argparse
import csv
import os.path
import os
import subprocess
import datetime

#Parse arguments; check input file
parser=argparse.ArgumentParser()
parser.add_argument('-f',action='store',dest='file',help='Enter the name of a binned CNV segments_CN*.bed file for annotation')
args=parser.parse_args()
if not args.file:
	print("Error: please enter a binned bed file name (of form */segments_CN*.bed) with -f")
	exit()
if not os.path.isfile(args.file):
	print("Error: the file name entered does not exist")
	exit()

#Pick a mode to run AnnotSV in, based on type of CNV
#gainsAmpsNeutral mode = whole gene must be covered by a segment to be counted as "included"
#delsLosses mode = half the gene must be covered by a segment to be counted as "included"
gainsAmpsNeutral_list=["CN2-3","CN4+","CN4-5","CN6+"]
delsLosses_list=["CN0","CN0-1","CN1"]
overlap=None

#check the file name for what type of CNVs it holds --> will dictate the call to AnnotSV
if any(substring in args.file for substring in gainsAmpsNeutral_list):
	print(args.file+" running in gainsAmpsNeutral mode...")
	overlap=100
elif any(substring in args.file for substring in delsLosses_list):
	print(args.file+ " running in delsLosses mode...")
	overlap=50
else:
	print("Incorrect file format: file must be a binned bed file of type */segments_CN*.bed. Exiting...")
	exit()

#make a tempIn file with the first 3 cols of input file (cat $ file | cut -f1-3 > temp.bed)
print("creating temporary file of input...")
subprocess.run("cat %s | cut -f1-3 > tempIn.bed" %(args.file),shell=True)

#name the directory and file for results
outputDir=((datetime.datetime.now()).strftime("%m%d"))+".AnnotSV.results/"
outputFile=os.path.splitext(args.file)[0]+".annotsv_out.tsv"
logFile=os.path.splitext(args.file)[0]+".annotsv_out.log.tsv"

#call to AnnotSV 
print("running annotation...")
#If AnnotSV is attempted on an empty file (just header), it will throw an error (that will be reported in log file for that file's annotation), then keep going
subprocess.run("$ANNOTSV/bin/AnnotSV -SVinputfile tempIn.bed -outputDir $(dirname %s)/%s -outputFile %s -overlap %s -typeOfAnnotation split -genomeBuild GRCh37 > %s" %(args.file,outputDir,outputFile,overlap,logFile),shell=True) 
#Get rid of the temporary file we created for annotation
print("Annotation finished, deleting temp file")
subprocess.run("rm tempIn.bed", shell=True)
#Move log file into results
print("Moving log file to the results folder "+outputDir+"...")
subprocess.run("mv %s $(dirname %s)/%s" % (logFile,args.file,outputDir), shell=True)

#If AnnotSV output has contents, we want to annotate each line with patient ID and tumor name 
infile=(os.path.dirname(args.file))+"/"+outputDir+os.path.basename(outputFile)
if os.path.getsize(infile)> 0:
	print("annotating "+infile+" with tumor and patient name...")
	outfile=(os.path.dirname(args.file))+"/"+outputDir+(os.path.splitext(os.path.basename(args.file))[0])+".genes.tsv" #final outfile name
	print("creating outfile named "+outfile)
	tumor=os.path.dirname(args.file) #directory of input segments.bed file should be the tumor name
	patient=(os.path.dirname(args.file)).split("-")[0] #name of the patient should be the first part of the directory name string (before first "-")
	with open(infile,'r') as input, open(outfile,'w') as output:
		reader=csv.reader(input, delimiter='\t')
		writer=csv.writer(output, delimiter='\t', lineterminator='\n')
		next(reader) #skip header
		writer.writerow(['Patient','Tumor','AnnotSV ID','SV chrom','SV start','SV end','SV length','AnnotSV type','Gene name','NM','CDS length','tx length','location','location2','intersectStart','intersectEnd','DGV_GAIN_n_samples_tested','DGV_GAIN_Frequency','DGV_LOSS_n_samples_tested','DGV_LOSS_Frequency','promoters','dbVar_event','dbVar_variant','dbVar_status','ACMG','morbidGenesCandidates','morbidGenes','pLI_ExAC','HI_CGscore','TriS_CGscore'])
		for row in reader:
			writer.writerow([patient] + [tumor] + row)
else:
	print("No genes encompassed by "+args.file+", annotation finished")