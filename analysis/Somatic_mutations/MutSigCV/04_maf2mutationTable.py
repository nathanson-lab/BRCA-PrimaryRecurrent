#Format a maf file of variants (from vcf2maf) into mutation table format for use with MutSigCV
#Required input:
	#1. Argument -i : infile with (at least) fields for Huge_Symbol (gene name), Tumor_Sample_Barcode (the tumor name), and Variant_Classification
#Required organization:
	#To be run out of a directory in which each tumor has a subdirectory. 
	#Maf infile should live in the subdirectory for its tumor (subdirectory name will be used for Tumor_Sample_Barcode column in outfile)
#Outfile will have all the same fields as input maf, plus an added column "effect" 
	#"effect" is derived here from Variant_Classification field, and must be one of "nonsilent", "silent", "noncoding"
	#Will also fix the Tumor_Sample_Barcode column from input and make it the tumor name (name of the directory of the infile)
	#Naming convention for outfile: <tumor subdirectory name>.all.somatic.variants.filtered.MutSigCV.table.maf
#Relevant documentation for MutSigCV input format: https://software.broadinstitute.org/cancer/cga/mutsig_run 

import csv
import os
import argparse

#Take in a maf file (can be one tumors' variants or many tumors' variants concatenated into one file)
parser=argparse.ArgumentParser()
parser.add_argument('-i',action='store',dest='input',help='Enter the maf filepath')
args=parser.parse_args()

#Check input
if args.input==None:
	print("Error: please specify an input maf file of variants with -i")
	exit()
elif not os.path.exists(args.input):
	print("Input file does not exist")
	exit()
else:
	()

#Name output file named based on input file; want tumor in the name (taken from directory of the file)
output_path=((os.path.dirname(args.input))+'/'+(os.path.dirname(args.input))+'.all.somatic.variants.filtered.MutSigCV.table.maf')

#List the fields in the header for input file 
#infields=['Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chromosome','Start_Position','End_Position','Strand','Variant_Classification','Variant_Type','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','dbSNP_RS','dbSNP_Val_Status','Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode','Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2','Tumor_Validation_Allele1','Tumor_Validation_Allele2','Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2','Verification_Status','Validation_Status','Mutation_Status','Sequencing_Phase','Sequence_Source','Validation_Method','Score','BAM_File','Sequencer','Tumor_Sample_UUID','Matched_Norm_Sample_UUID','HGVSc','HGVSp','HGVSp_Short','Transcript_ID','Exon_Number','t_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count','all_effects','Allele','Gene','Feature','Feature_type','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','ALLELE_NUM','DISTANCE','STRAND_VEP','SYMBOL','SYMBOL_SOURCE','HGNC_ID','BIOTYPE','CANONICAL','CCDS','ENSP','SWISSPROT','TREMBL','UNIPARC','RefSeq','SIFT','PolyPhen','EXON','INTRON','DOMAINS','AF','AFR_AF','AMR_AF','ASN_AF','EAS_AF','EUR_AF','SAS_AF','AA_AF','EA_AF','CLIN_SIG','SOMATIC','PUBMED','MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','IMPACT','PICK','VARIANT_CLASS','TSL','HGVS_OFFSET','PHENO','MINIMISED','ExAC_AF','ExAC_AF_AFR','ExAC_AF_AMR','ExAC_AF_EAS','ExAC_AF_FIN','ExAC_AF_NFE','ExAC_AF_OTH','ExAC_AF_SAS','GENE_PHENO','FILTER','flanking_bps','vcf_id','vcf_qual','ExAC_AF_Adj','ExAC_AC_AN_Adj','ExAC_AC_AN','ExAC_AC_AN_AFR','ExAC_AC_AN_AMR','ExAC_AC_AN_EAS','ExAC_AC_AN_FIN','ExAC_AC_AN_NFE','ExAC_AC_AN_OTH','ExAC_AC_AN_SAS','ExAC_FILTER','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF','vcf_pos']

#Construct a header for the output file, adding in a field for 'effect'
#outfields=infields.copy()
#outfields.append('effect')

#Make a list of noncoding variant classifications to check for when populating 'effect' field
noncoding=["3'Flank","3'UTR","5'Flank","5'UTR","Intron","IGR"]

#read in variants line by line from input file, and write line by line to output file 
with open(args.input,'r') as infile, open(output_path,'w') as outfile:
	#skip first line of input (specifying vcf version)
	reader=csv.reader(infile,delimiter='\t')
	next(reader)
	reader=csv.DictReader(infile,delimiter='\t')
	outfields=list(reader.fieldnames)
	outfields.append('effect')
	writer=csv.DictWriter(outfile,delimiter='\t',fieldnames=outfields)
	#add header to outfile
	writer.writeheader()
	for row in reader:
		#set a value for effect field in output, based on Variant_Classification field in input
			#Missense mutations are all coded as nonsilent, so check the Consequence column to see if synonymous_variant is the only entry
		if row["Variant_Classification"]=="Silent" or row["Consequence"]=="synonymous_variant":
			effect='silent'
		elif row["Variant_Classification"] in noncoding:
			effect='noncoding'
		else:
			effect='nonsilent'
		row['effect']=effect
		row['Tumor_Sample_Barcode']=os.path.dirname(args.input)
		writer.writerow(row)