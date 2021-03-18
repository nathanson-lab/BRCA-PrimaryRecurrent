#Script to run vcf2maf.pl (MSKCC program) given a vcf file 
#Required input:
	#Argument: only argument should be the file path of your input vcf. 
		#Should have tumor ID in the name of the file; will be used for outfile naming
	#OPTIONAL: vep may fail for very large (>1000 lines) files if input vcf is not sorted. 
		#To sort while keeping the header intact, run command below:
		##$for file in $(ls */*MutSigCV.vcf); do head -3 $file > temp1.tsv; tail -n+4 $file > temp2.tsv; sort -k1,1n -k2,2n temp2.tsv > temp2.sorted.tsv; cat temp1.tsv temp2.sorted.tsv > $file; rm temp1.tsv; rm temp2.tsv; rm temp2.sorted.tsv; done
#Required organization
	# --ref-fasta points to fasta file to use with vep
	# --vep-path points to location of ensembl-vep installation
	# --tumor-id points to the name of the vep file that will be created mid-script; don't change
	#To be run out of directory in which each tumor has its own subdirectory (with infile inside)
#Outfiles:
	#1. File of variants in maf format, placed in the same directory and named <infile_string>.maf
	#2. File of vep-annotated variants in vcf format, placed in same directory and named <infile_string>.vep.vcf
#Relevant documentation for vcf2maf: https://github.com/mskcc/vcf2maf

#!/bin/bash
perl ~/software/mskcc-vcf2maf/vcf2maf.pl \
--input-vcf $1 \
--output-maf $(echo $(dirname $1)/$(basename $1 .vcf).maf) \
--tumor-id $(echo $(basename $1 _all_variants.vcf)) \
--ref-fasta ~/.vep/homo_sapiens/99_GRCh37/Homo_*.fa \
--vep-path ~/software/ensembl-vep \
/