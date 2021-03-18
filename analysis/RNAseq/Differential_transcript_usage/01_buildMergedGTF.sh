#Run Stringtie for a set of gtf files (from RNAseq bams; the results of runStringtie.sh) to build a merged gtf file  (merged mode)
	#This merged gtf will be used to run Stringtie again, specifically for differential transcript usage (DTU) in next step
#Required 
	# gtf files generated from runStringtie.sh
	# reference gtf (should be the same reference file used in runStringtie.sh and also in the initial STAR alignment)
	# Stringtie installed 
#Required organization
	# To be run out of a directory in which each tumor has a subdirectory
	# Each tumor subdirectory should contain a file named stringtie.abundance.gtf (resulting from runStringtie.sh with same reference gtf)
	# $Stringtie should point to the location of your Stringtie installation
	# Reference gtf should sit in working directory 
#Output will be a merged gtf file placed in working directory (named merged.reference.gtf)
	#we will use this as a reference transcriptome to re-generate abundance files specifically for differential transcript usage
	#see Stringtie documentation for details on this 
#Relevant docs
	#Stringtie manual: https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
	#GTF reference file: https://www.gencodegenes.org/human/release_19.html

#!/bin/bash
$Stringtie --merge \
-G gencode.v19.annotation.gtf \
-o merged.reference.gtf \
-p 8 \
$(ls */stringtie.abundance.gtf)