#Generate HRD scores from sequenza results for a cohort of samples
#Required input (modify in lines 13,17):
  #file path to txt file of primary tumor names (one tumor per line)
  #file path to txt file of recurrent tumor names (one tumor per line)
#Required organization
  #all segments.txt files (Sequenza) in "segments.txt.files" directory
  #all segments.txt files named in convention of <tumor name>_segments.txt
#Output
  #tab-delimited outfile named <tumor name>_HRDstats.txt with HRD stats for each tumor

library(readr)
#Source required function
#make sure you have modified (not HRDex) referenceData.R and preprocessHRD.R 
  #these scripts will be called by getHRDstats.R
source("getHRDstats.R")

#Read in list of primary tumors; reconfigure the data structure to character vector 
primaries=read.delim(file="../tumor_lists/all_PR_WES_primaries.txt",header=FALSE,sep="\n")
primaries_list<-as.character(primaries[[1]])

#Do the same for recurrences
recurrences=read.delim(file="../tumor_lists/all_PR_WES_recs.txt",header=FALSE,sep="\n")
recs_list<-as.character(recurrences[[1]])

for (i in 1:(length(primaries_list))) {
  segments_file<-paste0("segments.txt.files/",primaries_list[i],"_segments.txt")
  results<-getHRDstats(segments_file)
  write_delim(results,path=paste0("./",primaries_list[i],"_HRDstats.txt"),append=FALSE,delim="\t")
}

for (i in 1:(length(recs_list))) {
  segments_file<-paste0("segments.txt.files/",recs_list[i],"_segments.txt")
  results<-getHRDstats(segments_file)
  write_delim(results,path=paste0("./",recs_list[i],"_HRDstats.txt"),append=FALSE,delim="\t")
}