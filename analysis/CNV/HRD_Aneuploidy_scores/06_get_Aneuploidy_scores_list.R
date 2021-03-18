#Assess for a list of tumors:
  #1. aneuploidy status by chromosome arm 
  #2. aneuploidy score for the tumor, defined as the number of chromosome arms
    #with aneuploidy not equal to that of the tumor (aneuploidy status=TRUE)
#Required inputs:
  #Txt file of primary tumor names; one tumor per line
  #Txt file of recurrent tumor names; one tumor per line
  #Segments.txt files (Sequenza) for each tumor, named <tumor name>_segments.txt
  #Ploidy estimates (Sequenza) for each tumor, named <tumor name>_confints.CP.txt
#Required organization
  #segments.txt files in "segments.txt.files/" directory
  #ploidy files in "confints_CP.txt.files/" directory
#Output 
  #1. File of aneuploidy status (TRUE/FALSE) by chromosome arm for each tumor
    #each named <tumor name>_Aneuploidy.txt
  #2. One txt file containing aneuploidy scores for primary tumor cohort
  #3. One txt file containing aneuploidy scores for recurrent tumor cohort

library(dplyr)
source("get_Aneuploidy_scores.R")

#Read in list of primary tumors; reconfigure the data structure to character vector 
primaries=read.delim(file="../tumor_lists/all_PR_WES_primaries.txt",header=FALSE,sep="\n")
primaries_list<-as.character(primaries[[1]])

#Make an empty dataframe in which to store aneuploidy scores for this list of tumors
batch_scores<-data.frame("tumor"=character(),"aneuploidy_score"=numeric())

for (i in 1:(length(primaries_list))) {
  segments_file<-paste0("segments.txt.files/",primaries_list[i],"_segments.txt")
  ploidy_file<-paste0("confints_CP.txt.files/",primaries_list[i],"_confints_CP.txt")
  #Make a dataframe of aneuploidy status by chromosome arm; write to a file  
  results<-get_Aneuploidy_scores(segments_file,ploidy_file)
  write_delim(results,path=paste0("./",primaries_list[i],"_Aneuploidy.txt"),append=FALSE,delim="\t")
  #Add the aneuploidy score for this tumor into the batch_scores
  aneuploidy_score<-as.numeric(sum(results$p.aneuploid == TRUE ) + sum(results$q.aneuploid == TRUE)) 
  batch_scores<-bind_rows(batch_scores,data.frame("tumor"=primaries_list[i],"aneuploidy_score"=aneuploidy_score))
}

#Write aneuploidy scores to a file for the entire batch 
write_delim(batch_scores,path="Aneuploidy_scores_primaries.txt",append=FALSE,delim="\t")

#same workflow for recurrences
recurrences=read.delim(file="../tumor_lists/all_PR_WES_recs.txt",header=FALSE,sep="\n")
recs_list<-as.character(recurrences[[1]])
batch_scores<-data.frame("tumor"=character(),"aneuploidy_score"=numeric())

for (i in 1:(length(recs_list))) {
  segments_file<-paste0("segments.txt.files/",recs_list[i],"_segments.txt")
  ploidy_file<-paste0("confints_CP.txt.files/",recs_list[i],"_confints_CP.txt")
  results<-get_Aneuploidy_scores(segments_file,ploidy_file)
  write_delim(results,path=paste0("./",recs_list[i],"_Aneuploidy.txt"),append=FALSE,delim="\t")
  aneuploidy_score<-as.numeric(sum(results$p.aneuploid == TRUE ) + sum(results$q.aneuploid == TRUE)) 
  batch_scores<-bind_rows(batch_scores,data.frame("tumor"=recs_list[i],"aneuploidy_score"=aneuploidy_score))
}

write_delim(batch_scores,path="Aneuploidy_scores_recs.txt",append=FALSE,delim="\t")
