#Script to generate aneuploidy status (TRUE/FALSE) by chromosome from sequenza output
#Required inputs:
  #segments_file = path to Sequenza segments.txt file
  #ploidy_file = path to Sequenza ploidy estimate file (.txt)

get_Aneuploidy_scores<- function(segments_file,ploidy_file)
{

library(readr)
library(HRDex)
library(dplyr)
source("referenceData.R")

segments_df<-read.delim(segments_file,sep='\t',header=TRUE)
ploidy_df<-read.delim(ploidy_file,sep='\t',header=TRUE)

#Make empty dataframe to store the aneuploidy info for each chromosome 
results<-data.frame("p.aneuploid"=logical(),"q.aneuploid"=logical(),"chr"=character())
#Find aneuploidy status for each chromosome and save to dataframe
for (i in 1:22)
  {
  chr_results<-getAneuploidy(segments_df,ploidy_df,chr=i,max.brk.len=3e06)
  results<-bind_rows(results,chr_results[[1]])
  }

#Find chrX aneuploidy status and save to dataframe 
chr_results<-getAneuploidy(segments_df,ploidy_df,chr="X",max.brk.len = 3e06)
results<-bind_rows(results,chr_results[[1]])

return(results)

}