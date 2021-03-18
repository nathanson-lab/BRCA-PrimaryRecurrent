#To be run out of the terminal on a list of files
#Required input: segments_file = the segments.txt file of Sequenza segments

getHRDstats<-function(segments_file)

{

library(HRDex)
source("referenceData.R") #modified chromosome col (from "chr1" to "1" format and factored variable)
source("preprocessHRD.R") #modified from source

#Read in raw sequencing data from file   
  seq_df<-read.delim(file=segments_file,sep = "\t",header=TRUE)
  seq_df$chromosome<-factor(seq_df$chromosome,levels=levels(grch37.ref.dat$chromosome))
  
#Preprocess it for HRD score calculation
  seq.dat<-preprocessHRD(seq_df, ref="grch37")
  CN.dat<-getCNt(seq.dat)
  
#Compute the scores
  HRD.Score.sum<-getHRD.Score(seq.dat,CN.dat,type = "sum")
  HRD.Score.average<-getHRD.Score(seq.dat,CN.dat,type = "average")
  LST.raw<-getLST(seq.dat)
  LOH.raw<-getLOH(seq.dat)
  NTAI.norm<-getNTAI.norm(seq.dat,CN.dat)
  NTAI.raw<-getNTAI.raw(seq.dat)
  
#Put results together
  hrd.stats<-data.frame(
    HRD.Score.sum=HRD.Score.sum,
    HRD.Score.average=HRD.Score.average,
    LST.raw=LST.raw,
    LOH.raw=LOH.raw,
    NTAI.norm=NTAI.norm,
    NTAI.raw=NTAI.raw
  )
  
  return(hrd.stats)
}