#Script to make bed files out of the significantly altered regions from GISTIC2.0
#infile should be amp_genes.conf_90.txt or del_genes.conf_90.txt files
#outfile is the desired filepath/name of your bed file
  #Resulting bed file will be UNSORTED; make sure to sort downstream as needed
  #Intended for use as input to compareBeds.sh
#Will select for significantly altered regions with residual q<0.05.

compileBeds<- function(infile,outfile){
  library(tidyverse)
  in.df<-read.delim(infile,sep="\t",header=FALSE)
  #restructure the df to only contain columns we want, for easier filtering
  out.df<-data.frame(
    residual.q=as.numeric(in.df[3,2:ncol(in.df)]),
    wide.peak.boundaries=as.character(in.df[4,2:ncol(in.df)]))
  #pull out only peaks with residual q<0.05
  peaks<-as.character(out.df %>% 
                        dplyr::filter(residual.q <0.05) %>% 
                        dplyr::pull(wide.peak.boundaries))
  #build the df for the bed file around your significant peaks
  bed.df<-data.frame("peak"=peaks)
  #split the peak column into 3 cols for bed format (chromosome, start, end)
  bed.df<-separate(bed.df,col=peak,sep=":",into=c("chromosome","start"))
  bed.df<-separate(bed.df,col=start,sep="-",into=c("start","end"))
  write_tsv(bed.df,path=outfile,col_names=FALSE,append=FALSE)
}
