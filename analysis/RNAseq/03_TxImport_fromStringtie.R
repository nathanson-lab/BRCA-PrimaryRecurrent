#Import abundance files and plot read counts
#Required inputs
  #excel spreadsheet containing sample metadata (study_design.xlsx)
  #t_data.ctab files for each tumor (each the result of runStringtie.sh) 
  #tsv file with read counts by tumor (the result of getRNAreadcount.sh)
#Required organization
  #tumors should be named in convention of <patient ID>-<tumor ID>, but note that 
    #here "wsl_id" refers to the patient ID and "knm_id" refers to tumor ID
  #abundance files should be located in ./abundance_files/<tumor name>/t_data.ctab
  #Should start an R project with this script; will want to save objects generated
#Outputs
  #plot of read counts (bar graph)
  #targets dataframe (sample metadata to use for rest of pipeline)
  #imported transcript abundance for use in rest of pipeline

#load packages----
library(tidyverse)
library(tximport) 
library(readxl)
library(scales)

#Import study design file with all covariates of interest ----
targets <- read_xlsx("study_design.xlsx",col_names=TRUE,trim_ws=TRUE) #read in MOST UP TO DATE study design file
#sample column takes form "WSL"_"KNM" whereas directory names are "WSL"-"KNM", which is why we have to use paste below 
# set file paths to your t_data.ctab files
path <- file.path("abundance_files",paste0(targets$wsl_id,"-",targets$knm_id), "t_data.ctab") 
#check if files exist for all samples in study design file
all(file.exists(path))

#Import abundance files (generated with Stringtie, */t_data.ctab) using Tximport ----
#grab gene/transcript annotations from first .ctab file in the list
tmp <- read_tsv(path[1])
tx2gene <- tmp[, c("t_name", "gene_name")]

#import all files using the annotation above
Txi_gene <- tximport(path, 
                type = "stringtie", 
                tx2gene = tx2gene,
                txOut=FALSE, #want genes, not transcripts
                countsFromAbundance = "lengthScaledTPM") #counts scaled by avg transcript length + library size

# Add in read count data and plot ----
#import metrics (the output of getRNAreadcount.sh); cut off last few rows with averages, min, max
readCounts<-(read_tsv("../RNAseq.readCounts.tsv"))[1:66,]
#format tumor name for R (with "_" as wsl/KNM separator)
readCounts<-separate(readCounts,Bam_file,into=c("patient","tumor"),sep="-",remove=FALSE)
readCounts<-unite(readCounts,sample,patient:tumor,sep="_")
readCounts<-select(readCounts,-Bam_file)
#add read counts into targets df
targets<-inner_join(targets,readCounts,by="sample")

#Plot read counts for the cohort
ggplot(targets,aes(x=sample,y=Mapped_reads,fill=tissue_sample_type))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  scale_y_continuous(name="Number of Mapped Reads",labels=scientific,breaks = seq(0,3e+08,25e+06))+
  xlab("Sample")+
  ggtitle("Mapped Reads Across RNAseq Cohort")+
  geom_hline(yintercept = 25e06, linetype="longdash", colour="black", size=1)+
  scale_fill_discrete(name="Sample Type",labels=c("Normal Breast","Primary Breast Tumor",
                                                  "Recurrent Breast Tumor","Normal Fallopian Tube",
                                                  "Primary Ovarian Tumor","Recurrent Ovarian Tumor"))
  
  

