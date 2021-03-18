# Objectives -----------
#Analyze differential transcript usage (DTU) with IsoformSwitchAnalyzeR
#Identify top switches
#Make switch plots depicting switches of interest 
#Identify which samples express each isoform of BRCA2
#Required input
  # tsv file containing sample metadata for RNAseq samples ("../PR_RNAseq/targets/updated/tsv")
  # t_data.ctab files generated from runStringtieDTU.sh
    #File hierarchy: abundance_files_forDTU/<sample name>/forDTU/t_data.ctab
    #Sample name takes the form of <patient ID>-<tumor ID>
      #but here, Patient ID = "wsl_id" and Tumor ID = "knm_id"
  # merged.reference.gtf (from buildMergedGTF.sh)
  # merged.reference.fa (from buildMergedFasta.sh)
#Output
  #Top isoform switches by q-value (write to file)
  #Switch plots for switches of interest (grouped or individual)
  #SwitchList objects 
  #Dataframes of samples expressing each isoform of BRCA2

# Load packages -----
library(tidyverse) 
library(limma) 
library(edgeR)
library(gt) 
library(DT) 
library(plotly) 
library(IsoformSwitchAnalyzeR)
# https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html

# Import sample metadata and transcript data ----
# The IsoformSwitchAnalyzeR package looks for certain column headers in our study design -->
# unique sample IDs must be contained in column called 'sampleID'
# covariate(s) of interest must be in column labeled 'condition'
  #Our conditions for comparison here are normal tissue, primary tumor, and recurrence 
# remove extraneous columns
targets<-read_tsv("../PR_RNAseq/targets.updated.tsv",col_names = TRUE)
sampleLabels<-targets$sample
targets.mod <- targets %>%
  dplyr::rename(sampleID = sample, condition = sample_type) %>%
  dplyr::select(sampleID, condition)

# import transcript data, using new path variable
  #readLength=150 because we did 150 paired end sequencing 
path.DTU<-file.path("abundance_files_forDTU/",paste0(targets$wsl_id,"-",targets$knm_id),"forDTU/t_data.ctab")
Txi_trans <- importIsoformExpression(readLength=150, sampleVector = path.DTU)

# fix column headers of abundance and counts data to match sampleID in targets.mod
colnames(Txi_trans$abundance) <- c("isoform_id", sampleLabels) 
colnames(Txi_trans$counts) <- c("isoform_id", sampleLabels) 

# import data into SwitchList for analysis ----
mySwitchList <- importRdata(
  isoformCountMatrix   = Txi_trans$counts,
  isoformRepExpression = Txi_trans$abundance,
  designMatrix         = targets.mod,
  removeNonConvensionalChr = TRUE,
  addAnnotatedORFs=TRUE,
  isoformExonAnnoation = "merged.reference.gtf",
  isoformNtFasta       = "merged.reference.fa",
  showProgress = TRUE)

# Isoform switch analysis ----

mySwitchList<-isoformSwitchAnalysisPart1( #combined analysis function did not work, so just doing Part1
  switchAnalyzeRlist = mySwitchList,
  outputSequences = FALSE)

#Step by step annotation of mySwitchList using other functions in the package
mySwitchList.analyzed<-isoformSwitchTestDEXSeq(mySwitchList) #statistical comparison
mySwitchList.analyzed<-analyzeORF(mySwitchList.analyzed) #annotate ORFs; needed for plotting below 
mySwitchList.analyzed<-analyzeIntronRetention(mySwitchList.analyzed) #annotate alternative splicing
mySwitchList.func<-analyzeSwitchConsequences(mySwitchList.analyzed,
                                             consequencesToAnalyze=c('intron_retention',
                                                                     'ORF_seq_similarity',
                                                                     'NMD_status'))
#summarize number of switches found
extractSwitchSummary(mySwitchList.analyzed)
extractSwitchSummary(mySwitchList.func)

# Extract top switches ----
# extract the top n isoform switching events by q value
topSwitches<-extractTopSwitches(
  mySwitchList.analyzed, 
  inEachComparison = TRUE,
  filterForConsequences = FALSE, 
  n = 200, #take top 200 results for each comparison
  sortByQvals = TRUE) 
write_tsv(topSwitches,"top.200.switches.tsv")

topFuncSwitches<-extractTopSwitches( #same as above, but filtering for consequences AND q value
  mySwitchList.func, 
  inEachComparison = TRUE,
  filterForConsequences = TRUE, 
  n = 200, #take top 200 results for each comparison
  sortByQvals = TRUE) 
write_tsv(topFuncSwitches,"top.200.switchesByConsequence.tsv")

#extract all switches for all comparisons
topSwitches.all<-extractTopSwitches(
  mySwitchList.analyzed,
  inEachComparison = TRUE,
  n=Inf, #all switches
  sortByQvals = TRUE)
write_tsv(topSwitches.all,"all.switchesByQval.tsv")

# visualize switches by making a 'switch plot' ----
#switchPlot function is good for browsing multiple plots, but the isoform IDs get cut off
#For any results of interest, make the plots individually using the other switchPlot functions

#Make and save switch plot for BRCA2 isoform switch, normal vs. recurrent 
#This plot may have labels cut off; if it does, make the constituent plots individually (see below)
pdf("BRCA2.switchplot.pdf",
    width=17,height=8.5)
switchPlot( 
  mySwitchList.analyzed,
  gene='BRCA2',
  condition1 = 'normal',
  condition2 = 'recurrence',
  localTheme = theme_bw(base_size = 15))
dev.off()

#example of isoform usage plot alone (can also do other constituent plots alone)
switchPlotIsoUsage(mySwitchList.analyzed,
           gene="CDKN2A",
           condition1 = "normal",
           condition2 = "recurrence",
           localTheme = theme_bw())

# Interrogate which tumors express the NMD-sensitive vs. NMD-resistant isoforms of BRCA2 ----
#Extract the isoform abundance data from SwitchList
BRCA2_isoforms<-mySwitchList.func@.Data[[7]] %>% #the df with abundance by isoform for each sample
  dplyr::filter(isoform_id=="ENST00000380152.3" | isoform_id=="ENST00000544455.1") %>%
  pivot_longer(cols=2:67,names_to="sample",values_to="abundance") #pivot the df
BRCA2_isoforms<-BRCA2_isoforms[,c("sample","isoform_id","abundance")] #rearrange columns to be prettier

#make dataframes for just the samples expressing NMD-sensitive isoform
BRCA2_NMDsensitive<-BRCA2_isoforms %>% 
  filter(isoform_id=="ENST00000544455.1" & abundance>0)
#Repeat for samples expressing NMD_insensitive isoform
BRCA2_NMDinsensitive<-BRCA2_isoforms %>%
  filter(isoform_id=="ENST00000380152.3" & abundance>0)
#Annotate with relevant clinical info
BRCA2_NMDsensitive<-left_join(BRCA2_NMDsensitive,targets,by="sample")
BRCA2_NMDinsensitive<-left_join(BRCA2_NMDinsensitive,targets,by="sample")
