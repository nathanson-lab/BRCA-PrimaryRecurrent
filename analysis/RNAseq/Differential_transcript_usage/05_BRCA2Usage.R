#Plot expression of BRCA2 ENST00000544455.1 and ENST00000380152.3 across cohort
#Assess relationship between isoform usage and sample type via CMH statistical test
#Required
  #BRCA2_isoforms dataframe (from DTU.R)
  #sample metadata (targets df from DTU.R)
#Output 
  #Heatmap depicting expression of BRCA2, grouped by sample type
  #Cochran-Mantel-Haenszel (CMH) test results
  #Heatmap depicting clinical characteristics of just those patients expressing ENST00000380152.3

#load packages ----
library(tidyverse)
library(reshape2)

#Make a heatmap of BRCA2 isoform usage across cohort -----
#Annotate BRCA2_isoforms df with sample metadata for plotting 
BRCA2_isoforms.toPlot<-BRCA2_isoforms
BRCA2_isoforms.toPlot<-left_join(BRCA2_isoforms.toPlot,targets,by="sample") %>% select(sample,isoform_id,abundance,wsl_id,
                                                                                       knm_id,sample_type,post_platinums,post_PARPi,
                                                                                       post_hormonal_chemo,post_nonhormonal_chemo,
                                                                                       post_radiation,loh_fromTumor,loh_fromPrimary)
#fix sample column from "_" to "-" format; fix sample_type column for prettier label
BRCA2_isoforms.toPlot <- BRCA2_isoforms.toPlot %>% mutate(sample=paste0(wsl_id,"-",knm_id))
BRCA2_isoforms.toPlot$sample_type[BRCA2_isoforms.toPlot$sample_type=="normal"]="Normal"
BRCA2_isoforms.toPlot$sample_type[BRCA2_isoforms.toPlot$sample_type=="primary"]="Primary"
BRCA2_isoforms.toPlot$sample_type[BRCA2_isoforms.toPlot$sample_type=="recurrence"]="Recurrence"
BRCA2_isoforms.toPlot$isoform_id<-factor(BRCA2_isoforms.toPlot$isoform_id,levels=c("ENST00000544455.1","ENST00000380152.3"))

#set 0's to NAs just for the purpose of not including them on our color scale for the plot
BRCA2_isoforms.toPlot$abundance[BRCA2_isoforms.toPlot$abundance==0]=NA
#arrange by sample type so recurrences, primaries, and normals will be grouped
BRCA2_isoforms.toPlot<-BRCA2_isoforms.toPlot %>% arrange(sample_type)
#order the isoforms so NMD-sensitive is first
BRCA2_isoforms.toPlot$sample<-factor(BRCA2_isoforms.toPlot$sample,levels=BRCA2_isoforms.toPlot$sample)

#Plot the heatmap of BRCA2 isoform usage by sample
ggplot(BRCA2_isoforms.toPlot,aes(x=isoform_id,y=sample,fill=abundance))+
  geom_tile(colour="black")+
  theme_bw()+
  xlab("")+
  ylab("Sample")+
  theme(axis.text.x.bottom=element_text(hjust=1,angle=45))+
  scale_fill_gradient(name="Expression",low="lightblue1",high="dodgerblue3",na.value = "white")+
  scale_x_discrete(expand=c(0,0),labels=c("ENST00000544455.1" = "ENST00000544455.1\nNMD-sensitive",
                                          "ENST00000380152.3" = "ENST00000380152.3\nNMD-insensitive"))+
  scale_y_discrete(expand=c(0,0))+
  facet_grid(sample_type ~.,scales="free_y",space="free_y")

#Do a Cochran-Mantel-Haenszel statistical test to compare presence of NMD-sensitive vs. NMD-insensitive BRCA2 by group----
#Build a df with count data for presence of each isoform by group
CMH.df<-data.frame(
  Group=c(rep("Normal",4),rep("Primary",4),rep("Recurrence",4)),
  Isoform=c(rep(c("NMD-sensitive","NMD-sensitive","NMD-insensitive","NMD-insensitive"),3)),
  Expressed=c(rep(c("yes","no"),6)),
  Count=c(12,0,0,12,16,4,4,16,16,18,20,14))

#make variables into factors to prevent re-ordering
CMH.df<-CMH.df %>% mutate(
  Group=factor(Group,levels=unique(Group)),
  Isoform=factor(Isoform,levels=unique(Isoform)),
  Expressed=factor(Expressed,levels=unique(Expressed)))

#cross-tabulate the data; view as flattened table
CMH.xtabs=xtabs(Count ~ Isoform + Expressed + Group, data=CMH.df)
ftable(CMH.xtabs)

#run the Mantel-Haenszel chi-squared test to get a p-value
mantelhaen.test(CMH.xtabs)

#Make a heatmap with clinical characteristics JUST for tumors expressing the NMD-insensitive isoform----
NMDins.cohort<-data.frame(
  express.ENST00000380152.3=TRUE,
  express.ENST00000544455.1=FALSE,
  BRCA2_NMDinsensitive[,4:20]
)
NMDins.cohort$sample<-paste0(NMDins.cohort$wsl_id,"-",NMDins.cohort$knm_id)
#Fill in two tumors that express both isoforms
NMDins.cohort$express.ENST00000544455.1[NMDins.cohort$sample=="7650-Brca1Br168"]<-TRUE
NMDins.cohort$express.ENST00000544455.1[NMDins.cohort$sample=="6489-Brca1Ov84"]<-TRUE
#Make logical (TRUE/FALSE) columns for plotting 
NMDins.cohort$post_platinums<-NMDins.cohort$post_platinums=="yes"
NMDins.cohort$post_PARPi<-NMDins.cohort$post_PARPi=="yes"
NMDins.cohort$post_nonhormonal_chemo<-NMDins.cohort$post_nonhormonal_chemo=="yes"
NMDins.cohort$post_hormonal_chemo<-NMDins.cohort$post_hormonal_chemo=="yes"
NMDins.cohort$post_radiation<-NMDins.cohort$post_radiation=="yes"
NMDins.cohort$loh_fromTumor<-NMDins.cohort$loh_fromTumor=="LOH"
NMDins.cohort$loh_fromPrimary<-NMDins.cohort$loh_fromPrimary=="LOH"
NMDins.cohort$recurrence<-NMDins.cohort$sample_type=="recurrence"

#melt df for plotting
NMDins.cohort<-NMDins.cohort %>% select(-c(wsl_id,knm_id,tissue_origin,tissue_sample_type,run_id,run_lane,prepped_by,site,sample_type)) #get rid of extraneous columns
melted<-melt(NMDins.cohort,id.vars = c("sample","gene"),measure.vars = c("recurrence","express.ENST00000380152.3","express.ENST00000544455.1","post_platinums","post_PARPi",
                                                                  "post_nonhormonal_chemo","post_hormonal_chemo","post_radiation","loh_fromTumor","loh_fromPrimary"))
colors<-c("TRUE"="lightgreen","FALSE"="lightcoral")

#Plot the heatmap of clinical characteristics and BRCA2 isoform usage
ggplot(melted,aes(x=variable,y=sample,fill=value))+
  geom_tile(colour="black")+
  theme_bw()+
  xlab("")+
  ylab("Tumor")+
  theme(axis.text.x.bottom=element_text(hjust=1,angle=45))+
  scale_fill_manual(name=NULL,values=colors,breaks=c("TRUE","FALSE"),labels=c("Yes","No"))+
  scale_x_discrete(expand=c(0,0),labels=c("Recurrence","Express ENST00000380152.3\n(NMD-insensitive)","Express ENST00000544455.1\n(NMD-sensitive)",
                                          "Received Platinums","Received PARPi","Received Nonhormonal Chemo","Received Hormonal Chemo","Received Radiation",
                                          "LOH in Tumor","LOH in Patient's Primary Tumor"))+
  facet_grid(gene ~ .,scales = "free_y",space="free_y")
