# Objectives -----------
# Identify modules of differentially expressed genes that are co-expressed
# Required inputs:
  #3 diffGenes matrices generated in diffGenes.R 
  #targets dataframe from TxImport_fromStringtie.R
# Outputs:
  # heatmaps depicting gene modules by color (hclustering of genes and samples)
  # gene modules for downstream analysis (for example, module of genes down in tumor)
    # number of genes per module depicted in bar graphs

# Load packages -----
library(tidyverse)
library(limma) #we only use limma in this script for the 'avearrays' function
library(RColorBrewer) #need colors to make heatmaps
library(gplots) #the heatmap2 function in this package is a primary tool for making heatmaps

# Choose your colorblind-friendly palette ----
myheatcolors3 <- c("#fed976", "#268f9c")

# Cluster DEGs ----

#subset diffGenes to tumor only matrices, by tumor type
diffGenes.breastTumor<-diffGenes.breast[, 1:20] #first 20 columns are tumor
diffGenes.ovaryTumor<-diffGenes.ovary[,1:34] #first 34 columns are tumor

#Make heatmaps out of diffGenes matrices 
#begin by clustering the genes (rows) in each set of differentially expressed genes
# we use the 'cor' function and the pearson method for finding all pairwise correlations of genes
# '1-cor' converts this to a 0-2 scale for each of these correlations, which can then be used to calculate a distance matrix using 'as.dist'
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") 

#now cluster samples (columns)
#we may not acutally use this clustering result, but it's good to have just in case
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete") #cluster columns by spearman correlation
#note: we use Spearman, instead of Pearson, for clustering samples because it gives equal weight to highly vs lowly expressed transcripts or genes

#Repeat for breast (all) and for breast tumor
clustRows.breast <- hclust(as.dist(1-cor(t(diffGenes.breast), method="pearson")), method="complete") 
clustColumns.breast <- hclust(as.dist(1-cor(diffGenes.breast, method="spearman")), method="complete") 
clustRows.breastTumor <- hclust(as.dist(1-cor(t(diffGenes.breastTumor), method="pearson")), method="complete") 
clustColumns.breastTumor <- hclust(as.dist(1-cor(diffGenes.breastTumor, method="spearman")), method="complete") 

#Repeat for ovary (all) and for ovary tumor 
clustRows.ovary <- hclust(as.dist(1-cor(t(diffGenes.ovary), method="pearson")), method="complete") 
clustColumns.ovary <- hclust(as.dist(1-cor(diffGenes.ovary, method="spearman")), method="complete")
clustRows.ovaryTumor <- hclust(as.dist(1-cor(t(diffGenes.ovaryTumor), method="pearson")), method="complete") 
clustColumns.ovaryTumor <- hclust(as.dist(1-cor(diffGenes.ovaryTumor, method="spearman")), method="complete")

#Cut the resulting tree and create color vector for clusters.  Can vary cut height as desired.
module.assign <- cutree(clustRows, k=2)
module.assign.breast <- cutree(clustRows.breast, k=2)
module.assign.breastTumor <- cutree(clustRows.breastTumor, k=2)
module.assign.ovary <- cutree(clustRows.ovary, k=2)
module.assign.ovaryTumor <- cutree(clustRows.ovaryTumor, k=2)

#now assign a color to each module (makes it easy to identify and manipulate)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
module.color.breast <- rainbow(length(unique(module.assign.breast)), start=0.1, end=0.9) 
module.color.breast <-module.color.breast[as.vector(module.assign.breast)]
module.color.breastTumor <- rainbow(length(unique(module.assign.breastTumor)), start=0.1, end=0.9) 
module.color.breastTumor <-module.color.breastTumor[as.vector(module.assign.breastTumor)]
module.color.ovary <- rainbow(length(unique(module.assign.ovary)), start=0.1, end=0.9) 
module.color.ovary <- module.color.ovary[as.vector(module.assign.ovary)]
module.color.ovaryTumor <- rainbow(length(unique(module.assign.ovaryTumor)), start=0.1, end=0.9) 
module.color.ovaryTumor <- module.color.ovaryTumor[as.vector(module.assign.ovaryTumor)]

# Produce a static heatmap of DEGs ----
#plot the hclust results as a heatmap

#all samples
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(10,26)) 

#breast only 
heatmap.2(diffGenes.breast,
          Rowv=as.dendrogram(clustRows.breast), 
          Colv=as.dendrogram(clustColumns.breast),
          RowSideColors=module.color.breast,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(10,26)) 

#breast tumor only
heatmap.2(diffGenes.breastTumor,
          Rowv=as.dendrogram(clustRows.breastTumor), 
          Colv=as.dendrogram(clustColumns.breastTumor),
          RowSideColors=module.color.breast,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(10,26)) 

#ovary only 
heatmap.2(diffGenes.ovary,
          Rowv=as.dendrogram(clustRows.ovary), 
          Colv=as.dendrogram(clustColumns.ovary),
          RowSideColors=module.color.ovary,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(10,26)) 

#ovary tumor only
heatmap.2(diffGenes.ovaryTumor,
          Rowv=as.dendrogram(clustRows.ovaryTumor), 
          Colv=as.dendrogram(clustColumns.ovaryTumor),
          RowSideColors=module.color.ovary,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(10,26)) 

# Simplify heatmaps into average ----

#a useful way to simplify heatmaps, especially when there are many conditions, is to average your biological replicates and display only one column per condition
#rerun the heatmap script above using diffGenes.AVG as input instead of diffGenes
colnames(diffGenes) <- targets$sample_type
#function from the limma package to average your replicates 
diffGenes.AVG <- avearrays(diffGenes)
heatmap.2(diffGenes.AVG, 
          Rowv=as.dendrogram(clustRows), 
          RowSideColors=module.color,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(10,26)) 

#breast averages
colnames(diffGenes.breast) <- targets$sample_type[tissue_origin=="breast"]
diffGenes.breast.AVG <- avearrays(diffGenes.breast)
heatmap.2(diffGenes.breast.AVG, 
          Rowv=as.dendrogram(clustRows.breast), 
          RowSideColors=module.color.breast,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(10,26)) 

#breast tumor averages
colnames(diffGenes.breastTumor) <- targets$sample_type[targets$tissue_sample_type %in% c("breast_primary","breast_recurrence")]
diffGenes.breastTumor.AVG <- avearrays(diffGenes.breastTumor)
heatmap.2(diffGenes.breastTumor.AVG, 
          Rowv=as.dendrogram(clustRows.breastTumor), 
          RowSideColors=module.color.breastTumor,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(10,26)) 

#ovary averages
colnames(diffGenes.ovary) <- targets$sample_type[tissue_origin=="ovary"]
diffGenes.ovary.AVG <- avearrays(diffGenes.ovary)
heatmap.2(diffGenes.ovary.AVG, 
          Rowv=as.dendrogram(clustRows.ovary), 
          RowSideColors=module.color.ovary,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(10,26)) 

#ovary tumor averages
colnames(diffGenes.ovaryTumor) <- targets$sample_type[targets$tissue_sample_type %in% c("ovary_primary","ovary_recurrence")]
diffGenes.ovaryTumor.AVG <- avearrays(diffGenes.ovaryTumor)
heatmap.2(diffGenes.ovaryTumor.AVG, 
          Rowv=as.dendrogram(clustRows.ovaryTumor), 
          RowSideColors=module.color.ovaryTumor,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(10,26)) 

# View modules of co-regulated genes ----
# view your color assignments for the different clusters
#won't look at breast tumor and ovarian tumor modules any further because they aren't showing between-group differences
names(module.color) <- names(module.assign) 
names(module.color.breast)<-names(module.assign.breast)
names(module.color.ovary)<-names(module.assign.ovary)

module.assign.df <- as_tibble(as.list(module.assign))
module.assign.breast.df<-as_tibble(as.list(module.assign.breast))
module.assign.ovary.df<-as_tibble(as.list(module.assign.ovary))

module.assign.pivot <- pivot_longer(module.assign.df, # dataframe to be pivoted
                          cols = 1:ncol(module.assign.df), # column names to be stored as a SINGLE variable
                          names_to = "geneID", # name of that new variable (column)
                          values_to = "module") # name of new variable (column) storing all the values (data)
module.assign.breast.pivot <- pivot_longer(module.assign.breast.df,
                                           cols = 1:ncol(module.assign.breast.df), 
                                           names_to = "geneID",
                                           values_to = "module") 
module.assign.ovary.pivot <- pivot_longer(module.assign.ovary.df,
                                          cols = 1:ncol(module.assign.ovary.df), 
                                          names_to = "geneID",
                                          values_to = "module")

#Make a bargraph depicting number of genes in each module
module.assign.pivot <- module.assign.pivot %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "#FF9900",
    module == 2 ~ "#FF0099"))
ggplot(module.assign.pivot) +
  aes(factor(module)) +
  geom_bar(aes(fill=moduleColor)) +
  scale_fill_discrete(name="Module",labels=c("Up in Tumor","Down in Tumor"))+
  ylab("Number of Genes")+
  xlab("Module Number")+
  theme_bw()

module.assign.breast.pivot <- module.assign.breast.pivot %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "#FF9900",
    module == 2 ~ "#FF0099"))
ggplot(module.assign.breast.pivot) +
  aes(factor(module)) +
  geom_bar(aes(fill=moduleColor)) +
  scale_fill_discrete(name="Module",labels=c("Up in Breast Tumor","Down in Breast Tumor"))+
  ylab("Number of Genes")+
  xlab("Module Number")+
  theme_bw()

module.assign.ovary.pivot <- module.assign.ovary.pivot %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "#FF9900",
    module == 2 ~ "#FF0099"))
ggplot(module.assign.ovary.pivot) +
  aes(factor(module)) +
  geom_bar(aes(fill=moduleColor)) +
  scale_fill_discrete(name="Module",labels=c("Up in Ovarian Tumor","Down in Ovarian Tumor"))+
  ylab("Number of Genes")+
  xlab("Module Number")+
theme_bw()

# Export modules for downstream analysis ----
#choose a cluster(s) of interest by selecting the corresponding number based on the previous graph
#IMPORTANT: edit line below to go between modules for writing to file, etc.
modulePick <- 1 
#now we pull out the genes from this module using a fancy subsetting operation on a named vector
myModule <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub <- hclust(as.dist(1-cor(t(myModule), method="pearson")), method="complete") 

#prints out genes in the order you see them in the cluster
moduleSymbols <- tibble(geneID = rev(hrsub$labels[hrsub$order]))
moduleData <- diffGenes[moduleSymbols$geneID,]
moduleData.df <- as_tibble(moduleData, rownames = "geneSymbol")
write_tsv(moduleData.df,"all.module.DownInTumor.tsv")

#Repeat for breast modules
modulePick <- 2 
myModule <- diffGenes.breast[names(module.assign.breast[module.assign.breast %in% modulePick]),] 
hrsub <- hclust(as.dist(1-cor(t(myModule), method="pearson")), method="complete") 
moduleSymbols <- tibble(geneID = rev(hrsub$labels[hrsub$order]))
moduleData <- diffGenes.breast[moduleSymbols$geneID,]
moduleData.df <- as_tibble(moduleData, rownames = "geneSymbol")
write_tsv(moduleData.df,"all.module.UpInBreastTumor.tsv")

#Repeat for ovarian samples
modulePick <- 2 
myModule <- diffGenes.ovary[names(module.assign.ovary[module.assign.ovary %in% modulePick]),] 
hrsub <- hclust(as.dist(1-cor(t(myModule), method="pearson")), method="complete") 
moduleSymbols <- tibble(geneID = rev(hrsub$labels[hrsub$order]))
moduleData <- diffGenes.ovary[moduleSymbols$geneID,]
moduleData.df <- as_tibble(moduleData, rownames = "geneSymbol")
write_tsv(moduleData.df,"all.module.DownInOvaryTumor.tsv")
