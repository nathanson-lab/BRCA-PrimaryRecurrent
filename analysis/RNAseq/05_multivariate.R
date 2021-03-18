# Objectives  -----------
# Assess sample relatedness by PCA plots and hierarchical clustering (dendrograms)
# Identify top altered genes for various comparisons, based on average between-group fold change
# Required inputs:
  # log2.cpm.filtered.norm for all tumors (generated in dataWrangling.R)
  # targets dataframe from TxImport_fromStringtie.R
# Outputs:
  # PCA plots stratified by any variable of interest (ex: ER status)
  # Dendrograms labeled by any grouping variable of interest (or sample names)
  # mydata.df = log2.cpm.filtered.norm df with new columns for between-group fold change differences

# Load packages ------
library(tidyverse) 
library(plotly) # for making interactive plots

# Identify variables of interest in study design file ----
targets
sample_type <- factor(targets$sample_type)
tissue_sample_type <- factor(targets$tissue_sample_type)
gene <- factor(targets$gene)
tissue_origin<-factor(targets$tissue_origin)

# Subset out log2.cpm.filtered.norm by tissue origin, so we can do PCA plots separately by tissue type (& together for whole cohort)----

#subset columns off the log2.cpm.filtered.norm.df by name --> make new dataframes by tissue origin
ovary.log2.cpm.filtered.norm.df<-select(log2.cpm.filtered.norm.df,one_of(ovary))
breast.log2.cpm.filtered.norm.df<-select(log2.cpm.filtered.norm.df,one_of(breast))
#coerce each df to matrix for use in hierarchical clustering, etc. Add geneID for rownames.
ovary.log2.cpm.filtered.norm<-data.matrix(ovary.log2.cpm.filtered.norm.df)
rownames(ovary.log2.cpm.filtered.norm)<-log2.cpm.filtered.norm.df$geneID
breast.log2.cpm.filtered.norm<-data.matrix(breast.log2.cpm.filtered.norm.df)
rownames(breast.log2.cpm.filtered.norm)<-log2.cpm.filtered.norm.df$geneID

# Hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame

#For whole cohort
distance <- dist(t(log2.cpm.filtered.norm), method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clusters, labels=sampleLabels)
plot(clusters, labels=targets$ER_status)
plot(clusters, labels=targets$TNBC)

#For ovary
distance <- dist(t(ovary.log2.cpm.filtered.norm), method = "maximum") 
clusters <- hclust(distance, method = "average") 
plot(clusters, labels=ovary)

#For breast
distance <- dist(t(breast.log2.cpm.filtered.norm), method = "maximum") 
clusters <- hclust(distance, method = "average") 
plot(clusters, labels=breast)
plot(clusters, labels=(targets$ER_status[targets$sample %in% breast]))
plot(clusters, labels=(targets$TNBC[targets$sample %in% breast]))

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
#look at the PCA result (pca.res) that you just created
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA
pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per


#PCA for ovary
ovary.pca.res <- prcomp(t(ovary.log2.cpm.filtered.norm), scale.=F, retx=T)
ovary.pc.var<-ovary.pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
ovary.pc.per<-round(ovary.pc.var/sum(ovary.pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
ovary.pc.per

#PCA for breast
breast.pca.res <- prcomp(t(breast.log2.cpm.filtered.norm), scale.=F, retx=T)
breast.pc.var<-breast.pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
breast.pc.per<-round(breast.pc.var/sum(breast.pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
breast.pc.per

# Visualize PCA results ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)

#Make a PCA plot of all samples - an overall plotting by sample type and gene 
pca.all<-ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels,color=tissue_sample_type,shape=gene) +
  geom_point(size=5,alpha=0.7) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot: All Samples",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  scale_shape_manual(name="Gene",values=c(16,17))+
  scale_color_discrete(name="Sample Type",labels=c("Normal Breast","Primary Breast Tumor","Recurrent Breast Tumor",
                                                  "Normal Fallopian tube","Primary Ovarian Tumor","Recurrent Ovarian Tumor"))
ggplotly(pca.all)

pca.all<-ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels,color=sample_type,shape=gene) +
  geom_point(size=5,alpha=0.7) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot: All Samples",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  scale_shape_manual(name="Gene",values=c(16,17))+
  scale_color_discrete(name="Sample Type",labels=c("Normal","Primary Tumor","Recurrent Tumor"))

#Try plotting other PCs against each other (PC2 vs. PC3, etc.) --> adapt code as needed for variable of interest
pca.all<-ggplot(pca.res.df) +
  aes(x=PC2, y=PC3, label=sampleLabels,color=factor(targets$TNBC)) +
  geom_point(size=5,alpha=0.7) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC2 (",pc.per[2],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[3],"%",")")) +
  labs(title="PCA plot: All Samples",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  scale_color_discrete(name="TNBC Status")
pca.all

#PCA plot color coded based on one variable (change out color and the labels for the legend to adapt)
pca.all<-ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels,color=targets$ER_status) +
  geom_point(size=5,alpha=0.7) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot: All Samples",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  scale_color_discrete(name="ER Status",labels=c("ER-","ER+","ER+ Ovary","N/A"))
pca.all


#Make a PCA plot of ovary samples 
ovary.pca.res.df <- as_tibble(ovary.pca.res$x)
pca.ovary<-ggplot(ovary.pca.res.df) +
  aes(x=PC1, y=PC2, label=ovary,color=tissue_sample_type[tissue_origin=="ovary"],shape=gene[tissue_origin=="ovary"]) +
  geom_point(size=5,alpha=0.7) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",ovary.pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",ovary.pc.per[2],"%",")")) +
  labs(title="PCA plot: Ovarian Samples",
       caption=paste0("produced on ", Sys.time())) +
  # coord_fixed() +
  theme_bw()+
  scale_shape_manual(name="Gene",values=c(16,17))+
  scale_color_discrete(name="Sample Type",labels=c("Normal Fallopian tube","Primary Ovarian Tumor","Recurrent Ovarian Tumor"))
pca.ovary

#PCA plot color coded based on one variable for ovary (change out color and the labels for the legend to adapt)
pca.ovary<-ggplot(ovary.pca.res.df) +
  aes(x=PC1, y=PC2, label=ovary,color=targets$sample_type[tissue_origin=="ovary"]) +
  geom_point(size=5,alpha=0.7) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",ovary.pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",ovary.pc.per[2],"%",")")) +
  labs(title="PCA plot: Ovarian Samples",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  scale_color_discrete(name="Sample Type")
pca.ovary

#PC2 vs. PC3, ovarian
pca.ovary<-ggplot(ovary.pca.res.df) +
  aes(x=PC2, y=PC3, label=ovary,color=targets$sample_type[tissue_origin=="ovary"]) +
  geom_point(size=5,alpha=0.7) +
  geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC2 (",ovary.pc.per[2],"%",")")) + 
  ylab(paste0("PC3 (",ovary.pc.per[3],"%",")")) +
  labs(title="PCA plot: Ovarian Samples",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  scale_color_discrete(name="Patient ID")
pca.ovary
ggplotly(pca.ovary)

#Make a PCA plot of breast samples
breast.pca.res.df <- as_tibble(breast.pca.res$x)
pca.breast<-ggplot(breast.pca.res.df) +
  aes(x=PC1, y=PC2, label=breast,color=tissue_sample_type[tissue_origin=="breast"],shape=gene[tissue_origin=="breast"]) +
  geom_point(size=5,alpha=0.7) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",breast.pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",breast.pc.per[2],"%",")")) +
  labs(title="PCA plot: Breast Samples",
       caption=paste0("produced on ", Sys.time())) +
  # coord_fixed() +
  theme_bw()+
  scale_shape_manual(name="Gene",values=c(16,17))+
  scale_color_discrete(name="Sample Type",labels=c("Normal Breast","Primary Breast Tumor","Recurrent Breast Tumor"))
pca.breast

#PCA plot color coded based on one variable for ovary (change out color and the labels for the legend to adapt)
pca.breast<-ggplot(breast.pca.res.df) +
  aes(x=PC1, y=PC2, label=breast.wsl,color=targets$ER_status[tissue_origin=="breast"]) +
  geom_point(size=5,alpha=0.7) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",breast.pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",breast.pc.per[2],"%",")")) +
  labs(title="PCA plot: Breast Samples",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  scale_color_discrete(name="ER Status",labels=c("ER-","ER+","N/A"))
pca.breast

#PC2 vs. PC3, breast
pca.breast<-ggplot(breast.pca.res.df) +
  aes(x=PC2, y=PC3, label=breast,color=targets$ER_status[tissue_origin=="breast"]) +
  geom_point(size=5,alpha=0.7) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC2 (",breast.pc.per[2],"%",")")) + 
  ylab(paste0("PC3 (",breast.pc.per[3],"%",")")) +
  labs(title="PCA plot: Breast Samples",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  scale_color_discrete(name="ER Status")
ggplotly(pca.breast)

# Add new fold change columns to log2.cpm.filtered.norm dataframes ----

#Make new columns for averages from each group of interest, for comparisons of interest
normal<-targets$sample[targets$sample_type=="normal"]
tumor<-targets$sample[targets$sample_type %in% c("primary","recurrence")]
breast.normal<-targets$sample[targets$tissue_sample_type=="breast_normal"]
ovary.normal<-targets$sample[targets$tissue_sample_type=="ovary_normal"]
breast.tumor<-targets$sample[targets$tissue_sample_type %in% c("breast_primary","breast_recurrence")]
ovary.tumor<-targets$sample[targets$tissue_sample_type %in% c("ovary_primary","ovary_recurrence")]
primary<-targets$sample[targets$sample_type=="primary"]
recurrence<-targets$sample[targets$sample_type=="recurrence"]
breast.primary<-targets$sample[targets$tissue_sample_type=="breast_primary"]
breast.recurrence<-targets$sample[targets$tissue_sample_type=="breast_recurrence"]
ovary.primary<-targets$sample[targets$tissue_sample_type=="ovary_primary"]
ovary.recurrence<-targets$sample[targets$tissue_sample_type=="ovary_recurrence"]

mydata.df <- log2.cpm.filtered.norm.df %>% 
  mutate(normal.AVG = rowMeans(.[normal]),
         tumor.AVG = rowMeans(.[tumor]),
         breast.normal.AVG = rowMeans(.[breast.normal]),
         ovary.normal.AVG = rowMeans(.[ovary.normal]),
         breast.tumor.AVG = rowMeans(.[breast.tumor]),
         ovary.tumor.AVG = rowMeans(.[ovary.tumor]),
         primary.AVG = rowMeans(.[primary]),
         recurrence.AVG = rowMeans(.[recurrence]),
         breast.primary.AVG = rowMeans(.[breast.primary]),
         breast.recurrence.AVG = rowMeans(.[breast.recurrence]),
         ovary.primary.AVG = rowMeans(.[ovary.primary]),
         ovary.recurrence.AVG = rowMeans(.[ovary.recurrence]),
         #now make columns comparing each of the averages above that you're interested in (for all tumors + by tumor type)
         LogFC.TumorVNormal = (tumor.AVG - normal.AVG), #up in tumor
         LogFC.PrimaryVNormal = (primary.AVG - normal.AVG), #up in primary
         LogFC.RecurrenceVNormal = (recurrence.AVG - normal.AVG), #up in recurrence
         LogFC.RecurrenceVPrimary = (recurrence.AVG - primary.AVG), #up in recurrence compared to primary
         LogFC.BreastTumorVNormal = (breast.tumor.AVG - breast.normal.AVG),
         LogFC.BreastPrimaryVNormal = (breast.primary.AVG - breast.normal.AVG),
         LogFC.BreastRecurrenceVNormal = (breast.recurrence.AVG - breast.normal.AVG),
         LogFC.BreastRecurrenceVBreastPrimary = (breast.recurrence.AVG - breast.primary.AVG),
         LogFC.OvaryTumorVNormal = (ovary.tumor.AVG - ovary.normal.AVG),
         LogFC.OvaryPrimaryVNormal = (ovary.primary.AVG - ovary.normal.AVG),
         LogFC.OvaryRecurrenceVNormal = (ovary.recurrence.AVG - ovary.normal.AVG),
         LogFC.OvaryRecurrenceVOvaryPrimary = (ovary.recurrence.AVG - ovary.primary.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

# Explore primary vs. recurrent gene expression in interactive scatterplot -----

myplot <- ggplot(mydata.df) + 
  aes(x=primary.AVG, y=recurrence.AVG, text = paste("Symbol:", geneID)) +
  geom_point(shape=16, size=1,alpha=0.2) +
  ggtitle("Gene Expression in Primary vs. Recurrent BRCA1/2 Breast and Ovarian Cancer") +
  xlab("Abundance in Primary Tumors \nLog2(counts per million)")+
  ylab("Abundance in Recurrences \nLog2(counts per million)")+
  theme_bw()

#now mouse over any points of interest to see what gene it is
ggplotly(myplot)


