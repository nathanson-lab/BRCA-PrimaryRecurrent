# Objectives ----
  # Identify differentially expressed genes (DEGs) for comparisons of interest
# Required input
  # 3 filtered, normalized DGELists (all-tumor, breast, and ovarian) from dataWrangling.R
  # targets df from TxImport_fromStringtie.R
# Outputs
  # a df for each comparison of interest (ex:RecurrenceVPrimary.df) 
    #containing Bayesian statistics for each gene based on limma modeling mean-variance
  # volcano plots for each comparison (adj. p value and logFC by gene)
    #can be made interactive for data exploration
  # diffGenes dataframes for all-tumor, breast-only, and ovarian-only comparisons
    #containing only the genes with adj.p<0.05 and |logFC|>1
    #depicted by Venn diagrams

# Load packages -----
library(tidyverse) 
library(limma) 
library(edgeR)
library(gt) 
library(DT) 
library(plotly) 

# Set up your design matrices ----

#Set up design matrix based on sample_type column (comparisons for primary/recurrent/normal, entire cohort)
group <- factor(targets$sample_type)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

#Set up design matrix for breast-only comparisons
breast.group <-factor(targets$tissue_sample_type[targets$sample %in% breast])
breast.design <- model.matrix(~0 + breast.group)
colnames(breast.design) <- levels(breast.group)

#Set up design matrix for ovarian-only comparisons
ovary.group <- factor(targets$tissue_sample_type[targets$sample %in% ovary])
ovary.design <- model.matrix(~0 + ovary.group)
colnames(ovary.design) <- levels(ovary.group)

# NOTE: if you need a paired analysis (a.k.a.'blocking' design) or have a batch effect, the following design is useful
# design <- model.matrix(~block + treatment) 
# this is just an example. 'block' and 'treatment' would need to be objects in your environment

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)

#Model mean-variance trend and fit linear model for breast-only and ovarian-only data
v.DEGList.filtered.norm.breast <- voom(breast.DGEList.filtered.norm, breast.design, plot = FALSE)
breast.fit <- lmFit(v.DEGList.filtered.norm.breast, breast.design)
v.DEGList.filtered.norm.ovary <- voom(ovary.DGEList.filtered.norm, ovary.design, plot = FALSE)
ovary.fit <- lmFit(v.DEGList.filtered.norm.ovary, ovary.design)

# Contrast matrix ----
#For sample_type comparisons, entire cohort
contrast.matrix <- makeContrasts(primary - normal,
                                 recurrence - normal,
                                 recurrence - primary,
                                 levels=design)
#breast-only comparisons
breast.contrast.matrix <-makeContrasts(breast_primary - breast_normal,
                                       breast_recurrence - breast_normal,
                                       breast_recurrence - breast_primary,
                                       levels=breast.design)
#ovarian-only comparisons
ovary.contrast.matrix <- makeContrasts(ovary_primary - ovary_normal,
                                       ovary_recurrence - ovary_normal,
                                       ovary_recurrence - ovary_primary,
                                       levels=ovary.design)

# extract the linear model fit -----

#For sample_type contrasts, entire cohort 
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
#write.fit(ebFit, file="lmfit_results.txt")

#For breast-only contrasts
breast.fits <- contrasts.fit(breast.fit, breast.contrast.matrix)
breast.ebFit <- eBayes(breast.fits)

#For ovarian-only contrasts
ovary.fits <- contrasts.fit(ovary.fit, ovary.contrast.matrix)
ovary.ebFit <- eBayes(ovary.fits)

# TopTable to view DEGs -----
#For sample_type contrast matrix, entire cohort
PrimaryVNormal.TopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
RecurrenceVNormal.TopHits <- topTable(ebFit, adjust ="BH", coef=2, number=40000, sort.by="logFC")
RecurrenceVPrimary.TopHits <- topTable(ebFit, adjust ="BH", coef=3, number=40000, sort.by="logFC")

#For breast-only contrast matrix
BreastPrimaryVNormal.TopHits <- topTable(breast.ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
BreastRecurrenceVNormal.TopHits <- topTable(breast.ebFit, adjust ="BH", coef=2, number=40000, sort.by="logFC")
BreastRecurrenceVPrimary.TopHits <- topTable(breast.ebFit, adjust ="BH", coef=3, number=40000, sort.by="logFC")

#For ovarian-only contrast matrix
OvaryPrimaryVNormal.TopHits <- topTable(ovary.ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
OvaryRecurrenceVNormal.TopHits <- topTable(ovary.ebFit, adjust ="BH", coef=2, number=40000, sort.by="logFC")
OvaryRecurrenceVPrimary.TopHits <- topTable(ovary.ebFit, adjust ="BH", coef=3, number=40000, sort.by="logFC")

# convert each df to a tibble
PrimaryVNormal.df <- PrimaryVNormal.TopHits %>%
  as_tibble(rownames = "geneID")
RecurrenceVNormal.df <- RecurrenceVNormal.TopHits %>%
  as_tibble(rownames = "geneID")
RecurrenceVPrimary.df <- RecurrenceVPrimary.TopHits %>%
  as_tibble(rownames = "geneID")
BreastPrimaryVNormal.df <- BreastPrimaryVNormal.TopHits %>%
  as_tibble(rownames = "geneID")
BreastRecurrenceVNormal.df <-BreastRecurrenceVNormal.TopHits %>%
  as_tibble(rownames = "geneID")
BreastRecurrenceVPrimary.df <-BreastRecurrenceVPrimary.TopHits %>%
  as_tibble(rownames = "geneID")
OvaryPrimaryVNormal.df <-OvaryPrimaryVNormal.TopHits %>%
  as_tibble(rownames = "geneID")
OvaryRecurrenceVNormal.df <-OvaryRecurrenceVNormal.TopHits %>%
  as_tibble(rownames = "geneID")
OvaryRecurrenceVPrimary.df <-OvaryRecurrenceVPrimary.TopHits %>%
  as_tibble(rownames = "geneID")

# Volcano Plots ----
# in topTable function above, set 'number=40000' to capture all genes

#To label specific genes : add annotation for the point representing the gene, 
  #then play around with positioning of the label to obscure the data as little as possible

# make a volcano plot for each comparison
vplot <- ggplot(PrimaryVNormal.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("point",x=PrimaryVNormal.df %>% filter(geneID=="PARP1") %>% pull(logFC),
           y=-log10(PrimaryVNormal.df %>% filter(geneID=="PARP1") %>% pull(adj.P.Val)),
           color="red")+
  annotate("label",x=1.7,y=7.2,label="PARP1",color="red")+
  labs(title="Differential Gene Expression in BRCA1/2 Primary Breast and Ovarian Primary Tumors",
       subtitle = "Compared to Normal Breast and Fallopian tube from BRCA1/2 Mutation Carriers",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  xlab("Log2(Fold change)")+
  scale_x_continuous(breaks=seq(-12,10,by=2))+
  scale_y_continuous(breaks=seq(0,16,by=2))

vplot <- ggplot(RecurrenceVNormal.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("point",x=RecurrenceVNormal.df %>% filter(geneID=="MYC") %>% pull(logFC),
           y=-log10(RecurrenceVNormal.df %>% filter(geneID=="MYC") %>% pull(adj.P.Val)),
           color="red")+
  annotate("label",x=-1,y=.7,label="MYC",color="red")+
  annotate("point",x=RecurrenceVNormal.df %>% filter(geneID=="MYCN") %>% pull(logFC),
           y=-log10(RecurrenceVNormal.df %>% filter(geneID=="MYCN") %>% pull(adj.P.Val)),
           color="red")+
  annotate("label",x=2.3,y=2.2,label="MYCN",color="red")+
  annotate("point",x=RecurrenceVNormal.df %>% filter(geneID=="MYCL") %>% pull(logFC),
           y=-log10(RecurrenceVNormal.df %>% filter(geneID=="MYCL") %>% pull(adj.P.Val)),
           color="red")+
  annotate("label",x=4.1,y=2.2,label="MYCL",color="red")+
  labs(title="Differential Gene Expression in BRCA1/2 Breast and Ovarian Tumor Recurrences",
       subtitle = "Compared to Normal Breast and Fallopian tube from BRCA1/2 Mutation Carriers",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  xlab("Log2(Fold change)")+
  scale_x_continuous(breaks=seq(-12,10,by=2))+
  scale_y_continuous(breaks=seq(0,24,by=2))

vplot <- ggplot(RecurrenceVPrimary.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  labs(title="Differential Gene Expression in BRCA1/2 Breast and Ovarian Tumor Recurrences",
       subtitle = "Compared to BRCA1/2 Primary Breast and Ovarian Tumors",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  xlab("Log2(Fold change)")+
  scale_x_continuous(breaks=seq(-12,10,by=2))+
  scale_y_continuous(breaks=seq(0,24,by=2))

vplot <- ggplot(BreastPrimaryVNormal.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 13, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -14, ymin = -log10(0.01), ymax = 13, alpha=.2, fill="#2C467A") +
  labs(title="Differential Gene Expression in BRCA1/2 Breast Primary Tumors",
       subtitle = "Compared to Normal Breast from BRCA1/2 Mutation Carrier",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  xlab("Log2(Fold change)")+
  scale_x_continuous(breaks=seq(-12,10,by=2))+
  scale_y_continuous(breaks=seq(0,12,by=2))

vplot <- ggplot(BreastRecurrenceVNormal.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 13, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -14, ymin = -log10(0.01), ymax = 13, alpha=.2, fill="#2C467A") +
  labs(title="Differential Gene Expression in BRCA1/2 Breast Tumor Recurrences",
       subtitle = "Compared to Normal Breast from BRCA1/2 Mutation Carrier",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  xlab("Log2(Fold change)")+
  scale_x_continuous(breaks=seq(-12,10,by=2))+
  scale_y_continuous(breaks=seq(0,12,by=2))

vplot <- ggplot(BreastRecurrenceVPrimary.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 13, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -14, ymin = -log10(0.01), ymax = 13, alpha=.2, fill="#2C467A") +
  labs(title="Differential Gene Expression in BRCA1/2 Breast Tumor Recurrences",
       subtitle = "Compared to BRCA1/2 Primary Breast Tumors",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  xlab("Log2(Fold change)")+
  scale_x_continuous(breaks=seq(-12,10,by=2))+
  scale_y_continuous(breaks=seq(0,12,by=2))

vplot <- ggplot(OvaryPrimaryVNormal.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 13, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -14, ymin = -log10(0.01), ymax = 13, alpha=.2, fill="#2C467A") +
  labs(title="Differential Gene Expression in BRCA1/2 Ovarian Primary Tumors",
       subtitle = "Compared to Normal Fallopian Tube from BRCA1/2 Mutation Carrier",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  xlab("Log2(Fold change)")+
  scale_x_continuous(breaks=seq(-12,10,by=2))+
  scale_y_continuous(breaks=seq(0,14,by=2))

vplot <- ggplot(OvaryRecurrenceVNormal.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 13, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -14, ymin = -log10(0.01), ymax = 13, alpha=.2, fill="#2C467A") +
  labs(title="Differential Gene Expression in BRCA1/2 Ovarian Tumor Recurrences",
       subtitle = "Compared to Normal Fallopian Tube from BRCA1/2 Mutation Carrier",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  xlab("Log2(Fold change)")+
  scale_x_continuous(breaks=seq(-14,10,by=2))+
  scale_y_continuous(breaks=seq(0,20,by=2))

vplot <- ggplot(OvaryRecurrenceVPrimary.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 13, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -14, ymin = -log10(0.01), ymax = 13, alpha=.2, fill="#2C467A") +
  labs(title="Differential Gene Expression in BRCA1/2 Ovarian Tumor Recurrences",
       subtitle = "Compared to BRCA1/2 Primary Ovarian Tumors",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  xlab("Log2(Fold change)")+
  scale_x_continuous(breaks=seq(-12,10,by=2))+
  scale_y_continuous(breaks=seq(0,14,by=2))

# Make any of these volcano plots above interactive with plotly
ggplotly(vplot)

# decideTests to pull out the DEGs and make Venn Diagram ----

#Pull out only the genes that meet our pvalue and log fold change thresholds
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)
breast.results <- decideTests(breast.ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)
ovary.results <- decideTests(ovary.ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)

# take a look at what the results of decideTests looks like
vennDiagram(results, include="up") 
vennDiagram(results, include="down")
vennDiagram(breast.results, include="up")
vennDiagram(breast.results, include="down")
vennDiagram(ovary.results, include="up")
vennDiagram(ovary.results, include="down")

# retrieve expression data for your DEGs ----
#Add tumor names to your expression data
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
colnames(v.DEGList.filtered.norm.breast$E) <- breast
colnames(v.DEGList.filtered.norm.ovary$E) <- ovary

#pull out diffGenes: keeping all genes for which >=1 comparison column has a non-zero value in results
diffGenes<- v.DEGList.filtered.norm$E[((results[,1]|results[,2]|results[,3])!=0),]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
  
#Repeat for breast-only and ovarian-only diffGenes
diffGenes.breast <- v.DEGList.filtered.norm.breast$E[((breast.results[,1]|breast.results[,2]|breast.results[,3])!=0),]
diffGenes.breast.df <- as_tibble(diffGenes.breast, rownames = "geneID")
diffGenes.ovary <- v.DEGList.filtered.norm.ovary$E[((ovary.results[,1]|ovary.results[,2]|ovary.results[,3])!=0),]
diffGenes.ovary.df <- as_tibble(diffGenes.ovary, rownames = "geneID")

#write DEGs to a file
write_tsv(diffGenes.df,"DiffGenes.all.txt") #NOTE: this .txt file can be directly used for input into other clustering or network analysis tools (e.g., String, Clust (https://github.com/BaselAbujamous/clust, etc.)
write_tsv(diffGenes.breast.df,"DiffGenes.breast.txt")
write_tsv(diffGenes.ovary.df,"DiffGenes.ovary.txt")