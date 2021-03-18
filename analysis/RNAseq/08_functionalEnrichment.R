# Objectives ----
#Run Gene Ontology (GO) analysis on 
  #TopHits df (from diffGenes.R)
  #Individual gene modules (from modules.R)
#Run GSEA based on logFC columns in the mydata.df dataframe (from multivariate.R)
#Required inputs:
  # mydata.df from multivariate.R
  # targets dataframe from TxImport_fromStringtie.R
  # modules from modules.R
  # Here we pull gene sets from msigdbr, but can load from elsewhere
#Outputs
  # GSEA bubble plots
  # GSEA results in tables (interactive and not)
  # GSEA Enrichment plots for individual gene set(s) of interest
  # GO Manhattan plots (interactive and not)
  # GO results written out to file

# Load packages ----
library(tidyverse)
library(limma)
library(DT) #interactive table for browsing GSEA results
library(gt) #for making tables
library(gplots) #for heatmaps
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots

# Carry out GO enrichment using gProfiler2 ----
# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis 

#Run gost on any of the TopHits dfs created in diffGenes.R script
#run gost on the top x genes by magnitude of logFC (so this would capture positive or negative) 
gost.res <- gost(rownames(OvaryRecurrenceVPrimary.TopHits[1:1000,]), 
                 organism = "hsapiens", correction_method = "fdr")

# produce an interactive manhattan plot of enriched GO terms
gostplot(gost.res, interactive = T, capped = T)

# produce a static manhattan plot of enriched GO terms (same as above but with highlighted terms)
mygostplot<-gostplot(gost.res, interactive = F, capped = T) 
publish_gostplot(
  mygostplot, 
  highlight_terms = c("GO:0002376","GO:0006955","GO:0002682","GO:0050900","GO:0002449","GO:0050776",
                      "GO:0002250","GO:0002443","GO:0002684","GO:0045321","GO:0002440","GO:0050778",
                      "MIRNA:hsa-miR-335-5p"),
  filename = NULL,
  width = NA,
  height = NA)

#Run gost on gene modules from modules.R
  #Go back to modules.R and edit "myModule" to point to module of choice in comparison of choice
gost.res <- gost(rownames(myModule), 
                 organism = "hsapiens", correction_method = "fdr")
write_tsv(gost.res$result,"plots/Gost_plots/all.DownInTumor.mod1.tsv")
mygostplot<-gostplot(gost.res, interactive = F, capped = F) 
publish_gostplot(
  mygostplot, 
  highlight_terms = c("GO:0000278","GO:0007059","GO:0000280","GO:000070","GO:0022402","GO:0007049",
                      "GO:0098813","GO:0140014","GO:0044772","GO:0042119","GO:0036230","GO:0043312",
                      "GO:0034097","GO:0002275","GO:0002283","GO:0002274","GO:0002366","GO:0006261",
                      "GO:0006260","REAC:R-HSA-1640170","REAC:R-HSA-69278","MIRNA:hsa-miR-193b-3p"),
  filename = NULL,
  width = NA,
  height = NA)

# Perform GSEA using clusterProfiler ----

# choose specific msigdb collection/subcollections from msigdbr
  #only want the gene set and gene symbol columns
#Hallmark gene sets
hs_gsea_h <-msigdbr(species = "Homo sapiens",
                    category = "H") %>%
  dplyr::select(gs_name,gene_symbol)
#C6 (Oncogenic) gene sets
hs_gsea_c6 <- msigdbr(species = "Homo sapiens",
                      category ="C6") %>%
  dplyr::select(gs_name, gene_symbol)
#C7 (Immune) gene sets
hs_gsea_c7 <-msigdbr(species = "Homo sapiens",
                     category = "C7") %>%
  dplyr::select(gs_name,gene_symbol)

# Prepare data for GSEA
# Using mydata.df, pull out columns with gene symbols and LogFC for at least one pairwise comparison

##workflow for ONE GSEA run on ONE df from above (see below to loop through several runs at once)
# construct a named vector
mydata.gsea <- mydata.df.TumorVNormal$LogFC.TumorVNormal
names(mydata.gsea) <- as.character(mydata.df.TumorVNormal$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)
# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=TRUE) #change TERM2GENE arg to your collection of interest
myGSEA.df <- as_tibble(myGSEA.res@result)

##workflow to run several GSEA at once (multiple collections + all LogFC comparison columns)
#will run GSEA using hallmark, c6, and c7 gene sets; for each LogFC comparison in mydata.df
#results are appropriately-labeled GSEA objects
for (i in (grep("LogFC",names(mydata.df),value=TRUE))){
  colIndex<-match(i,names(mydata.df)) #save the index of your LogFC column
  temp.gsea<-pull(mydata.df[colIndex]) #save logfc values as a vector
  names(temp.gsea)<-as.character(mydata.df$geneID) #store geneIDs as names for the named vector
  temp.gsea<-sort(temp.gsea,decreasing=TRUE) #sort the named vector in prep for GSEA
  #run GSEA with hallmark gene sets (h), c6 gene set (oncogenic), and c7 gene set (immune)
  assign(paste0(i,".h.GSEA.res"),GSEA(temp.gsea, TERM2GENE=hs_gsea_h, verbose=TRUE))
  assign(paste0(i,".c6.GSEA.res"),GSEA(temp.gsea, TERM2GENE=hs_gsea_c6, verbose=TRUE))
  assign(paste0(i,".c7.GSEA.res"),GSEA(temp.gsea, TERM2GENE=hs_gsea_c7, verbose=TRUE))
  rm(temp.gsea)
}

# view GSEA results as an interactive table (input results need to be a df) --> use this to browse
datatable(as_tibble(LogFC.RecurrenceVPrimary.h.GSEA.res), 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Hallmark Gene Sets enriched in BRCA1/2 Tumors Compared to Normals',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(4:8), digits=5)

#make and save gt tables of results for hallmark and c6 gene sets; output for c7 is too large for table to work
LogFC.OvaryTumorVNormal.h.GSEA.res %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  select(-c("Description","enrichmentScore","pvalue","core_enrichment")) %>%
  gt() %>%
  fmt_number(columns=3,decimals=2) %>%
  fmt_scientific(columns=c(4:5),decimals=2) %>%
  tab_header(title = md("**Gene Set Enrichment Analysis Results: BRCA1/2 Primary+Recurrent Ovarian Tumors vs. BRCA1/2 Normal Fallopian Tube**"),
             subtitle = md("*Hallmark Gene Sets*")) %>%
  gtsave("tables/GSEA/OvaryTumorVNormal.GSEA.h.pdf")

# create enrichment plots using the enrichplot package (input needs to be a GSEA result object)
gseaplot2(LogFC.RecurrenceVPrimary.h.GSEA.res, 
          geneSetID = c(1), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          title = "Hallmark Myogenesis Enrichment in BRCA1/2 Recurrences Compared to Primary Tumors (adj.p=0.00001)") #or can turn off

# add a variable to this result that matches enrichment direction with phenotype
#change to GSEA result of interest; also change labels to reflect conditions
myGSEA.df <- as_tibble(LogFC.RecurrenceVPrimary.h.GSEA.res)%>%
  arrange(p.adjust) %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "Enriched in Recurrence",
    NES < 0 ~ "Enriched in Primary"))

#no gene sets enriched in primary, so instead can run
myGSEA.df <- as_tibble(LogFC.RecurrenceVPrimary.h.GSEA.res)%>%
  arrange(p.adjust) %>%
  mutate(Phenotype = "Enrichment in Recurrence")

# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.df[1:16,], aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  xlab("")+
  ylab("Gene Set")+
  ggtitle("Hallmark Gene Sets Enriched in BRCA1/2 Recurrent Tumors \nCompared to BRCA1/2 Primary Tumors")+
  theme_bw()
