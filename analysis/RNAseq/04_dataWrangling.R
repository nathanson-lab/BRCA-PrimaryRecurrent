#Objectives ----
# 1 - Filter and normalize the imported abundance data
# 2 - Visualize (violin plots) the impact of filtering and normalization on our data
  # Can also assess data for outlying samples here 
# 3 - Generate filtered, normalized DGELists for entire cohort
  #as well as for breast and ovarian samples separately
#Required inputs
  #Txi_gene object generated in Tximport_fromStringtie.R
  #targets dataframe generated in Tximport_fromStringtie.R

# Load packages -----
library(tidyverse) 
library(edgeR) 
library(matrixStats) 
library(cowplot) 

# Store counts and TPM from Tximport ----
#Note: Abundance data are TPM (from Tximport), 
#while the counts are read counts mapping to each gene or transcript
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts

# capture sample labels from the study design file
sampleLabels <- targets$sample

# Generate summary stats for your data ----
# Calculate summary stats for each gene using matrixStats package, then add these to your data matrix
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))

# Produce a scatter plot of the transformed data ----
scatter<-ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=25, size=2) +
  geom_smooth(method=lm) +
  labs(y="Median", x = "Standard deviation",
       title="Transcripts per million (TPM) by Gene, across Cohort",
       subtitle="Unfiltered, Non-normalized")+
  theme_classic() +
  theme_dark() + 
  theme_bw()
scatter

# Make a DGElist from your counts, and plot ----
myDGEList <- DGEList(myCounts)
#save DGEList object to working directory;can be easily shared and loaded with load(file = "myDGEList)
save(myDGEList, file = "myDGEList")

# use the 'cpm' function from EdgeR to get counts per million
cpm <- cpm(myDGEList) 
colSums(cpm) #each column of counts per million should add up to 1e+06 (1 million)
log2.cpm <- cpm(myDGEList, log=TRUE)

# 'coerce' your data matrix to a dataframe so that you can use tidyverse tools on it
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
log2.cpm.df
# add your sample names to this dataframe (we lost these when we read our data in with tximport)
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
# use the tidy package to 'pivot' your dataframe (from wide to long)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = 2:67, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

# plot this pivoted data
p1<-ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="\n\n Log2 expression", x = "Sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="Unfiltered, Non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))
  
# can also try using coord_flip() at the end of the ggplot code

# Filter your data ----
#first, take a look at how many genes or transcripts have no read counts at all
# 1st - 'myDGEList$counts==0' returns a new 'logical matrix' where each observation (gene) is evaluated (TRUE/FALSE) for each variable (sample) as to whether it has zero counts
# 2nd - passing this logical matrix to 'rowsums' allows you to sum the total number of times an observation was 'TRUE' across all samples
# 3rd - adding the '==66' is a simple way of asking how many genes had 0 counts (TRUE) for all samples in our dataset
table(rowSums(myDGEList$counts==0)==66)

# now set some cut-off to get rid of genes/transcripts with low counts
# again using rowSums to tally up the 'TRUE' results of a simple evaluation
# how many genes had more than 1 CPM (TRUE) in at least 3 samples 
keepers <- rowSums(cpm>1)>=3 #using a general n=3 cutoff here 
# subset DGEList to just include the keepers
myDGEList.filtered <- myDGEList[keepers,]
#check size of filtered DGEList
dim(myDGEList.filtered)

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
# pivot this FILTERED data, just as you did earlier
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = 2:67, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)


p2<-ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="\n\n Log2 expression", x = "Sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="Filtered, Non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))

# Normalize your data ----
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")

# use the 'cpm' function from EdgeR to get counts per million from your normalized data
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
# pivot this NORMALIZED data, just as you did earlier
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = 2:67, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)


p3<-ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="\n\n Log2 expression", x = "Sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="Filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))

# we'll use the 'plot_grid' function from the cowplot package to put these together in a figure
plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

#Make separate DGELists for breast and ovarian samples separately ----

#make a list of column indices pertaining to breast samples
breast<-targets$sample[targets$tissue_origin=="breast"]
breast.cols<-which(targets$sample %in% breast)
#make a DGEList for breast samples alone, based on only the breast columns in Txi_gene$counts
breast.DGEList <- DGEList(Txi_gene$counts[,breast.cols])
save(breast.DGEList, file = "breast.DGEList")
#get cpm from DGEList
breast.cpm <- cpm(breast.DGEList)
#decide on rows (genes) to keep
breast.keepers <- rowSums(breast.cpm>1)>=3 #based on smallest breast group : BRCA2 normal breast, n=3
#filter and normalize DGEList
breast.DGEList.filtered <- breast.DGEList[breast.keepers,]
breast.DGEList.filtered.norm <- calcNormFactors(breast.DGEList.filtered, method = "TMM")

#Repeat for ovary
ovary <- targets$sample[targets$tissue_origin=="ovary"]
ovary.cols<-which(targets$sample %in% ovary)
ovary.DGEList <- DGEList(Txi_gene$counts[,ovary.cols])
save(ovary.DGEList, file = "ovary.DGEList")
ovary.cpm <- cpm(ovary.DGEList)
ovary.keepers <- rowSums(ovary.cpm>1)>=2 #based on smallest ovary group : BRCA2 normal fallopian tube, n=2
ovary.DGEList.filtered <- ovary.DGEList[ovary.keepers,]
ovary.DGEList.filtered.norm <- calcNormFactors(ovary.DGEList.filtered, method = "TMM")
