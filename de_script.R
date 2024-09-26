## Gene-level differential expression analysis using DESeq2

install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggnewscale")

install.packages("DEGreport")

library(tidyverse)
library(RColorBrewer)
install.packages("RSQLite")
library(DESeq2)
library(pheatmap)
library(DEGreport)
install.packages("readxl")
library(readxl)

data <- read.table("data/raw_counts_xxx.txt", header = T, row.names = 1)
meta <- read_excel("meta/RNA-seq samples for data analysisnew.xlsx")


class(data)
class(meta)

View(data)
View(meta)

## firstly normalization of raw counts data. To this, we need to determine the appropriate statistical model. We will need information about the distribution of counts.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggplot")

library(ggplot2)

ggplot(data) +
  geom_histogram(aes(x = EUP), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

#majority of genes with counts of zero

ggplot(data) +
  geom_histogram(aes(x = EUP), stat = "bin", bins = 200) +
  xlim(-5, 500) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

#these are for one gene, it would be better for me to see general 
#how our counts are distributed but I could not find the function for more than one genes. 

#this plot shows us the low number of counts are related to large number of genes. 
#And the line that goes to the right unlimitless indicates that there is no upper limit for expression. Note that RNA-seq data does not give normal distribution unlike microarray data

#for modelling our count data, deciding our models. 
#If it is count data (as our data) then it should be negative binomial model. This means that mean < variance
mean_counts <- apply(data[, 3:5], 1, mean)
variance_counts <- apply(data[, 3:5], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x = mean_counts, y = variance_counts)) +
  geom_line(aes(x = mean_counts, y = mean_counts, color = "red")) +
  scale_y_log10() +
  scale_x_log10()
  
vignette("DESeq2")

head(meta)

head(data)

## for DE analysis, 1- Count normalization 

#create DESeq2Dataset object
library(DESeq2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")

install.packages("SummarizedExperiment")
install.packages("pasilla")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pasilla")

library("pasilla")

#bringing meta and raw count data for DESEq2 and make readable for DESEq2
#firstly, I have tried to do with tximport. But then I noticed that with FASTQ raw count data, I cannot use tximport method. 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")

install.packages("readr")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximportData")

library("tximport")
library("readr")
library("tximportData")
dir <- system.file("data", package="tximportData")
meta <- read_excel(file.path(dir,"/Users/asudeberber/Desktop/DEAnalysis/pathxxxx/meta/RNA-seq samples for data analysisnew.xlsx"))
getwd()

meta$treatment <- factor(rep(c("MCS", "ActA+IL6", "MCS", "ActA+IL6", "MCS", "ActA+IL6", "MCS", "ActA+IL6", "MCS", "ActA+IL6", "MCS", "ActA", "MCS", "ActA", "MCS", "ActA", "MCS", "ActA", "MCS", "ActA", "MCS", "IL6", "MCS", "IL6", "MCS", "IL6", "MCS", "IL6", "MCS", "IL6", "MCS", "IL6")))

library(tibble)

#then for the count normalization, I applied count matrix input way. 
#By this I mean that to ensure the row names of the metadata dataframe are present and in the same order as the column names of the counts dataframe.
#(matching meta data and our data was quite challenge for me) 
#Then, Create a DESeqDataSet object. Then, Generate the normalized counts. 

meta %>% as_tibble(meta = "id")
colnames(meta)[1] <- "id"
has_rownames(id)
has_rownames(mouse)



remove_rownames(id) %>% has_rownames()
meta <- rownames_to_column(id, var = "ID CRIBI") %>% as_tibble()

rownames(meta) <- meta$'ID CRIBI'


meta[,c("id", "MY ID", "mouse", "TA", "treatment", "...6", "...7")]

#tximport (Transcript abundance files)
##firstly, I have tried to do with tximport. 
#But then I noticed that with FASTQ raw count data, I cannot use tximport method. 
files <- file.path(dir,"salmon", meta$id, "/Users/asudeberber/Desktop/DEAnalysispathpath/data/raw_counts.....txt")
names(files) <- meta$id
tx2gene <- read.table(file.path(dir, "/Users/asudeberber/Desktop/DEAnalysispathxxxx/data/raw_counts...txt"), header=TRUE, row.names = 1)
txi <- tximport(files, type="salmon", tx2gene=tx2gene)



x <- rownames(meta)
y <- colnames(data)
rownames(meta) <- c("EUP", "EUQ", "EUR", "EUS", "EUT", "EUU", "EUV", "EUW", "EUX", "EUY", "EUZ", "EVA", "EVD", "EVE", "C43", "C44", "C45", "C46", "C47", "C48", "C49", "C50", "C51", "C52", "C53", "C54", "C55", "C56", "C57", "C58", "C59", "C60")
y <- c("EUP", "EUQ", "EUR", "EUS", "EUT", "EUU", "EUV", "EUW", "EUX", "EUY", "EUZ", "EVA", "EVD", "EVE", "C43", "C44", "C45", "C46", "C47", "C48", "C49", "C50", "C51", "C52", "C53", "C54", "C55", "C56", "C57", "C58", "C59", "C60")
all(x %in% y)
match(rownames(meta), colnames(data))

rownames(meta)
colnames(data)
genomic_idx <- match(rownames(meta), colnames(data))
data_ordered <- data[,genomic_idx]
head(data_ordered)

ncol(data)
nrow(meta)

length(data)
length(meta)

all(colnames(data) %in% rownames(meta)) ##
all(colnames(data) == rownames(meta)) ##

## same variables but in different order. 
#Matching problem is solved.


has.row.names = function(df) {
  !all(row.names(df)==seq(1, nrow(df)))
}

vector1 <- rownames(meta)
vector2 <- colnames(data)

vector1 <- meta$`ID CRIBI`
match(vector1, colnames(data))

match(vector1, colnames(data))
all(colnames(data) == rownames(meta))

meta$rownames <- NULL
View(meta)

install.packages("write.xlsx")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("write.xlsx")
library(WriteXLS)
library(WriteXLS)

WriteXLS(meta, files/meta, col.names=TRUE, row.names=TRUE, append = FALSE)


match(colnames(data), rownames(meta))


all(colnames(data) == rownames(reordered_metadata)) 

for(i in ((length(rownamesmeta)+1):length(colnamesdata))) 
  +{rownamesmeta = c(colnamesdata, 0)}
#I have tried to some functions to reorder them but it did not work with functions.
#At the end I got them in order by fixing meta data in excel file. 


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")

library(SummarizedExperiment)

installed.packages("magrittr")
library(magrittr)


# check the order
rownames(meta) <- meta$`ID CRIBI`

rownames(meta)
colnames(data)

match(rownames(meta), colnames(data))
all(colnames(data) %in% rownames(meta)) ##
all(colnames(data) == rownames(meta)) ##

## our rownames(meta) and colnames(data) are same. 
#They should have matched easily and all function should be true.
## hard to understand why it gives FALSE result even if they are matching and same 

idx <- match(rownames(meta), colnames(data))
reordered_meta <- meta[idx,]
View(reordered_meta)
rownames(reordered_meta) #it was not necessary anymore to use reordered_meta file 

#Even if they are in same order, with 'all() %in%' and 'all() == ' functions, 
#I got FALSE result but I had to continue to create DESeq2 object in order to complete DEGs analysis. Also, there was no problem with my meta and data file, rownames and colnames were matching. 
#Matching problem is solved.

### Check that sample names match in both files

library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)

countData <- data
colData <- meta
#it is needed to create DESeqmatrix
ds <- DESeqDataSetFromMatrix(countData = round(data),
                             colData = meta,
                             design = ~ treatment)

View(counts(ds))
##normalization
dds <- estimateSizeFactors(ds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
##to save normalized data matrix to file 
write.table(normalized_counts, file="data/normalized_countsforDESeq.txt", sep="\t", quote=F, col.names=NA)
#this normalizedcount file will be useful for downstream visualization of results. 
#But we have already normalized count file which was provided to me.
#So for PCA and Heatmap analysis, I used that normalized count file. 

#However, After read count, for the quality control of count data: Sample level QC = also determined outliers and it will be used PCA and hierarchical clustering methods (unsupervised)
#to moderate the variance across, to improve PCA and clustering visualization methods. 
#Apply rlog transformation to moderate variance across the mean
rld <- rlog(dds, blind=TRUE)
rld <- vst(dds, blind=TRUE) #vst() is faster way

## PCA for DESeq2 = for Quality control/Quality Assessment 
plotPCA(rld, intgroup="treatment")
#for only ActA+IL6 treatment there is an outlier and some dispersedness

##hierarchical clustering = Quality Control/Quality Assessment for DESeq2
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat) #compute pairwise correlation values
head(rld_cor)
pheatmap(rld_cor)

heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)
?pheatmap
#I observed high correlations across the board (generally between >0.99-0.97 range)
#These plots are indicated that we are using the good quality of data and 
#it can be proceeded with differential expression analysis.

#For differential expression analysis 
design <- ~ treatment
dds <- DESeqDataSetFromMatrix(countData = round(data), colData = meta, design = ~ treatment)
dds <- DESeq(dds)
#DESeq() function completed all the steps for differential expression analysis
#1-estimating size factors
#2-estimating dispersions
#3-gene-wise dispersion estimates
#4-mean-dispersion relationship
#5-final dispersion estimates
#6-fitting model and testing

# 1- estimate size factor #DESeq2 automatically calculate sizefactors 
sizeFactors(dds)
#total numbers of raw counts per sample. Reconciled, how do numbers correlate with sizefactors
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T)) 
#How do the values across samples compare with the total counts taken for each sample?

# 2- estimate gene wise dispersion
#DESeq dispersion = count dispersion = black dot plot

# 3- Plot gene-wise dispersion estimates (shrinkage)
plotDispEsts(dds)

#log2 foldchanges estimate = gives calculation or measure between conditions. 
#LFC (shrunken log2 foldchanges) gives information about dispersion of gene. For example two genes can have similar normalized count but differ LFC degree shrinkage.


#hypothesis testing and model fitting = I used Wald test. Because I want to compare two groups. Wald test bases upon p-value.
## I defined contrasts, and extracted results table, and shrinked the log2 fold changes.

contrast_ctrlvstrtdboth <- c("treatment", "MCS", "ActA+IL6")

res_tablectrlvstrtdboth_unshrunken <- results(dds, contrast=contrast_ctrlvstrtdboth, alpha = 0.05)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")
library(apeglm)
res_tablectrlvstrtdboth <- lfcShrink(dds, contrast = contrast_ctrlvstrtdboth, res = res_tablectrlvstrtdboth_unshrunken)
res_tablectrlvstrtdboth <- lfcShrink(dds, coef = 3, type = "apeglm")
plotMA(res_tablectrlvstrtdboth_unshrunken, ylim=c(-2,2))
#The genes that are significantly DE are colored to be easily identified. 
#This is also a great way to illustrate the effect of LFC shrinkage

#for shrunken result
plotMA(res_tablectrlvstrtdboth, ylim=c(-2,2))
#to evaluate the magnitude of fold changes and how they are distributed relative to mean expression. Generally, we would expect to see significant genes across the full range of expression levels.

#to explore results
class(res_tablectrlvstrtdboth)
mcols(res_tablectrlvstrtdboth, use.names=T)

res_tablectrlvstrtdboth %>% data.frame() %>% View()
write.table(res_tablectrlvstrtdboth, file="results/res_tablectrlvstrtdboth.txt", sep="\t", quote=F, col.names=NA)
write.table(res_tablectrlvstrtdboth, file="results/res_tablectrlvstrtdboth.xlsx")
#FDR < 0.05 mean  the proportion of false positives we expect amongst our differentially expressed genes is 5%.
#NOTE: on p-values set to NA
#If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p-value and adjusted p-value will all be set to NA.
#If a row contains a sample with an extreme count outlier then the p-value and adjusted p-value will be set to NA. These outlier counts are detected by Cook’s distance.
#If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p-value will be set to NA.

## Define contrasts, extract results table and shrink log2 fold changes
contrast_ctrlvstrtdActA <-  c("treatment", "MCS", "ActA")

res_tablectrlvstrtdActA_unshrunken <- results(dds, contrast=contrast_ctrlvstrtdActA, alpha = 0.05)

res_tablectrlvstrtdActA <- lfcShrink(dds, coef = 3, type = "apeglm")

plotMA(res_tablectrlvstrtdActA_unshrunken, ylim=c(-2,2))
plotMA(res_tablectrlvstrtdActA, ylim=c(-2,2))

class(res_tablectrlvstrtdActA)
mcols(res_tablectrlvstrtdActA, use.names=T)
res_tablectrlvstrtdActA %>% data.frame() %>% View() #same!!
write.table(res_tablectrlvstrtdboth, file="results/res_tablectrlvstrtdboth.txt", sep="\t", quote=F, col.names=NA)

## Define contrasts, extract results table and shrink log2 fold changes
contrast_ctrlvstrtdIL6 <-  c("treatment", "MCS", "IL6")

res_tablectrlvstrtdIL6_unshrunken <- results(dds, contrast=contrast_ctrlvstrtdIL6, alpha = 0.05)

res_tablectrlvstrtdIL6 <- lfcShrink(dds, coef = 3, type = "apeglm")

plotMA(res_tablectrlvstrtdIL6_unshrunken, ylim=c(-2,2))
plotMA(res_tablectrlvstrtdIL6, ylim=c(-2,2))

summary(res_tablectrlvstrtdboth)
# the function also reports the number of genes that were tested (genes with non-zero total read count), 
#and the number of genes not included in multiple test correction due to a low mean count.

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

res_tablectrlvstrtdboth_tb <- res_tablectrlvstrtdboth %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#By the way, I defined the first column name as 'gene'

sigctrlvstrtdboth <- res_tablectrlvstrtdboth_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
sigctrlvstrtdboth
#total 789 genes are significantly differentially expressed in MCS vs ActA+IL6 group

res_tablectrlvstrtdActA_tb <- res_tablectrlvstrtdActA %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigctrlvstrtdActA <- res_tablectrlvstrtdActA_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
sigctrlvstrtdActA
#like above again total 789 genes are significantly differentially expressed in MCS vs ActA group. Also same genes are sig DE again here.

res_tablectrlvstrtdIL6_tb <- res_tablectrlvstrtdIL6 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigctrlvstrtdIL6 <- res_tablectrlvstrtdIL6_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
sigctrlvstrtdIL6
##like above again total 789 genes are significantly differentially expressed in MCS vs IL6 group. Also same genes are sig DE again here.


#A tibble: 789 × 6
#gene               baseMean log2FoldChange lfcSE   pvalue     padj
#<chr>                 <dbl>          <dbl> <dbl>    <dbl>    <dbl>
#1 ENSMUSG00000000031  2812.            0.810 0.229 2.27e- 5 7.63e- 4
#2 ENSMUSG00000000318    73.7           3.44  0.326 4.29e-27 4.79e-24
#3 ENSMUSG00000000916    49.6          -0.886 0.487 2.82e- 3 3.29e- 2
#4 ENSMUSG00000000958     6.57          1.76  0.712 4.89e- 4 9.05e- 3
#5 ENSMUSG00000001128    29.0           1.64  0.335 6.07e- 8 4.57e- 6
#6 ENSMUSG00000001247     9.94          1.03  0.488 1.55e- 3 2.15e- 2
#7 ENSMUSG00000001627   553.           -0.846 0.162 1.30e- 8 1.20e- 6
#8 ENSMUSG00000001750    48.2           0.590 0.224 7.59e- 4 1.27e- 2
#9 ENSMUSG00000001865     9.53          1.58  0.458 2.45e- 5 8.09e- 4
#10 ENSMUSG00000001930    77.4           0.632 0.208 2.07e- 4 4.63e- 3
# … with 779 more rows

#same result for all three conditions. for all three conditions = MCS vs ActA+IL6, MCS vs ActA, MCS vs IL6


#alternative approach to see in separate table file
results(dds, contrast = contrast_ctrlvstrtdboth, alpha = 0.05, lfcThreshold = 0.58) %>% data.frame() %>% View()
results(dds, contrast = contrast_ctrlvstrtdboth, alpha = 0.05, lfcThreshold = 0.58)
results(dds, contrast = contrast_ctrlvstrtdActA, alpha = 0.05, lfcThreshold = 0.58) %>% data.frame() %>% View()
results(dds, contrast = contrast_ctrlvstrtdActA, alpha = 0.05, lfcThreshold = 0.58)
results(dds, contrast = contrast_ctrlvstrtdIL6, alpha = 0.05, lfcThreshold = 0.58) %>% data.frame() %>% View()
results(dds, contrast = contrast_ctrlvstrtdIL6, alpha = 0.05, lfcThreshold = 0.58)

#visualizing results
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)




normalized_counts <- read.table("data/normalized_countsforDESeq.txt", header = T, row.names = 1)
#I used normalized count table which I created in the beginning of DESeq analysis process above). 
#But here, I might have used the (ready) normalized count file
View(normalized_counts)
# Create tibbles including row names
ourRNAseq_meta <- meta %>% 
  rownames_to_column(var="MCS") %>% 
  as_tibble()

ourRNAseq_meta <- meta %>% 
  rownames_to_column(var="treatment") %>% 
  as_tibble()
normalized_counts <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
#it should be for single gene 

# Plot expression for single gene for example
plotCounts(dds, gene="ENSMUSG00000000001", intgroup="treatment")

# Save plotcounts to a data frame object
d <- plotCounts(dds, gene="ENSMUSG00000000001", intgroup="treatment", returnData=TRUE)

# Plotting the gene normalized counts, using the samplenames (rownames of d as labels)
ggplot(d, aes(x = treatment, y = count, color = treatment )) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("ENSMUSG00000000001") +
  theme(plot.title = element_text(hjust = 0.5))
#these functions are creating are useful plots for the single gene. 

#I want to see the table and volcano plot table for top 20 DEGs, use ggplot2 to plot multiple genes
## Order results by padj values
top20_sigctrlvstrtdboth_genes <- res_tablectrlvstrtdboth_tb %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=20) 		#Extract the first 20 genes
View(top20_sigctrlvstrtdboth_genes)

#normalized counts for top 20 sig 
top20_sigctrlvstrtdboth_norm <- normalized_counts %>%
  filter(gene %in% top20_sigctrlvstrtdboth_genes)
top20_sigctrlvstrtdboth_norm %>% data.frame %>% View()
write.table(top20_sigctrlvstrtdboth_norm, file="results/top20genesigDEctrlvsActA+IL6.txt", sep="\t", quote=F, col.names=NA)

# Gathering the columns to have normalized counts to a single column
gathered_top20_sigctrlvstrtdboth <- top20_sigctrlvstrtdboth_norm %>%
  gather(colnames(top20_sigctrlvstrtdboth_norm)[3:34], key = "ID CRIBI", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(gathered_top20_sigctrlvstrtdboth)
write.table(gathered_top20_sigctrlvstrtdboth, file="results/gathered_top20genesigDEctrlvsActA+IL6.txt", sep="\t", quote=F, col.names=NA)

gathered_top20_sigctrlvstrtdboth <- inner_join(ourRNAseq_meta, gathered_top20_sigctrlvstrtdboth)

#Actually I used tibble version of DESeq result and I used the MCS vs. ActA+IL6 one. 
#Because results of all groups from DESeq gave the same result table.
ggplot(gathered_top20_sigctrlvstrtdboth) +
  geom_point(aes(x = gene, y = normalized_counts, color = meta$treatment)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))


#volcanoplot
res_tablectrlvstrtdboth_tb <- res_tablectrlvstrtdboth_tb %>% 
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58) #ctrlvstrtdboth means all genes from our data


ggplot(res_tablectrlvstrtdboth_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("Volcano Plot table (All genes)") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

#to see the top 10 genes (lowest padj)
## Create a column to indicate which genes to label
res_tablectrlvstrtdboth_tb <- res_tablectrlvstrtdboth_tb %>% arrange(padj) %>% mutate(genelabels = "")

res_tablectrlvstrtdboth_tb$genelabels[1:10] <- res_tablectrlvstrtdboth_tb$gene[1:10]

View(res_tablectrlvstrtdboth_tb)
write.table(res_tablectrlvstrtdboth_tb, file="results/VolcanoplottableDEGbothmeansallgenes.txt", sep="\t", quote=F, col.names=NA)

ggplot(res_tablectrlvstrtdboth_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_OE)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("Volcano Plot (All genes)") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

#GSEA functional analyses GO and KEGG
library(org.Mm.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)

#set threshold
padj <- 0.05

allourgenes <- as.character(res_tablectrlvstrtdboth_tb$gene)
siggenesofallourgenes <- stats::filter(res_tablectrlvstrtdboth_tb, padj < 0.05) #extract significant ones
siggenesofallourgenes <- dplyr::filter(res_tablectrlvstrtdboth_tb, padj < 0.05)
siggenesofallourgenes <- as.character(siggenesofallourgenes$gene)
siggenesofallourgenes

#I used all genes from differentially expressed and significantly differentially expressed results table 
## Run GO enrichment analysis 
ego <- enrichGO(gene = siggenesofallourgenes, 
                universe = allourgenes,
                keyType = "ENSEMBL",
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.csv(cluster_summary, "results/clusterProfiler_fromallourtosiggenes.csv")

#The dotplot shows the number of genes associated with the first 50 terms (size) 
#and the p-adjusted values for these terms (color). 
#This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
## Dotplot 
#Top 50 and top 20 have been tried respectively
dotplot(ego, showCategory=20)

## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
ego2 <- pairwise_termsim(ego, method="JC", semData = NULL, showCategory = 50)
?pairwise_termsim
library(GOSemSim)
library(enrichplot)
library(clusterProfiler)
library(cluster)
emapplot(ego2)
library(ggnewscale)
emapplot_cluster(ego)
library(dplyr)

## Remove any NA values
sartori_entrez <- dplyr::filter(res_tablectrlvstrtdboth_tb, entrez != "NA")
rlang::last_error()
dplyr::filter()

## Remove any NA values
siggenesofallourgenes <- stats::filter(siggenesofallourgenes, gene != "NA") ##none

## Remove any Entrez duplicates
entrezIDourgenetoKEGG <- entrezIDourgenetoKEGG[which(duplicated(sartori_entrez$gene) == F), ]

## Extract the foldchanges
foldchanges <- res_tablectrlvstrtdboth_tb$log2FoldChange
## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- entrezIDourgenetoKEGG$kegg
## Sort fold changes in decreasing order

foldchanges <- sort(foldchanges, decreasing = TRUE)
gene <- sort(gene, decreasing = TRUE)

gene <- c("80911", "56348", "72535", "14187", "13382", "13807")
entrezIDourgenetoKEGG <- bitr_kegg(gene, fromType = "kegg", toType = "Module", organism = "mmu")

head(foldchanges)
## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(gene = entrezIDourgenetoKEGG, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "mmu", # supported organisms listed below
                    nPerm = 1000, # default number permutations
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result

#KEGG pathway trials

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGGREST")

library(KEGGREST)
library(clusterProfiler)
search_kegg_organism('mmu', by='kegg_code')
mouse <- search_kegg_organism('Mus musculus', by='scientific_name')
dim(mouse)
head(mouse)

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
head(kk)

bitr_kegg(gene, fromType = "kegg", toType = "Path", organism = "mmu")
bitr_kegg(gene, fromType = "kegg", toType = "Module", organism = "mmu")

k <- keys(org.Mm.eg.db, keytype = "ENTREZID")
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)
library(xlsx)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(openxlsx)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("patchwork")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggalluvial")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Seurat")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("openxlsx")

gene.list <- select(org.Mm.eg.db, keys = k, columns = c("SYMBOL", "ENSEMBL"), keytype = "ENTREZID")
head(gene.list)
entrez <- filter(gene.list, ENTREZID != "NA")
entrez <- entrez[which(duplicated(entrez$ENTREZID) == F), ]
## Extract the foldchanges
foldchanges <- entrez$log2FoldChange

## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- entrez$ENTREZID
## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)

head(foldchanges)

## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = entrez, 
                    organism = "mmu", # supported organisms listed below
                    pvalueCutoff = 0.05) # padj cutoff value
                    

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result

#another method
gsets_list <- get_gene_sets_list(source = "KEGG",
                                 org_code = "mmu")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathfindR")

install.packages("sos")
require("sos")
findFn("get_gene_sets_list")
library(pathfindR)
gsets_list <- get_gene_sets_list(source = "KEGG",
                                 org_code = "mmu")

mmu_kegg_genes <- gsets_list$gene_sets
mmu_kegg_descriptions <- gsets_list$descriptions

## Save both as RDS files for later use
saveRDS(mmu_kegg_genes, "mmu_kegg_genes.RDS")
saveRDS(mmu_kegg_descriptions, "mmu_kegg_descriptions.RDS")

mmu_kegg_genes <- readRDS("mmu_kegg_genes.RDS")
mmu_kegg_descriptions <- readRDS("mmu_kegg_descriptions.RDS")

class(res_tablectrlvstrtdboth)
diffexpresulttable <- read.csv(file = 'results/diffexpr-results.csv')
View(diffexpresulttable)

knitr::kable(head(res_tablectrlvstrtdboth_tb))
KEGG_output <- run_pathfindR(input = res_tablectrlvstrtdboth_tb,
                                convert2alias = FALSE,
                                gene_sets = "Custom",
                                custom_genes = mmu_kegg_genes,
                                custom_descriptions = mmu_kegg_descriptions)

## for the KEGG pathway, I found the genes but I am getting an error to create KEGG pathway plot
