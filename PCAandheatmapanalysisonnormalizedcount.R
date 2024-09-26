normalizedcount <- read.table("data/normalized_counts.txt", header = T, row.names = 1)
data <- read.table("data/raw_counts.txt", header = T, row.names = 1)
View(normalizedcount)

View(normalizedcount)
class(normalizedcount)
library(readxl)
meta <- read_excel("meta/RNA-seq samples for data analysisnew.xlsx")
PCA <- prcomp(normalizedcount)
summary(PCA)
#we are able to retain about %75 of the total variance in the data
biplot(PCA)
var <- PCA$sdev^2
pve <- var/sum(var)
cumsum(pve)       

#Another trial
pca <- prcomp(normalizedcount, scale = TRUE, center = TRUE, retx = T)
names(pca)
summary(pca)
pca$rotation
dim(pca$x)
pca$x
biplot(pca, main = "Biplot", scale = 0)
par(mfrow = c(1, 2))
plot(pca$rotation[,1:2], col = meta$treatment, pch = 19)
legend("bottomleft", legend = unique(meta$treatment),
       col = 1:2, pch = 19, bty = "n")
       
#for the visualization analyses, matching control
rownames(meta)
rownames(meta) <- meta$'ID CRIBI'
remove_rownames(id) %>% has_rownames()
meta <- rownames_to_column(id, var = "ID CRIBI") %>% as_tibble()
meta[,c("id", "MY ID", "mouse", "TA", "treatment", "...6", "...7")]
rownames(meta) <- c("EUP", "EUQ", "EUR", "EUS", "EUT", "EUU", "EUV", "EUW", "EUX", "EUY", "EUZ", "EVA", "EVD", "EVE", "C43", "C44", "C45", "C46", "C47", "C48", "C49", "C50", "C51", "C52", "C53", "C54", "C55", "C56", "C57", "C58", "C59", "C60")
colnames(normalizedcount)
all(colnames(data) %in% rownames(meta)) ##
all(colnames(data) == rownames(meta)) ##

dds_DS <- DESeqDataSetFromMatrix(countData = round(normalizedcount),
                                 colData = meta,
                                 design = ~ treatment)
normalizedcount <- as.matrix(normalizedcount)
head(normalizedcount)

conditiontreatment <- factor(rep(c("MCS", "ActA+IL6", "MCS", "ActA+IL6", "MCS", "ActA+IL6", "MCS", "ActA+IL6", "MCS", "ActA+IL6", "MCS", "ActA", "MCS", "ActA", "MCS", "ActA", "MCS", "ActA", "MCS", "ActA", "MCS", "IL6", "MCS", "IL6", "MCS", "IL6", "MCS", "IL6", "MCS", "IL6", "MCS", "IL6")))
coldata <- data.frame(row.names=colnames(normalizedcount), conditiontreatment)
head(coldata)
dds <- DESeqDataSetFromMatrix(countData = round(normalizedcount), colData = coldata, design = ~ conditiontreatment)
dds <- DESeq(dds)

# Plot Dispersions:
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

library(pcaMethods)
# Principal Components Analysis on normalized count file, and at the end I also created separate PCA plot with clustering by ellipses
plotPCA(rld, intgroup = c("conditiontreatment")) + stat_ellipse()
install.packages("devtools")
library(devtools)
devtools::install_github("cmartin/ggConvexHull")
library(ggConvexHull)
ggplot(PCA_plot_sartori_ellipses, aes(PC1, PC2, color=sample(), shape=conditiontreatment, group=treatment)) +
  geom_point(size=3) +
  stat_ellipse()
            
vignette("extending-ggplot2", package = "ggplot2")
plotPCA(rld, intgroup = c("ID CRIBI"))

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(conditiontreatment))])


# Sample distance heatmap on treatment group 
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "blue", "pink"),
          ColSideColors=mycols[conditiontreatment], RowSideColors=mycols[conditiontreatment],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Get differential expression results #this can't be done on normalized data, it was a trial of something 
#res <- results(dds)
#table(res$padj<0.05)
#it is observed that 622 differentially expressed genes with adjusted p value <= 0.05

## Order by adjusted p-value
#res <- res[order(res$padj), ]
## Merge with normalized count data
#resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
#names(resdata)[1] <- "Gene"
#head(resdata)
## Write results
#write.csv(resdata, file="diffexpr-results.csv",quote = FALSE,row.names = F)

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

library(pheatmap)
library(DESeq2)
library(RColorBrewer)
library(RSQLite)

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## This is Stephen Turner's code:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), points(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

install.packages("calibrate")
## Plots to Examine Results:

## Volcano plot with "significant" genes labeled on normalized data (this was also another trial)
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()

#createheatmap more advance = by clustering and adding labels
heatmap(MSnbase::exprs(normalized_counts)) #didnotuse 
install.packages("MSnbase")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")
yes
library(AnnotationDbi)

### Annotate our heatmap (optional)
annotation <- data %>% 
  dplyr::select(rownames, colnames) %>% 
  data.frame(row.names = "rownames")

### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap (this belongs the latest heat map that I did (the latest in the slide as well))
pheatmap(normalizedcount, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

install.packages("heatmaply")
library(heatmaply)

heatmaply(normalizedcount)


#compute the variance of each gene across samples
V <- apply(normalizedcount, 1, var)
#sort the results by variance in decreasing order 
#and select the top 100 genes 
selectedGenes <- names(V[order(V, decreasing = T)][1:100])
pheatmap(normalizedcount[selectedGenes,], scale = 'row', show_rownames = FALSE)


pheatmap(normalizedcount[selectedGenes,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = coldata)

##PCA (This is just a trial. I wanted to create a PCA table in a different way)
library(DESeq2)
# extract normalized counts from the DESeqDataSet object
normalizedcount <- DESeq2::counts(dds, normalized = TRUE)

# select top 500 most variable genes
#selectedGenes <- names(sort(apply(rld, 1, var), 
                            #decreasing = TRUE)[1:500])

#plotPCA(rld[selectedGenes,], 
        #colData = as.numeric(colData$treatment), adj = 0.5, 
        #xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.6))

library(pcaMethods)

dds_DS <- DESeqDataSetFromMatrix(countData = round(normalizedcount),
                                 colData = meta,
                                 design = ~ treatment)
normalizedcount <- as.matrix(normalizedcount)
head(normalizedcount)

#extracting only MCS and ActA+IL6 group from meta data is not working. Giving an error. 
conditiontreatment <- factor(rep(c("MCS", "ActA+IL6", "MCS", "ActA+IL6", "MCS", "ActA+IL6", "MCS", "ActA+IL6", "MCS", "ActA+IL6")))
coldata <- data.frame(row.names=colnames(normalizedcount), conditiontreatment)
head(coldata)
dds <- DESeqDataSetFromMatrix(countData = round(normalizedcount), colData = coldata, design = ~ conditiontreatment)
dds <- DESeq(dds)

#I tried to analyse just two group comparison by creating a separate file and eliminate only IL6 and only ActA. And I tried to analyze solo MCS vs ActA+IL6
normalizedcountnewfixed <- read.table("data/normalized_counts.txt", header = T, row.names = 1)
metanewfixed <- read_excel("meta/RNA-seq samples for data analysisnewfixed.xlsx")
View(data)

View(normalizedcountnewfixed)
class(normalizedcountnewfixed)

View(metanewfixed)
class(metanewfixed)

PCA <- prcomp(normalizedcountnewfixed)
summary(PCA)
#we are able to retain about %75 of the total variance in the data
biplot(PCA)
var <- PCA$sdev^2
pve <- var/sum(var)
cumsum(pve)

data.matrix(normalizedcountnewfixed)

dds_DS <- DESeqDataSetFromMatrix(countData = round(normalizedcountnewfixed),
                                 colData = metanewfixed,
                                 design = ~ treatment)
normalizedcount <- as.matrix(normalizedcount)
head(normalizedcount)

#Note: Don't forget to change file if you want to analyze from scratch. 


