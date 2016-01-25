#!/usr/bin/env Rscript
library(optparse)
library(DESeq2)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)

source("utils.R")

option_list = list(
	make_option(c("-c", "--counts-file"), type="character"),
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


## Read in raw counts
countMtx <- read.table(opt$c)
countMtx <- round(countMtx)
eNames <- apply(as.matrix(colnames(countMtx)), 1,)

## Run DEseq (filter low count genes, require at least 2 reads)
ddsMF <- countMtx.deseq
ddsMF <- ddsMF[ rowSums(counts(ddsMF)) > 1, ]
ddsMF <- DESeq(ddsMF)


pdf("")

# ### PCA ###
rld <- rlog(ddsMF, blind=FALSE)
plotPCA(rld, intgroup = c(""))
plotPCA(rld, intgroup = c(""))

### Heatmaps top 100 genes based only on variance
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[c()])
pheatmap(mat, annotation_col=df)

dev.off()
