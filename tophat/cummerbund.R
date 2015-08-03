## Collect arguments
args <- commandArgs(TRUE)
dir <- args[1]

# load it up
library("cummeRbund")
cuff <- readCufflinks(dir)
pdf(file=paste(dir, "/", "cummeRbund.pdf", sep=""), width=10, height=10)

# generate the dispersion plot and add a title
d<-dispersionPlot(genes(cuff))
d$labels$title <- "Dispersion Plot\n-\nDispersion (deviation from a threshold) vs FPKM"
d$labels$title.theme <- d$theme
plot(d)

pBoxRep<-csBoxplot(genes(cuff),replicates=T)
plot(pBoxRep)

pDendro<-csDendro(genes(cuff),replicates=T)
plot(pDendro)

pBox<-csBoxplot(genes(cuff))
plot(pBox)

myDistHeat<-csDistHeat(genes(cuff))
plot(myDistHeat)

myDistHeat<-csDistHeat(genes(cuff), replicates=T)
plot(myDistHeat)

PCA<-PCAplot(genes(cuff),"PC1","PC2")
plot(PCA)

PCA<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
plot(PCA)
dev.off();