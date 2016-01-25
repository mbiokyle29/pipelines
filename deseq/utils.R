# FUNCTIONS #
convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

mergeIDList <- function(ensg, symbol) {

	if(length(ensg) != length(symbol)) {
		return(ensg)
	} else {
		for(i in 1:length(ensg)) {
			if(is.na(symbol[i])) {
				symbol[i] = ensg[i]
			}
		}
		return(symbol)
	}
}

plotUpDownSigGenes <- function(results, colNums, rld, title) {
	
	# make the lists
	upgenes <- rownames(head(results[ order( results$log2FoldChange ), ], n=20))
	downgenes <- rownames(head(results[ order( -results$log2FoldChange ), ], n=20))

	# this gives us the rows we want
	rows <- match(upgenes, row.names(rld))

	mat <- assay(rld)[rows,colNums]
	conv <- convertIDs(row.names(mat), "ENSEMBL", "SYMBOL", org.Hs.eg.db, ifMultiple="useFirst")
	row.names(mat) <- mergeIDList(row.names(mat), conv)
	mat <- mat - rowMeans(mat)

	df <- as.data.frame(colData(rld)[c("diffVar","ebvVar")])
	pheatmap(mat, fontsize=5, annotation_col=df, main=paste(title,"top 20 up genes"))


	# this gives us the rows we want
	rows <- match(downgenes, row.names(rld))

	mat <- assay(rld)[rows,colNums]
	conv <- convertIDs(row.names(mat), "ENSEMBL", "SYMBOL", org.Hs.eg.db, ifMultiple="useFirst")
	row.names(mat) <- mergeIDList(row.names(mat), conv)
	mat <- mat - rowMeans(mat)

	df <- as.data.frame(colData(rld)[c("diffVar","ebvVar")])
	pheatmap(mat, fontsize=5, annotation_col=df, main=paste(title,"top 20 down genes"))
}

cleanUpNames <- function(name) {
	strsplit(strsplit(x, "RNAseq.")[[1]][-1], ".genes.results")[[1]][1]
}