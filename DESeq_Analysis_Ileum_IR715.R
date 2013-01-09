setwd("C:/Dropbox/Macaque_Project/miRNA_analysis/Differential_expression/BLAST_against_Macaca_only")
myPath = file.path("C:/Dropbox/Macaque_Project/miRNA_analysis/Differential_expression/BLAST_against_Macaca_only/miRNA_Expression_Macaca_only_db.tsv")
counts = read.table(myPath, header=TRUE, row.names=1)
metadata = data.frame(
	row.names = colnames(counts),
	condition = c("NA", "NA", "NA", "NA", "NA", "treated10wk", "NA", "untreated", "NA", "NA"),
	libType = c("single-end", "single-end", "single-end", "single-end", "single-end", "single-end", "single-end", "single-end", "single-end", "single-end") )
singleSamples = metadata$libType == "single-end"
countTable = counts[,singleSamples]
condition = metadata$condition[singleSamples]
library("DESeq")
cds = newCountDataSet(countTable, condition)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
#head( counts( cds, normalized=TRUE) ) # Prints library-normalized read counts
cds = estimateDispersions(cds, method="blind", sharingMode="fit-only")
#plotDispEsts(cds) # Plot dispersion estimates

# Examine differential expression between 10wk infected vs non-infected ileum
res = nbinomTest(cds, "untreated", "treated10wk")
#head(res)
#plotMA(res)
#hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="") # Prints histogram of p-values
#hist(res$padj, breaks=100, col="skyblue", border="slateblue", main="") # Prints histogram of adjusted p-values
#resSig2wk = res[res$padj < 0.1,] # With only one biological sample, all gene-by-gene comparisons will have padj=1. Must use pval instead
resSig2wk = res[res$pval < 0.1,]
resSig2wk = na.omit(resSig2wk)
head(resSig2wk[order(resSig2wk$pval),]) # Most significantly differentially expressed genes
head(resSig2wk[order(resSig2wk$foldChange, -resSig2wk$baseMean ), ] ) # Most strongly down-regulated significant genes
head(resSig2wk[order(-resSig2wk$foldChange, -resSig2wk$baseMean ), ] ) # Most strongly up-regulated significant genes
write.table(resSig2wk[order(resSig2wk$pval), ], file="IR715_10WkSigDifExp.tsv", quote=FALSE, sep="\t", row.names = FALSE )
write.table(resSig2wk[order(resSig2wk$foldChange, -resSig2wk$baseMean ), ], file="IR715_10WkSigMostDownReg.tsv", quote=FALSE, sep="\t", row.names = FALSE )
write.table(resSig2wk[order(-resSig2wk$foldChange, -resSig2wk$baseMean ), ], file="IR715_10WkSigMostUpReg.tsv", quote=FALSE, sep="\t", row.names = FALSE )
