setwd("C:/Dropbox/Macaque_Project/miRNA_analysis/Differential_expression/BLAST_against_Macaca_only")
myPath = file.path("C:/Dropbox/Macaque_Project/miRNA_analysis/Differential_expression/BLAST_against_Macaca_only/miRNA_Expression_Macaca_only_db.tsv")
counts = read.table(myPath, header=TRUE, row.names=1)
metadata = data.frame(
	row.names = colnames(counts),
	condition = c("treated2wk", "treated27wk", "treated2wk", "treated27wk", "treated27wk", "NA", "NA", "NA", "NA", "untreated"),
	libType = c("single-end", "single-end", "single-end", "single-end", "single-end", "single-end", "single-end", "single-end", "single-end", "single-end") )
singleSamples = metadata$libType == "single-end"
countTable = counts[,singleSamples]
condition = metadata$condition[singleSamples]
library("DESeq")
cds = newCountDataSet(countTable, condition)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
#head( counts( cds, normalized=TRUE) ) # Prints library-normalized read counts
cds = estimateDispersions(cds)
#plotDispEsts(cds) # Plot dispersion estimates

# Examine differential expression between 2wk infected vs non-infected jejunum 
res2wk = nbinomTest(cds, "untreated", "treated2wk")
#head(res2wk)
#plotMA(res2wk)
#hist(res2wk$pval, breaks=100, col="skyblue", border="slateblue", main="") # Prints histogram of p-values
#hist(res2wk$padj, breaks=100, col="skyblue", border="slateblue", main="") # Prints histogram of adjusted p-values
resSig2wk = res2wk[res2wk$padj < 0.1,]
resSig2wk = na.omit(resSig2wk)
resSig2wkUp = resSig2wk[resSig2wk$foldChange > 1,]
resSig2wkDown = resSig2wk[resSig2wk$foldChange < 1,]
head(resSig2wk[order(resSig2wk$pval),]) # Most significantly differentially expressed genes
head(resSig2wkUp[order(-resSig2wkUp$foldChange, -resSig2wkUp$baseMean ), ] ) # Most strongly up-regulated significant genes
head(resSig2wkDown[order(resSig2wkDown$foldChange, -resSig2wkDown$baseMean ), ] ) # Most strongly down-regulated significant genes
write.table(resSig2wk[order(resSig2wk$pval), ], file="2WkSigDifExp.tsv", quote=FALSE, sep="\t", row.names = FALSE )
write.table(resSig2wkUp[order(-resSig2wkUp$foldChange, -resSig2wkUp$baseMean ), ], file="2WkSigMostUpReg.tsv", quote=FALSE, sep="\t", row.names = FALSE )
write.table(resSig2wkDown[order(resSig2wkDown$foldChange, -resSig2wkDown$baseMean ), ], file="2WkSigMostDownReg.tsv", quote=FALSE, sep="\t", row.names = FALSE )

# Examine differential expression between 2wk infected vs non-infected jejunum without applying a padj cutoff
resSig2wk = res2wk[res2wk$padj <= 1,]
resSig2wk = na.omit(resSig2wk)
resSig2wkUp = resSig2wk[resSig2wk$foldChange > 1,]
resSig2wkDown = resSig2wk[resSig2wk$foldChange < 1,]
head(resSig2wk[order(resSig2wk$pval),]) # Most significantly differentially expressed genes
head(resSig2wkUp[order(-resSig2wkUp$foldChange, -resSig2wkUp$baseMean ), ] ) # Most strongly up-regulated significant genes
head(resSig2wkDown[order(resSig2wkDown$foldChange, -resSig2wkDown$baseMean ), ] ) # Most strongly down-regulated significant genes
write.table(resSig2wk[order(resSig2wk$pval), ], file="2WkSigDifExp_NOCUTOFF.tsv", quote=FALSE, sep="\t", row.names = FALSE )
write.table(resSig2wkUp[order(-resSig2wkUp$foldChange, -resSig2wkUp$baseMean ), ], file="2WkSigMostUpReg_NOCUTOFF.tsv", quote=FALSE, sep="\t", row.names = FALSE )
write.table(resSig2wkDown[order(resSig2wkDown$foldChange, -resSig2wkDown$baseMean ), ], file="2WkSigMostDownReg_NOCUTOFF.tsv", quote=FALSE, sep="\t", row.names = FALSE )

# Examine differential expression between 27wk infected vs non-infected jejunum 
res27wk = nbinomTest(cds, "untreated", "treated27wk")
#head(res27wk)
#plotMA(res27wk)
#hist(res27wk$pval, breaks=100, col="skyblue", border="slateblue", main="") # Prints histogram of p-values
resSig27wk = res27wk[res27wk$padj < 0.1,]
resSig27wk = na.omit(resSig27wk)
resSig27wkUp = resSig27wk[resSig27wk$foldChange > 1,]
resSig27wkDown = resSig27wk[resSig27wk$foldChange < 1,]
head(resSig27wk[order(resSig27wk$pval),]) # Most significantly differentially expressed genes
head(resSig27wkUp[order(-resSig27wkUp$foldChange, -resSig27wkUp$baseMean ), ] ) # Most strongly up-regulated significant genes
head(resSig27wkDown[order(resSig27wkDown$foldChange, -resSig27wkDown$baseMean ), ] ) # Most strongly down-regulated significant genes
write.table(resSig27wk[order(resSig27wk$pval), ], file="27WkSigDifExp.tsv", quote=FALSE, sep="\t", row.names = FALSE )
write.table(resSig27wkUp[order(-resSig27wkUp$foldChange, -resSig27wkUp$baseMean ), ], file="27WkSigMostUpReg.tsv", quote=FALSE, sep="\t", row.names = FALSE )
write.table(resSig27wkDown[order(resSig27wkDown$foldChange, -resSig27wkDown$baseMean ), ], file="27WkSigMostDownReg.tsv", quote=FALSE, sep="\t", row.names = FALSE )

# Examine differential expression between 27wk infected vs non-infected jejunum without applying a padj cutoff
resSig27wk = res27wk[res27wk$padj <= 1,]
resSig27wk = na.omit(resSig27wk)
resSig27wkUp = resSig27wk[resSig27wk$foldChange > 1,]
resSig27wkDown = resSig27wk[resSig27wk$foldChange < 1,]
head(resSig27wk[order(resSig27wk$pval),]) # Most significantly differentially expressed genes
head(resSig27wkUp[order(-resSig27wkUp$foldChange, -resSig27wkUp$baseMean ), ] ) # Most strongly up-regulated significant genes
head(resSig27wkDown[order(resSig27wkDown$foldChange, -resSig27wkDown$baseMean ), ] ) # Most strongly down-regulated significant genes
write.table(resSig27wk[order(resSig27wk$pval), ], file="27WkSigDifExp_NOCUTOFF.tsv", quote=FALSE, sep="\t", row.names = FALSE )
write.table(resSig27wkUp[order(-resSig27wkUp$foldChange, -resSig27wkUp$baseMean ), ], file="27WkSigMostUpReg_NOCUTOFF.tsv", quote=FALSE, sep="\t", row.names = FALSE )
write.table(resSig27wkDown[order(resSig27wkDown$foldChange, -resSig27wkDown$baseMean ), ], file="27WkSigMostDownReg_NOCUTOFF.tsv", quote=FALSE, sep="\t", row.names = FALSE )
