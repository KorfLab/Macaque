setwd("C:/Dropbox/Macaque_Project/miRNA_analysis/Differential_expression/BLAST_against_Macaca_only")
myPath = file.path("C:/Dropbox/Macaque_Project/miRNA_analysis/Differential_expression/BLAST_against_Macaca_only/miRNA_Expression_Macaca_only_db.tsv")
counts = read.table(myPath, header=TRUE, row.names=1)
metadata = data.frame(
	row.names = colnames(counts),
	condition = c("treated", "treated", "treated", "treated", "treated", "NA", "NA", "NA", "NA", "untreated"),
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

# Examine differential expression between both 2wk and 27wk infected vs non-infected jejunum 
res = nbinomTest(cds, "untreated", "treated")
#head(res)
#plotMA(res)
#hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="") # Prints histogram of p-values
#hist(res$padj, breaks=100, col="skyblue", border="slateblue", main="") # Prints histogram of adjusted p-values
resSig = res[res$padj < 0.1,]
resSig = na.omit(resSig)
head(resSig[order(resSig$pval),]) # Most significantly differentially expressed genes
head(resSig[order(resSig$foldChange, -resSig$baseMean ), ] ) # Most strongly down-regulated significant genes
head(resSig[order(-resSig$foldChange, -resSig$baseMean ), ] ) # Most strongly up-regulated significant genes
write.table(resSig[order(resSig$pval), ], file="2Wk+27WkSigDifExp.tsv", quote=FALSE, sep="\t" )
write.table(resSig[order(resSig$foldChange, -resSig$baseMean ), ], file="2Wk+27WkSigMostDownReg.tsv", quote=FALSE, sep="\t" )
write.table(resSig[order(-resSig$foldChange, -resSig$baseMean ), ], file="2Wk+27WkSigMostUpReg.tsv", quote=FALSE, sep="\t" )
