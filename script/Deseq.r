#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq")
library(DESeq)
#setwd(“path/to/dir/containing/read/counts”)
countsTable <- read.csv('uniq_gname_read_count.tsv',header=TRUE,sep="\t")
#countsTable$miRNA = paste(countsTable$miRNA,1:nrow(countsTable),sep ="_")
rownames(countsTable) <- countsTable$gene_name	
countsTable <- countsTable[,-1]
head(countsTable)
#conditions <- factor(c("normal","normal","cancer","cancer"))
#OR
conditions <- factor(c(rep("control",3),rep("treatment",3)))
countDataSet <- newCountDataSet(countsTable,conditions)
countDataSet <-estimateSizeFactors(countDataSet)
sizeFactors(countDataSet)
head(counts(countDataSet))
head(counts(countDataSet,normalized=TRUE))
countDataSet <-estimateDispersions(countDataSet)
###### for no replicates :-
#countDataSet <-estimateDispersions(countDataSet, method="blind",sharingMode="fit-only")


pdf("CvsE_disp_DESeq.pdf",width=5,height=5)
plotDispEsts(countDataSet)
dev.off()


DEVal=nbinomTest(countDataSet,"control","treatment")
pdf("CvsE_plotMA_DESeq.pdf",width=5,height=5)
plotMA(DEVal,col = ifelse(DEVal$pval<=0.05,"gray32", "red"), linecol = "red")
dev.off()

#######To filter significant genes according to dalse discovery rate(FDR)
#DEValSig = DEVal [DEVal$padj < 0.1,]
#head (DEValSig [order (DEValSig$pval),])
#write.csv(DEValSig, file = "D12M12_padj.csv")

###### To filter most significantly differentially expressed genes
#diff = (DEValSig [order(DEValSig$pval),])
#head (diff)
#write.csv(diff, file = "D12M12_diffexp.csv")

####### To filter downregulated genes
#down = ( DEValSig [order(DEValSig$foldChange,-DEValSig$baseMean),])
#head (down)
#write.csv(down, file = "D12M12_downreg.csv")

###### To filter upregulated genes
#upreg = (DEValSig [order(-DEValSig$foldChange,-DEValSig$baseMean),])
#head (upreg)
#write.csv(upreg, file = "D12M12_upreg.csv")


#png("D1234_hist_pVal_DESeq.png")
#hist(DEVal$pval, breaks=100, col="skyblue",border="slateblue")
#dev.off()

#### To plot heatmap 
vsdFull = varianceStabilizingTransformation (countDataSet)
library ("RColorBrewer")
library ("gplots")
select = order (rowMeans(counts(countDataSet)), decreasing=TRUE) [1:30]   #### for first 30 genes in list
hmcol =colorRampPalette(brewer.pal(9,"GnBu"))(100)
pdf("CvsE_heatmap_first30.DESeq.pdf",width=5,height=10)
heatmap.2(exprs(vsdFull)[select,],col=hmcol,trace="none",margin=c(10,6))
dev.off()
 #### Heatmap for read counts to confirm the result
#png("BCL_heatma_forreadcounts_first30mirnas_DESeq.png")
#heatmap.2(counts(countDataSet)[select,],col=hmcol,trace="none",margin=c(10,6))
#dev.off()

###### Heatmap to show Euclidean distances between samples
#dists = dist( t(exprs(vsdFull)))
#mat=as.matrix (dists)
#rownames (mat) =colnames(mat)=with(pData(countDataSet), paste (condition, sep= ":"))
#png("BCL_distancebetsamples_DESeq.png")
#heatmap.2 (mat, trace ="none", col =rev (hmcol), margin =c (13,13))
#dev.off()

###### PCA plot of samples
pdf("CvsE_PCAplot_DESeq.pdf",width=5,height=5)
print (plotPCA(vsdFull, intgroup=c("condition")))
dev.off ()


write.csv(DEVal,file="CvsE_DEGs.csv")


#Plot heat map
#library(gplots)
#library(pheatmap)
#data=read.table("SRP002628_DEGs_top100_RPKM.txt",sep="\t",header=T)
#data_filtered <- as.matrix(data.frame(data$C02,data$C03,data$C11,data$C13,data$C15,data$C19,data$C23,data$N02,data$N03,data$N11,data$N13,data$N15,data$N19,data$N23))colnames(data_filtered) <- c("C02","C03","C11","C13","C15","C19","C23","N02","N03","N11","N13","N15","N19","N23")
#rownames(data_filtered) <- data$Gene_ID


#tiff("DEG_heatmap.tiff",height=10, width=8,units='in',res=500,compression = c ("lzw"))

#pheatmap(data_filtered ,scale = "none",cluster_cols=TRUE,color = colorRampPalette(c("green","red")) (100),margins=c(8,8),fontsize = 10)

#dev.off()


#########plotDispEsts <- function( countDataSet ){plot(rowMeans( counts( cds, normalized=TRUE ) ), fitInfo(countDataSet)$perGeneDispEsts, pch = ‘.’, log=”xy”, ylab=”dispersion”,xlab=”mean of normalized counts”)xg = 10^seq( -.5, 5, length.out=300 )lines( xg, fitInfo(countDataSet)$dispFun( xg ), col=”red” )}

##########plotDispEsts <- function(countDataSet){plot(rowMeans(counts(countDataSet,normalized=TRUE)),fitInfo(countDataSet)$perGeneDispEsts,pch = '.', log="xy", ylab="dispersion", xlab="mean of normalized counts",xg = 10^seq(-.5, 5, length.out=300 ), lines( xg, fitInfo(countDataSet)$dispFun( xg ), col="red" ))}







