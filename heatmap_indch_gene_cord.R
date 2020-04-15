library(pheatmap)
library(stats)
library(corrplot)
library(knitr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

# Reading File and making matrix dataframe
cdd <- as.data.frame(read.csv("IndV3s_gene_cord.tsv",header = F,sep="\t"))
cdd <- cdd[c(3,4,5,6,7,8,9,10,11,12,13,14,2,1)]
cdd$V15 <- cdd$V2-cdd$V1
cdd <- cdd[c(1,2,3,4,5,6,7,8,9,10,11,12,15)]
colnames(cdd) <- c("gene","EM0_1","EM2_4","EM4_8","EM8_12","Larva","Pupa","AdFe","AdMa","PBM","PEM","FeCa","length")

df <- data.frame(colSums(subset(cdd, select = c(-length,-gene))))
colnames(df)<-c("ReadCount")
df$Sample <- rownames(df)

p<-ggplot(data=df, aes(x=Sample, y=ReadCount,fill=Sample)) +
  geom_bar(stat="identity")+
  ggtitle("Raw Read Count")+
  theme_minimal()
p



# Taking gene lenghts as vector
geneLengths <- as.vector(subset(cdd, select = c(length)))
head(geneLengths)

# compute RPKM
rpkm <- apply(X = subset(cdd, select = c(-length,-gene)), 
              MARGIN = 2, 
              FUN = function(x) 10^9 * x / geneLengths / sum(as.numeric(x)))
# Converting RPKM as dataframe
rpkm <- data.frame(rpkm)

# Changing colname of RPKM
colnames(rpkm) <- c("EM0_1","EM2_4","EM4_8","EM8_12","Larva","Pupa","AdFe","AdMa","PBM","PEM","FeCa")

#Summary of RPKM
summary(rpkm)
ncol(rpkm)
rpkm_int<-as.data.frame(lapply(rpkm[,1:11], as.integer))

head(rpkm)
# Bar plot of colsums RPKM
#barplot(log2(as.integer(colSums(rpkm))))
#head(as.integer(colSums(rpkm)))
df<-data.frame(as.integer(colSums(rpkm)))
head(df)
colnames(df)<-c("Length")
df$Sample <- c("EM0_1","EM2_4","EM4_8","EM8_12","Larva","Pupa","AdFe","AdMa","PBM","PEM","FeCa")
df <- data.frame(df)
p<-ggplot(data=df, aes(x=Sample, y=Length,fill=Sample)) +
  geom_bar(stat="identity")+
  ggtitle("RPKM Normalized Read Count")+
  theme_minimal()
p

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### TPM not RPKM : Calculationg RPK value
rpk <- apply( subset(cdd, select = c(-length,-gene)), 
              2, 
              function(x) x/(geneLengths/1000))
#normalize by the sample size using rpk values
tpm <- apply(data.frame(rpk), 2, function(x) x / sum(as.numeric(x)) * 10^6)
head(tpm)
colnames(tpm)<-c("EM0_1","EM2_4","EM4_8","EM8_12","Larva","Pupa","AdFe","AdMa","PBM","PEM","FeCa")
#Plotting colusum Library Size
df<-data.frame(as.integer(colSums(tpm)))
colnames(df)<-c("Length")
df$Sample <- c("EM0_1","EM2_4","EM4_8","EM8_12","Larva","Pupa","AdFe","AdMa","PBM","PEM","FeCa")
head(df)

p<-ggplot(data=df, aes(x=Sample, y=Length,fill=Sample)) +
  geom_bar(stat="identity")+
  ggtitle("TPM Normalized Read Count")+
  theme_minimal()
p

##################### UPTO THIS ##### TPM and RPKM Is DONE #############
 #Place for boxplot

dat <- cdd[,2:12]
dat.m <- melt(dat, measure.vars=c("EM0_1","EM2_4","EM4_8","EM8_12","Larva","Pupa","AdFe","AdMa","PBM","PEM","FeCa"))

# dat.m$value<-log2(dat.m$value)
# stack(cdd)
# df <- dat.m %>%
#   ggplot(aes(x=variable, y=value, fill=variable, alpha=0.1)) + 
#   geom_boxplot() +
#   theme_ipsum() +
#   theme(legend.position = "left") +
#   xlab("")
# df

############################# HEATMAP PLOTTING #############################################

### Findig Std.Devaiation ( Variance )
rownames(rpkm)<-cdd$gene
V <- apply(rpkm, 1, var)

head(V)
rownames(tpm)<-cdd$gene
V1 <- apply(tpm, 1, var)

#sort the results by variance in decreasing order and select the top 100 genes (ROWNUMBER)
selectedGenes <- names(V1[order(V1, decreasing = T)][1:100])

head(selectedGenes)
### Pheatmap
colData <- read.table('coldata_file', header = T, sep = '\t', stringsAsFactors = TRUE)
rownames(colData) <- colnames(rpkm)
rownames(colData) <- colnames(tpm)
#write.table(rpkm,"RPKM_FILE",sep = "\t",quote = F)

pheatmap(tpm[selectedGenes,], scale = 'row', 
         show_rownames = T, 
         fontsize_row = 6,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = colData)

## Corrlatinal Matrix between Samples
correlationMatrix <- cor(tpm)
kable(correlationMatrix,booktabs = TRUE)

# The correlation plot order by the results of the hierarchical clustering/ PHEATMAP
corrplot(correlationMatrix, order = 'hclust', addCoef.col = 'grey50',title = 'Correlation Matrix Between Sample')
pheatmap(correlationMatrix,  annotation_col = colData, cutree_cols = 2)


################## MAKE IT AS FUNCTION ########################


## $$$ WHAT ma'am WANTS for PAPER $$$$ ###############

####################### TPM For Furhter Analysis #################### 

rownames(tpm)<-cdd$gene
rownames(rpkm)<-cdd$gene
tpm <- data.frame(tpm)
rpkm <- data.frame(rpkm)

head(tpm[order(-tpm$EM0_1),])
head(rpkm[order(-rpkm$EM0_1),])

countsTable<-data.frame(tpm[1:100,])
typeof(data.frame(countsTable))
countsTable[!is.finite(countsTable)] <- NA
colMeans(countsTable$V2, na.rm=TRUE)
countsTable<-countsTable[is.finite(rowSums(countsTable)),]

tpm_em01 <- tpm[order(tpm$EM0_1),]
head(tpm_em01)
tpm_em10 <- as.integer(tpm[order(-tpm$EM0_1,tpm$AdMa),])

length(tpm)
df <- data.frame(unlist(tpm),nrow=length(tpm), byrow=T)
tpm2 <- as.data.frame(tpm)
as.integer(as.data.frame(tpm2))


################## RPKM Heatmap ##############
colData <- read.table('coldata_file', header = T, sep = '\t', stringsAsFactors = TRUE)
rownames(colData) <- colnames(rpkm)
rownames(rpkm)<-cdd$gene
#rpkm <- data.frame(rpkm)

head(rpkm[order(-rpkm$EM0_1),])

countsTable<-data.frame(rpkm[1:100,])
countsTable<-countsTable[is.finite(rowSums(countsTable)),]

rpkm_em01 <- rpkm[order(rpkm$EM0_1),]
head(rpkm_em01)

df <- rpkm[order(-rpkm$EM4_8),]
head(df[1:100,])
pheatmap(df[1:100,], scale = 'row', fontsize_row = 6,
         show_rownames = T, 
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = colData,
         main = "TPM RC SORT DSEC EM 4-8hr")

  dev.off()

  geneofint<-NULL
 geneofint <-  rpkm[c("ID=g3704",
             "ID=g6054","ID=g6842",
             "ID=g6843","ID=g6908",
             "ID=g5128","ID=g18313",
             "ID=g18569","ID=g8216",
             "ID=g16696","ID=g12024",
             "ID=g12863"),]

 geneofint$Gene <- c("Spo11","Ace1","CYP6P9a/b","CYP6P4","Rad51","Cardinal","KH","FREP1","Mre11","KDR","Rab5","CYP6M10a")
 
 

Molten <- melt(geneofint, id.vars = "Gene")
head(Molten)
ggplot(Molten, aes(x = variable, y = value, colour = Gene, group=Gene)) + 
  geom_line() + 
  facet_wrap(~ Gene)+
  ggtitle("Genes of Intrest Expression : RPKM ")+
  labs(y="RPKM Expression", x = "Stage")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

