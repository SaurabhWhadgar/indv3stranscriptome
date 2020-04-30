library(pheatmap)
library(stats)
library(corrplot)
library(knitr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(stats)
library(ggplot2)
library(devtools)
library(ggfortify)
library(reshape)
library(wesanderson)

#Reading File and making matrix dataframe
cdd <- as.data.frame(read.csv("63_04_05_06_07_08_09_10_15_16_41.count.extra.gene",header = F,sep="\t"))
cdd<-cdd[,-c(1,2,3,6,7,8)]
cdd <- cdd[c(3,4,5,6,7,8,9,10,11,12,13,14,2,1)]
cdd <-  cdd[c(1,11,12,2,3,4,5,6,7,8,9,10,13,14)]
cdd$V21 <- cdd$V5-cdd$V4
cdd <- cdd[c(1,2,3,4,5,6,7,8,9,10,11,12,15)]
colnames(cdd) <- c("gene","EM0_1","EM2_4","EM4_8","EM8_12","Larva","Pupa","AdFe","AdMa","PBM","PEM","FeCa","length")


df <- data.frame(colSums(subset(cdd, select = c(-length,-gene))))

colnames(df)<-c("ReadCount")
df$Sample <- rownames(df)
options("scipen"=200000)
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
head(rpkm)
colnames(rpkm) <- c("EM0_1","EM2_4","EM4_8","EM8_12","Larva","Pupa","AdFe","AdMa","PBM","PEM","FeCa")
summary(rpkm)


df<-data.frame(as.integer(colSums(rpkm)))
head(df)
colnames(df)<-c("Length")
df$Sample <- c("EM0_1","EM2_4","EM4_8","EM8_12","Larva","Pupa","AdFe","AdMa","PBM","PEM","FeCa")
f <- data.frame(df)
p<-ggplot(data=df, aes(x=Sample, y=ReadCount,fill=Sample)) +
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


### Findig Std.Devaiation ( Variance )
V <- apply(rpkm, 1, var)
rownames(rpkm)<-cdd$gene
  head(V)
V <- apply(tpm, 1, var)

#sort the results by variance in decreasing order and select the top 100 genes (ROWNUMBER)
selectedGenes <- names(V[order(V, decreasing = T)][1:100])

head(selectedGenes)
### Pheatmap
colData <- read.table('coldata_file', header = T, sep = '\t', stringsAsFactors = TRUE)
rownames(colData) <- colnames(rpkm)
rownames(colData) <- colnames(tpm)
write.table(rpkm,"RPKM_FILE",sep = "\t",quote = F)

pheatmap(rpkm[selectedGenes,], scale = 'row', 
         show_rownames = T, 
         fontsize_row = 6,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = colData)

### PCA

# #transpose the matrix 
# M <- t(rpkm[selectedGenes,])
# transform the counts to log2 scale 
# M <- log2(M + 1)
# compute PCA 
# pcaResults <- prcomp(M)
# pcaResults$group<-c("EM0_1","EM2_4","EM4_8","EM8_12","Larva","Pupa","AdFe","AdMa","PBM","PEM","FeCa")
# autoplot(pcaResults,fill=pcaResults$group)
# head(pcaResults,fill = 'group')

## Corrlatinal Matrix between Samples
correlationMatrix <- cor(rpkm)
kable(correlationMatrix,booktabs = TRUE)

# The correlation plot order by the results of the hierarchical clustering/ PHEATMAP
corrplot(correlationMatrix, order = 'hclust', addCoef.col = 'grey50',title = 'Correlation Matrix Between Sample')
pheatmap(correlationMatrix,  annotation_col = colData, main = "Pairwise Correlation Matrix Between Sample")


rownames(tpm)<-cdd$gene
rownames(rpkm)<-cdd$gene
tpm <- data.frame(tpm)
rpkm <- data.frame(rpkm)


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

df <- rpkm[order(-rpkm$EM0_1,-rpkm$Larva),]
head(df[1:100,])
draw_rownames <- function(.data) .data %>% do(mutate(.,rownames=rownames(.)))
head(rpkm)
df <- rpkm %>% draw_rownames() %>% filter(AdFe> 200, AdMa < 50)
rownames(df)<-df$rownames
df
pheatmap(df[1:29,1:11], scale = 'row', fontsize_row = 6,
         show_rownames = T, 
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = colData,
         main = "Differentially Expressed Genes in AdFe")



geneofint <- rpkm[c(
             "ID=g3704","ID=g5128",
             "ID=g6054","ID=g6842",
             "ID=g6843","ID=g6908",
             "ID=g24450","ID=g18313",
             "ID=g12024","ID=g12863",
             "ID=g16696","ID=g18569")
             ,]

geneofint$Gene <- c("Spo11","Cardinal","ACE1","CYP6P9a/b","CYP6P4",
                     "RAD51","MRE11","KH","RAB5","CYP6M10a","KDR","FREP1")
 


Molten <- melt(geneofint, id.vars = "Gene")
ggplot(Molten, aes(x = variable, y = value, colour = Gene, group=Gene)) + 
  geom_line() + 
  facet_wrap(~ Gene)+
  ggtitle("Genes of Interest Expression : RPKM ")+
  labs(y="RPKM Expression", x = "Stage")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#################### GOI Separetly ###################

rpkm[c("ID=g3704"),1:11]
gene1<-rpkm[c("ID=g3704"),1:11]
gene1$Gene <- c("SPO11")
gene1
Molten <- melt(gene1, id.vars = "Gene")
ggplot(Molten, aes(x = variable, y = value, colour = Gene, group=Gene)) + 
  geom_line() + 
  facet_wrap(~ Gene)+
  ggtitle("Genes of Intrest Expression : SPO11 ")+
  labs(y="RPKM", x = "Stage")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

############## Finding Novel gene using SD approch in progress.... #####################
head(rpkm)
rpkm$sd = apply(rpkm, 1, sd)
head(rpkm)

geneofint_sd <-  rpkm[c("ID=g6842","ID=g6843","ID=g12863"),]

geneofint$Gene <- c("Spo11","Cardinal","ACE1","CYP6P9a/b","CYP6P4",
                    "RAD51","MRE11","KH","RAB5","CYP6M10a","KDR","FREP1")

head(rpkm)
pbm <- rpkm[order(-rpkm$PBM, rpkm$sd), ]
rpkm %>% filter(Larva > 4000, EM0_1<2000)


draw_rownames <- function(.data) .data %>% do(mutate(.,rownames=rownames(.)))
df12<-rpkm %>% draw_rownames() %>% filter(Larva > 800, EM0_1 < 800)
rownames(df12)<-df12$rownames
df12<-rapply(df12[,1:11], f = function(x) replace(x, is.infinite(x), 0), classes = "numeric", how = "replace")
round(df12[,1:11])
write.table(df12,"L800_EM01_800",sep="\t",quote = F,row.names = T,col.names = T)
  pheatmap(round(df12[,1:11]),scale = 'row', fontsize_row = 8,
         show_rownames = T, 
         cluster_rows = F,
         cluster_cols = F,
         main = "Differential Expressed Gene Larva")

head(df12)
melted_cormat <- melt(df12)
head(melted_cormat)
pal <- wes_palette("Zissou1", 10, type = "continuous")
ggplot(data = melted_cormat, aes(x=variable, y=X, fill=value)) + 
  geom_tile()+
  scale_fill_gradient(low = "white",  high = "steelblue")+
theme(axis.text.x = element_text(angle = 90, hjust = 1))

set.seed(42)

head(data.frame(apply(rpkm, 1, sd)))

data.frame(rowVars(as.matrix(rpkm$EM0_1), suma = NULL, std = TRUE))


rpkm$var <- data.frame(apply(rpkm, 1, var))

rpkm[,1:11] %>%  summarise_all(funs(sd(., na.rm=TRUE)))

sort(rpkm$var)

rm <- rowMeans(rpkm, na.rm=TRUE)
rv <- apply(rpkm, MARGIN=1, FUN=var, na.rm=TRUE)


rpkm <- rpkm[order(-rpkm$EM0_1,rpkm$Larva),]
pheatmap(rpkm[1:100,1:11],scale = 'row', fontsize_row = 8,
         show_rownames = T, 
         cluster_rows = F,
         cluster_cols = F,
         main = "Differential Expressed Gene Larva")


