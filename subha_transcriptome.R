library(DESeq)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(plotly)
print("Github Update file")
#setwd(“path/to/dir/containing/read/counts”)
countsTable <- read.csv('gene_matrxi.tsv',header=F,sep="\t")
#countsTable$miRNA = paste(countsTable$miRNA,1:nrow(countsTable),sep ="_")
ncol(countsTable)
colnames(countsTable)<- c("Gene","316","341","863","304","305","306","307","308","309","310","315")
rownames(countsTable) <- countsTable$Gene
countsTable <- countsTable[,-1]
head(countsTable)
conditions <- c('0-1','2-4','4-8','8-12','Larva','Pupa','AdFe','AdM','PBM','OPE','FeCa')
countDataSet <- newCountDataSet(countsTable,conditions)
countDataSet <-estimateSizeFactors(countDataSet)
sizeFactors(countDataSet)
head(counts(countDataSet))
normaliezed <- log2(counts(countDataSet,normalized=TRUE))
head(normaliezed)
countDataSet <-estimateDispersions(countDataSet,method='blind',sharingMode='fit-only')
#cds<-estimateDispersions(cds, method='blind',sharingMode='fit-only')
boxplot(log2(countsTable))
vsd <- vst(countDataSet, blind=FALSE)


countsTable %>%
  ggplot(aes(x=rownames(countsTable), y=countsTable, fill=countsTable)) +
  geom_boxplot() 
head(countsTable)

ncol(countsTable)
countsTable<-data.frame(log2(countsTable[1:100,2:12]))
countsTable$V2[!is.finite(countsTable$V2)] <- NA
colMeans(countsTable$V2, na.rm=TRUE)
countsTable<-countsTable[is.finite(rowSums(countsTable)),]


library(reshape2)
dat.m <- melt(countsTable, measure.vars=c('V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12'))
head(dat.m)
dat.m$value<-log2(dat.m$value)
ggplot(dat.m) + 
  geom_boxplot(aes(x=variable, y=value, fill=value)) + 
  geom_point() +
  labs(title="Log Transformed read Count",x="Sample", y = "Counts") +
  geom_jitter(color="black", size=0.4, alpha=0.9) 


ggplot(dat.m) +
  geom_boxplot(aes(x=variable, y=value, fill=value)) + 
 # geom_jitter(color="black", size=0.4, alpha=0.9) +
  geom_point(aes(x=value, y=value)) 

df <- dat.m %>%
ggplot(aes(x=variable, y=value, fill="orange", alpha=0.1)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("red", "grey")) +
  theme_ipsum() +
  theme(legend.position = "left") +
  xlab("")
  #geom_line(aes(group = value), colour = "grey")
df
df + geom_line(aes(group = value), colour = "grey")

boysbox <- countsTable + geom_boxplot()
boysbox + geom_line(aes(group = Subject), colour = "blue")

head(stack(countsTable))

ggplot(stack(countsTable), aes(x = ind, y = values)) +
  geom_boxplot()+ geom_line(aes(group = values), colour = "blue")

stack(countsTable)
ggplot(stack(countsTable), aes(ind, values)) +
  geom_boxplot(aes(fill = ind)) +
  theme_classic()+ 
  geom_line(aes(group = values), colour = "blue")
 # geom_dotplot(binaxis = 'y', dotsize = 0.1, stackdir = 'center')

dat.m + geom_line(aes(group = variable))
ggplot(melt(df), aes(variable, value)) + geom_boxplot()
dat.m %>%
  ggplot( aes(x=variable, y=value, fill=value)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle("A boxplot with jitter") +
  xlab("")

ggplot(gapminder,aes(x=continent, y=lifeExp, fill=continent)) +
  geom_boxplot() +
  geom_jitter(width=0.25, alpha=0.5)