## COMP90014 - Assignment 2
# Michael McLellan
# 665968

library(edgeR)
library(knitr)

## Task 1 ##

# Read in data
alldata = read.table("assignment2data_fullcounts.tsv", header=TRUE, sep="\t")
dim(alldata)
colnames(alldata)
alldata[1:3,1:20] # female columns
alldata[1:3,21:40] # male colums
alldata[1:2,41:44] # geneInfo

#remove genes with lots of 0 counts in both male and female replicates
low_counts = apply(alldata[,1:20],1,median)==0 & apply(alldata[,21:40],1,median)==0
filtered = alldata[!low_counts,]

#split data into expression counts and geneInfo
expression = filtered[,1:40]
geneInfo = filtered[,41:44]

#design matrix determining which samples are male 
ismale = c(rep(0,20), rep(1,20))

#create dge object
dge = DGEList(expression, group=ismale)

#normalise
dge = calcNormFactors(dge)

#plot multi-dimensional scaling 
plotMDS(dge, col=c( rep("red",20), rep("blue",20) ))

#DGE analysis, fit variance estimators
dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)
plotBCV(dge)
results = exactTest(dge)

#plot mean variance plot to show how data fits to the model
plotMeanVar(dge, show.raw.vars=TRUE, show.tagwise.vars=TRUE, show.binned.common.disp.vars=FALSE, show.ave.raw.vars=FALSE, dispersion.method="qcml", NBline=TRUE, nbins=100, pch=16 , xlab="Mean Expression (Log10 Scale)", ylab = "Variance (Log10 Scale)")

#get top ten hits
topGI = merge(geneInfo[rownames(topTen),], topTags(results, n=10), by=0, sort=FALSE)
del = c("start_position", "FDR") #Columns not needed in output
topGI[,!(colnames(topGI) %in% del), drop=FALSE]

#get bottom ten hits
low = topTags(results, n=nrow(results$table))$table
low_GI = merge(geneInfo[rownames(low),], tail(low, n=10), by=0, sort=FALSE)
low_GI[,!(colnames(low_GI) %in% del), drop=FALSE]

#Summary of up/down regulaton genes, tagwise results
de = decideTestsDGE(results)
summary(de)
detags = rownames(dge)[as.logical(de)]

### Visualisation ####
#heatmap
topgenes = rownames(topTags(results, n=20))
log.dge = cpm(dge, prior.count=2, log=TRUE)
topgene.log.expression = log.dge[topgenes,]
redtogreen = colorRampPalette(c("red", "yellow", "green"))(n=299)
heatmap(topgene.log.expression, col=redtogreen)

#histogram
hist(results$table$PValue, breaks=150)

#smear plot - shows the relationship between concentraiton and the fold-change across the genes. Differentially expressed genes are colored red and non-DE are black. 
#Blue line is at log-FC of 2, level of biological significance. 
resultsTbl = topTags(results, n=nrow(results$table))$table
de.genes = rownames(resultsTbl)[resultsTbl$PValue <= 0.05]
plotSmear(results, de.tags=de.genes)
abline(h=c(-2, 2), col="blue") 


### Numeric Questions ###
#1. Number of genes with p-values below 0.05, What proportion of the genes tested?
nrow(results$table[results$table$PValue < 0.05, ])
(nrow(results$table[results$table$PValue < 0.05, ]) / nrow(results$table))

#2. Out of all genes on Y chromosome, how many have a P-value below 0.05?
sig = results$table[results$table$PValue <= 0.05, ]
sigGI = geneInfo[rownames(sig),]
sigY = sigGI[sigGI$chromosome_name=="Y", ]
nrow(sigY)

#3. Out of the 100 genes with the lowest P-value, how many are from the X chromosome?
top100 = topTags(results, n=100)
sigX = geneInfo[rownames(top100), ]
sigX.GI = sigX[sigX$chromosome_name=="X", ]
nrow(sigX[sigX$chromosome_name=="X", ])

#4. What is the log-fold change for the gene XIST?
xist = geneInfo[geneInfo$gene_name=="XIST", ]
xistResult = results$table[rownames(results$table)==rownames(xist), ]
round(xistResult$logFC, digits=4)


#### Task 2 ############################################################################################
sub.Data = read.table("assignment2data_subsampled_665968.tsv", header=TRUE, sep="\t")

dim(sub.Data)

#remove genes with 0 counts in male and female
sub.low.counts = apply(sub.Data[,1:20],1,median)==0 & apply(sub.Data[,21:40],1,median)==0
sub.filtered = sub.Data[!sub.low.counts,]

#split data
sub.expression = sub.filtered[,1:40]
sub.geneInfo = sub.filtered[,41:44]

#create dge object, using same group vector from task 1
sub.dge = DGEList(sub.expression, group=ismale)

#normalise
sub.dge = calcNormFactors(sub.dge)

#plot MDS plot
plotMDS(sub.dge, col=c( rep("red",20), rep("blue",20) ))

#fit models
sub.dge = estimateCommonDisp(sub.dge)
sub.dge = estimateTagwiseDisp(sub.dge)
plotBCV(sub.dge)
sub.results = exactTest(sub.dge)

#plot mean variance plot to show how data fits to the model
plotMeanVar(sub.dge, show.raw.vars=TRUE, show.tagwise.vars=TRUE, show.binned.common.disp.vars=FALSE, show.ave.raw.vars=FALSE, dispersion.method="qcml", NBline=TRUE, nbins=100, pch=16 , xlab="Mean Expression (Log10 Scale)", ylab = "Variance (Log10 Scale)")

#get top ten hits
sub.topGI = merge(sub.geneInfo[rownames(sub.topTen),], topTags(sub.results, n=10), by=0, sort=FALSE)
sub.topGI[,!(colnames(sub.topGI) %in% del), drop=FALSE]

#bottom 10 hits
sub.low = topTags(sub.results, n=nrow(sub.results$table))$table
sub.low.GI = merge(sub.geneInfo[rownames(sub.low),], tail(sub.low, n=10), by=0, sort=FALSE)
sub.low.GI[,!(colnames(sub.low.GI) %in% del), drop=FALSE]

#Summary of up/down regulaton genes, tagwise results
sub.de = decideTestsDGE(sub.results)
summary(sub.de)
sub.detags = rownames(sub.dge)[as.logical(sub.de)]

## Visualisation ##

#heatmap
sub.topgenes = rownames(topTags(sub.results, n=20))
sub.log.dge = cpm(sub.dge, prior.count=2, log=TRUE)
sub.topgene.log.expression = sub.log.dge[sub.topgenes,]
heatmap(sub.topgene.log.expression, col = redtogreen)

#histogram
hist(sub.results$table$PValue, breaks=150)

#smear plot - shows the relationship between concentraiton and the fold-change across the genes. Differentially expressed genes are colored red and non-DE are black. 
#Blue line is at log-FC of 2, level of biological significance. 
sub.resultsTbl = topTags(sub.results, n=nrow(sub.results$table))$table
sub.de.genes = rownames(sub.resultsTbl)[sub.resultsTbl$PValue <= 0.05]
plotSmear(sub.results, de.tags=sub.de.genes)
abline(h=c(-2, 2), col="blue") 


### Numerical Questons ###
#1. Number of genes with p-values below 0.05, What proportion of the genes tested?
nrow(sub.results$table[sub.results$table$PValue < 0.05, ])
(nrow(sub.results$table[sub.results$table$PValue < 0.05, ]) / nrow(sub.results$table))

#2. Out of all genes on Y chromosome, how many have a P-value below 0.05?
sub.sig = sub.results$table[sub.results$table$PValue < 0.05, ]
sub.sigGI = sub.geneInfo[rownames(sub.sig),]
sub.sigY = sub.sigGI[sub.sigGI$chromosome_name=="Y", ]
nrow(sub.sigY)

#3. Out of the 100 genes with the lowest P-value, how many are from the X chromosome?
sub.top100 = topTags(sub.results, n=100)
sub.sigX = sub.geneInfo[rownames(sub.top100), ]
nrow(sub.sigX[sub.sigX$chromosome_name=="X", ])

#4. What is the log-fold change for the gene XIST?
subxist = sub.geneInfo[sub.geneInfo$gene_name=="XIST", ]
sub.xistResult = sub.results$table[rownames(sub.results$table)==rownames(xist), ]
round(sub.xistResult$logFC,digits=4)

#### Task 3 ########################################################################
#create vector of random 0 and 1
seed = 1234
set.seed(seed)
randomgroup = sample(c(0,1), 40, replace=TRUE)

#create dge object
rand.dge = DGEList(sub.expression, group=randomgroup)

#normalise
rand.dge = calcNormFactors(rand.dge)

#plot MDS plot
plotMDS(rand.dge, col=c( rep("red",20), rep("blue",20) ))

#fit models
rand.dge = estimateCommonDisp(rand.dge)
rand.dge = estimateTagwiseDisp(rand.dge)
plotBCV(rand.dge)
rand.results = exactTest(rand.dge)

#plot mean variance plot to show how data fits to the model
plotMeanVar(rand.dge, show.raw.vars=TRUE, show.tagwise.vars=TRUE, show.binned.common.disp.vars=FALSE, show.ave.raw.vars=FALSE, dispersion.method="qcml", NBline=TRUE, nbins=100, pch=16 , xlab="Mean Expression (Log10 Scale)", ylab = "Variance (Log10 Scale)")

#get top ten hits
rand.topGI = merge(sub.geneInfo[rownames(rand.topTen),], topTags(rand.results, n=10), by=0, sort=FALSE)
rand.topGI[,!(colnames(rand.topGI) %in% del), drop=FALSE]

#bottom 10 hits
rand.low = topTags(rand.results, n=nrow(rand.results$table))$table
rand.lowGI = merge(sub.geneInfo[rownames(rand.low),], tail(rand.low, n=10), by=0, sort=FALSE)
rand.lowGI[,!(colnames(rand.lowGI) %in% del), drop=FALSE]

#Summary of up/down regulaton genes, tagwise results
rand.de = decideTestsDGE(rand.results)
summary(rand.de)
rand.detags = rownames(rand.dge)[as.logical(rand.de)]

## Visualisation ##

#heatmap
rand.topgenes = rownames(topTags(rand.results, n=20))
rand.log.dge = cpm(rand.dge, prior.count=2, log=TRUE)
rand.topgene.log.expression <- rand.log.dge[rand.topgenes,]
heatmap(rand.topgene.log.expression)

#histogram
hist(rand.results$table$PValue, breaks=150)

#smear plot - shows the relationship between concentraiton and the fold-change across the genes. 
rand.resultsTbl = topTags(rand.results, n=nrow(rand.results$table))$table
rand.de.genes = rownames(rand.resultsTbl)[rand.resultsTbl$PValue <= 0.05]
plotSmear(sub.results, de.tags=rand.de.genes)
abline(h=c(-2, 2), col="blue") 

#1. Number of genes with p-values below 0.05, What proportion of the genes tested?
nrow(rand.results$table[rand.results$table$PValue < 0.05, ])
(nrow(rand.results$table[rand.results$table$PValue < 0.05, ]) / nrow(rand.results$table))

#2. Out of all genes on Y chromosome, how many have a P-value below 0.05?
rand.sig = rand.results$table[rand.results$table$PValue < 0.05, ]
rand.sigGI = sub.geneInfo[rownames(rand.sig),]
rand.sigY = rand.sigGI[rand.sigGI$chromosome_name=="Y", ]
nrow(rand.sigY)

#3. Out of the 100 genes with the lowest P-value, how many are from the X chromosome?
rand.top100 = topTags(sub.results, n=100)
rand.sigX = sub.geneInfo[rownames(rand.top100), ]
nrow(rand.sigX[rand.sigX$chromosome_name=="X", ])

#4. What is the log-fold change for the gene XIST?
rand.xist = sub.geneInfo[sub.geneInfo$gene_name=="XIST", ]
rand.xist.Result = rand.results$table[rownames(rand.results$table)==rownames(xist), ]
rand.xist.Result$logFC

#Session info
sessionInfo()

