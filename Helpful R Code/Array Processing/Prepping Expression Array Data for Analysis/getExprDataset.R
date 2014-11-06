rm(list=ls())
setwd("/data/kiemele/LXS/LXSpaper2014/DID")

#############
# Load Data #
#############
### load affy data Laura made where probesets are clustered by parental RNA-Seq data ### 
load("/data/Tabastore3/LauraS/LXS/RNA-Seq/totalRNA.24Oct13/arrayMask/Rdata/collapsed.18Jun14.Rdata")
ls()
dim(collapsed)
collapsed[1:5, 1:5]

#############################
# Check For Outlier Samples #
#############################
rm(list=ls())
library(WGCNA)
options(stringsAsFactors = FALSE)
setwd("/data/kiemele/LXS/LXSpaper2014/DID")
load("heritExprs.Rdata")

datExpr_t <- t(heritExprs)
sampleTree = flashClust(as.dist(1-cor(heritExprs)), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#get a better look at LXS 110
dim(heritExprs)
strain = sapply(strsplit(colnames(heritExprs), split="_", fixed=TRUE), "[[", 2)
test = heritExprs[,which(strain=="LXS110")]
test = heritExprs[,which(strain=="LXS52")]

## BAD samples;
badSamples = c("M0481_LXS110_run16.CEL", "M0483_LXS110_run16.CEL", "M0103_LXS52_Run4.CEL")

sig = c()
for(i in badSamples){
	sig = c(sig, which(colnames(heritExprs)==i))
}
length(sig)
cleanedExprs = heritExprs[,-sig]
save(cleanedExprs, file="cleanedExprsIndividualMouse.Rdata")

#######################
# Heritability Filter #
#######################
### Get the heritability of the clusters and filter only on those with heritability > 33%;
strain = sapply(strsplit(colnames(collapsed), split="_", fixed=TRUE), "[[", 2)

test = apply(collapsed, 1, function(a) summary(lm(a~as.factor(strain)))$r.squared)

herits = cbind(rownames(collapsed), test)
colnames(herits) = c("transcript", "heritability")
write.table(herits, file="herits.LSclusters.LXS.brain.txt", sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

heritable <- herits[herits[,2]>0.33,1]
heritablePS <- match(heritable,rownames(collapsed))
heritablePS <- heritablePS[!is.na(heritablePS)]
heritExprs <- collapsed[heritablePS,]
save(heritExprs, file="heritExprs.Rdata")

#############################
# Check For Outlier Samples #
#############################
rm(list=ls())
library(WGCNA)
options(stringsAsFactors = FALSE)
setwd("/data/kiemele/LXS/LXSpaper2014/DID")
load("heritExprs.Rdata")

datExpr_t <- t(heritExprs)
sampleTree = flashClust(as.dist(1-cor(heritExprs)), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#get a better look at LXS 110
dim(heritExprs)
strain = sapply(strsplit(colnames(heritExprs), split="_", fixed=TRUE), "[[", 2)
test = heritExprs[,which(strain=="LXS110")]
test = heritExprs[,which(strain=="LXS52")]

## BAD samples;
badSamples = c("M0481_LXS110_run16.CEL", "M0483_LXS110_run16.CEL", "M0103_LXS52_Run4.CEL")

sig = c()
for(i in badSamples){
	sig = c(sig, which(colnames(heritExprs)==i))
}
length(sig)
cleanedExprs = heritExprs[,-sig]
save(cleanedExprs, file="cleanedExprsIndividualMouse.Rdata")


#####################
#  Get Strain Means #
#####################
rm(list=ls())
load("cleanedExprsIndividualMouse.Rdata")

#first determine what array belongs to what strain;
strain = sapply(strsplit(colnames(cleanedExprs), split="_", fixed=TRUE), "[[", 2)

means = t(apply(cleanedExprs, 1, function(a) lm(a~as.factor(strain) -1)$coefficients))
colnames(means) = sapply(strsplit(colnames(means), split=")", fixed=TRUE), "[[", 2)

save(means, file="datExprStrainMeans.RData")

#############################
# Check for Outlier Strains #
#############################
sampleTree = flashClust(as.dist(1-cor(means)), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "strainClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Strain clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
dev.off()
