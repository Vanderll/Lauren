options(stringsAsFactors=FALSE)

origData <- read.table(file="rma.fullPS.HXB_BXH.brain.PhenoGen.txt",sep="\t",header=TRUE,row.names=1)
calls <- read.table(file="dabg.fullPS.HXB_BXH.brain.PhenoGen.txt",sep="\t",header=TRUE,row.names=1)

table(rownames(origData)==rownames(calls))

print(paste(nrow(origData),"probe sets on the array"))
print(paste(ncol(origData),"arrays"))

#729859 probe sets on the array, 108 arrays

#####################################################################
##  limit to probe sets that have >1% of samples above background  ##
#####################################################################

presentExprs <- origData[rowSums(calls<0.0001)>(ncol(origData)*0.01),]
presentCalls <- calls[rowSums(calls<0.0001)>(ncol(origData)*0.01),]
dim(presentExprs)
save(presentExprs,file="presentExprs.Rdata")

print(paste(nrow(presentExprs),"probe sets remained after DABG filter"))

#reduced to 206,552 probesets

#TEST 
load("presentExprs.Rdata")


#########################################################
##  limit to probe sets whose expression is heritable  ##
#########################################################

herits <- read.table(file="herits.fullPS.HXB_BXH.brain.txt",sep="\t",header=TRUE)
heritable <- herits[herits[,2]>0.33,1]
heritablePS <- match(heritable,rownames(presentExprs))
heritablePS <- heritablePS[!is.na(heritablePS)]
heritExprs <- presentExprs[heritablePS,]
save(heritExprs,file="heritExprs.Rdata")

print(paste(nrow(heritExprs),"probe sets remained after heritability filter"))
#85,048

##################################
## Clean the Data & Take Means  ##
##################################
rm(list=ls())
library(WGCNA)
options(stringsAsFactors = FALSE)

setwd("C:/Users/kiemelel/Documents/HXB/Brain")
load("heritExprs.Rdata")
objects()
dim(heritExprs)

#######################
##   Clean Samples   ##
#######################
datExpr_t <- t(heritExprs)
##Now Check for Outliers###;
#sampleTree = flashClust(dist(t(heritExprs)), method = "average");
sampleTree = flashClust(as.dist(1-cor(heritExprs)), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)

##Remove BAD Sample(s): PD_8604_02_Brain

Sig_list <- c()
Sig_list = c(Sig_list, colnames(heritExprs)=="PD_8604_02_Brain")

heritExprs[1:5, 85:87]

heritExprs_sampleclean <- heritExprs[,-86]
dim(heritExprs_sampleclean)

samples <- colnames(heritExprs_sampleclean)
write.csv(samples, file="clean_samples_HXBbrain.csv")
##########################
##   Get Strain Means   ##
##########################
BNLX_dat = heritExprs[, c("BN.LX_01_Brain", "BN.LX_02_Brain", "BN.LX_03_Brain", "BN.LX_04_Brain")]
BNLX <- as.matrix(apply(BNLX_dat, 1, mean))

BXH08_dat = heritExprs[, c("X08BXH_01_Brain", "X08BXH_03_Brain")]
BXH08 = as.matrix(apply(BXH08_dat, 1, mean))

BXH10_dat = heritExprs[, c("X10BXH_01_Brain", "X10BXH_02_Brain", "X10BXH_03_Brain", "X10BXH_04_Brain")]
BXH10 <- as.matrix(apply(BXH10_dat, 1, mean))

BXH11_dat = heritExprs[, c("X11BXH_01_Brain", "X11BXH_02_Brain", "X11BXH_03_Brain", "X11BXH_04_Brain")]
BXH11 <- as.matrix(apply(BXH11_dat, 1, mean))

BXH13_dat = heritExprs[, c("X13BXH_01_Brain", "X13BXH_02_Brain", "X13BXH_03_Brain", "X13BXH_04_Brain")]
BXH13 <- as.matrix(apply(BXH13_dat, 1, mean))

BXH6_dat = heritExprs[, c("X6BXH_01_Brain", "X6BXH_02_Brain", "X6BXH_03_Brain", "X6BXH_04_Brain")]
BXH6 <- as.matrix(apply(BXH6_dat, 1, mean))

HXB1_dat = heritExprs[, c("X1HXB_01_Brain", "X1HXB_02_Brain", "X1HXB_04_Brain")]
HXB1 <- as.matrix(apply(HXB1_dat, 1, mean))

HXB10_dat = heritExprs[, c("X10HXB_01_Brain", "X10HXB_02_Brain", "X10HXB_03_Brain", "X10HXB_04_Brain")]
HXB10 <- as.matrix(apply(HXB10_dat, 1, mean))

HXB13_dat = heritExprs[, c("X13HXB_01_Brain", "X13HXB_02_Brain", "X13HXB_03_Brain", "X13HXB_04_Brain")]
HXB13 <- as.matrix(apply(HXB13_dat, 1, mean))

HXB15_dat = heritExprs[, c("X15HXB_01_Brain", "X15HXB_02_Brain", "X15HXB_03_Brain", "X15HXB_04_Brain")]
HXB15 <- as.matrix(apply(HXB15_dat, 1, mean))

HXB17_dat = heritExprs[, c("X17HXB_01_Brain", "X17HXB_02_Brain", "X17HXB_03_Brain", "X17HXB_04_Brain")]
HXB17 <- as.matrix(apply(HXB17_dat, 1, mean))

HXB18_dat = heritExprs[, c("X18HXB_01_Brain", "X18HXB_02_Brain", "X18HXB_03_Brain", "X18HXB_04_Brain")]
HXB18 <- as.matrix(apply(HXB18_dat, 1, mean))

HXB2_dat = heritExprs[, c("X2HXB_01_Brain", "X2HXB_02_Brain", "X2HXB_03_Brain", "X2HXB_04_Brain")]
HXB2 <- as.matrix(apply(HXB2_dat, 1, mean))

HXB23_dat = heritExprs[, c("X23HXB_01_Brain", "X23HXB_02_Brain", "X23HXB_03_Brain", "X23HXB_04_Brain")]
HXB23 <- as.matrix(apply(HXB23_dat, 1, mean))

HXB25_dat = heritExprs[, c("X25HXB_01_Brain", "X25HXB_02_Brain", "X25HXB_03_Brain", "X25HXB_04_Brain")]
HXB25 <- as.matrix(apply(HXB25_dat, 1, mean))

HXB26_dat = heritExprs[, c("X26HXB_01_Brain", "X26HXB_02_Brain", "X26HXB_03_Brain", "X26HXB_04_Brain")]
HXB26 <- as.matrix(apply(HXB26_dat, 1, mean))

HXB27_dat = heritExprs[, c("X27HXB_01_Brain", "X27HXB_02_Brain", "X27HXB_03_Brain", "X27HXB_04_Brain")]
HXB27 <- as.matrix(apply(HXB27_dat, 1, mean))

HXB29_dat = heritExprs[, c("X29HXB_01_Brain", "X29HXB_02_Brain", "X29HXB_03_Brain", "X29HXB_04_Brain")]
HXB29 <- as.matrix(apply(HXB29_dat, 1, mean))

HXB3_dat = heritExprs[, c("X3HXB_01_Brain", "X3HXB_02_Brain", "X3HXB_03_Brain", "X3HXB_04_Brain")]
HXB3 <- as.matrix(apply(HXB3_dat, 1, mean))

HXB31_dat = heritExprs[, c("X31HXB_01_Brain", "X31HXB_02_Brain", "X31HXB_03_Brain", "X31HXB_04_Brain")]
HXB31 <- as.matrix(apply(HXB31_dat, 1, mean))

HXB4_dat = heritExprs[, c("X4HXB_01_Brain", "X4HXB_02_Brain", "X4HXB_03_Brain", "X4HXB_04_Brain")]
HXB4 <- as.matrix(apply(HXB4_dat, 1, mean))

HXB7_dat = heritExprs[, c("X7HXB_01_Brain", "X7HXB_02_Brain", "X7HXB_03_Brain", "X7HXB_04_Brain")]
HXB7 <- as.matrix(apply(HXB7_dat, 1, mean))

PD8604_dat = heritExprs[, c("PD_8604_03_Brain", "PD_8604_04_Brain")]
PD8604 <- as.matrix(apply(PD8604_dat, 1, mean))

PD8664_dat = heritExprs[, c("PD_8664_01_Brain", "PD_8664_02_Brain", "PD_8664_03_Brain", "PD_8664_04_Brain")]
PD8664 <- as.matrix(apply(PD8664_dat, 1, mean))

SHRH_dat = heritExprs[, c("SHR.H_01_Brain", "SHR.H_02_Brain", "SHR.H_03_Brain", "SHR.H_04_Brain")]
SHRH <- as.matrix(apply(SHRH_dat, 1, mean))

SHRlj_dat = heritExprs[, c("SHR.lj_01_Brain", "SHR.lj_02_Brain", "SHR.lj_03_Brain", "SHR.lj_04_Brain")]
SHRlj <- as.matrix(apply(SHRlj_dat, 1, mean))

SHRlx_dat = heritExprs[, c("SHR.Lx_01_Brain", "SHR.Lx_02_Brain", "SHR.Lx_03_Brain", "SHR.Lx_04_Brain")]
SHRlx <- as.matrix(apply(SHRlx_dat, 1, mean))

Wky_dat = heritExprs[, c("WKY.Lj_01_Brain", "WKY.Lj_02_Brain", "WKY.Lj_03_Brain", "WKY.Lj_04_Brain")]
Wky <- as.matrix(apply(Wky_dat, 1, mean))

################################################
means <- cbind(BNLX, BXH08, BXH10, BXH11, BXH13, BXH6, HXB1, HXB10, HXB13, HXB15, HXB17, HXB18, HXB2, HXB23,
		HXB25, HXB26, HXB27, HXB29, HXB3, HXB31, HXB4, HXB7, PD8604, PD8664, SHRH, SHRlj, SHRlx, Wky)

colnames(means) <- c("BNLX", "BXH08", "BXH10", "BXH11", "BXH13", "BXH6", "HXB1", "HXB10", "HXB13", "HXB15", 
			"HXB17", "HXB18", "HXB2", "HXB23", "HXB25", "HXB26", "HXB27", "HXB29", "HXB3", "HXB31", 
			"HXB4", "HXB7", "PD8604", "PD8664", "SHRH", "SHRlj", "SHRlx", "Wky")

means[1:5,]

save(means, file="datExpr_means_precollapse.RData")

##############clean means################
##Now Check for Outliers###;
sampleTree = flashClust(as.dist(1-cor(means)), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Strain clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)

###No BAD Strains 

############################################################################################
##  collapse probe sets that are associated with the same gene symbol and are correlated  ##
############################################################################################
library(WGCNA)
load("probeToCluster.ratExon.Rdata")

collapsingInfo <- probeToCluster[match(rownames(means),probeToCluster[,1]),]
collapsingInfo <- collapsingInfo[match(collapsingInfo[,1],rownames(means)),]
table(collapsingInfo[,1]==rownames(means))


geneSymbols = collapsingInfo[,2]
names(geneSymbols) = collapsingInfo[,1]
collapsedExprs <- collapseRows(datET=as.matrix(as.data.frame(means)),rowID=rownames(means),rowGroup=geneSymbols,thresholdCombine=0.5)

dim(collapsedExprs[[1]])
## 49,686 probe sets

save(collapsedExprs,file="collapsedExprs.Rdata")

