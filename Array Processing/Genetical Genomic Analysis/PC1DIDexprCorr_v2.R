rm(list=ls())
setwd("/data/kiemele/LXS/LXSpaper2014/DID/version2")

#############
# Load Data #
#############
#expression data
load("datExprStrainMeans.RData")
exprs = t(means)

#PC1 DID Data 
traits = read.csv(file="/data/kiemele/LXS/LXSpaper2014/DID/logEPCA_strainSexMeans.csv", header=TRUE)
traits

males = traits[which(traits$Sex=="M"),]
males = males[,c("Strain", "Prin1_Mean")]
colnames(males) = c("Strain", "Male_PC1DID")

females = traits[which(traits$Sex=="F"),]
females = females[,c("Strain", "Prin1_Mean")]
colnames(females) = c("Strain", "Female_PC1DID")

traits = merge(males, females, by="Strain", all=TRUE)
rownames(traits) = traits$Strain
traits = traits[,c(2,3)]

############################
# Correlate Traits & Exprs #
############################
dim(exprs)
dim(traits)
colnames(traits)

test = merge(traits, exprs, by.x = 0, by.y=0,all.y=TRUE)
dim(test)
test[1:5, 1:5]

traitNames <- colnames(traits)

phenoExprCorrs <- data.frame(trans=colnames(exprs))

for(i in traitNames){
  phenoExprCorrs <- cbind(phenoExprCorrs,apply(test[,as.character(phenoExprCorrs[,"trans"])],2,function(a) cor.test(test[,i],as.numeric(a))$estimate))
  phenoExprCorrs <- cbind(phenoExprCorrs,apply(test[,as.character(phenoExprCorrs[,"trans"])],2,function(a) cor.test(test[,i],as.numeric(a))$p.value))
}

colnames(phenoExprCorrs)<- c("trans",paste(rep(traitNames,each=2),rep(c(".coef",".pvalue"),length(traitNames)),sep=""))
dim(phenoExprCorrs)
head(phenoExprCorrs)

#Calculate FDR Values ###NOTE: need to change column name from pheno.pvalue to whatever your specific column name is for phenotype.pvalue;
phenoExprCorrs$maleFDR = p.adjust(phenoExprCorrs[,"Male_PC1DID.pvalue"], method="BH") 
phenoExprCorrs$femaleFDR = p.adjust(phenoExprCorrs[,"Female_PC1DID.pvalue"], method="BH") 

dim(phenoExprCorrs)
head(phenoExprCorrs)

#Write to CSV file
write.csv(phenoExprCorrs, file="PC1DID_Exprs.correlations.csv")

#########################
# Get Significant Exprs #
#########################
rm(list=ls())
setwd("/data/kiemele/LXS/LXSpaper2014/DID/version2")
phenoExprCorrs = read.csv(file="PC1DID_Exprs.correlations.csv",header=TRUE,row.names=1)
summary(phenoExprCorrs$maleFDR)
summary(phenoExprCorrs$femaleFDR)
summary(phenoExprCorrs$Male_PC1DID.pvalue)
summary(phenoExprCorrs$Female_PC1DID.pvalue)



Msig.results = phenoExprCorrs[-which(phenoExprCorrs[,"Male_PC1DID.pvalue"]>0.001),]
dim(Msig.results)

Fsig.results = phenoExprCorrs[-which(phenoExprCorrs[,"Female_PC1DID.pvalue"]>0.001),]
dim(Fsig.results)

BothSig.results = rbind(Msig.results, Fsig.results)
write.csv(BothSig.results, file="PC1DIDexpr.sigCorrs.csv", row.names=FALSE, quote=FALSE)

########################################
# Get Expr of Sig Transcripts for eQTL #
########################################

rm(list=ls())
setwd("/data/kiemele/LXS/LXSpaper2014/DID/version2")
BothSig.results = read.csv(file="PC1DIDexpr.sigCorrs.csv",header=TRUE,row.names=1)

load("datExprStrainMeans.RData")

want = as.matrix(rownames(BothSig.results))

sig = c()
for(i in want){
	sig = c(sig, which(rownames(means)==i))
}

length(sig)
length(want)

Expr.sigCorrTranscripts = t(means[sig,])
write.csv(Expr.sigCorrTranscripts, file="Expr.sigCorrTranscripts.csv")




