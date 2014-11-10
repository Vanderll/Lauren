####################
# Set Up Workspace #
####################

rm(list=ls())
setwd("C:/Users/kiemelel/Documents/BorisRandom/20140903_cAMPassay")

library(drc)
library(plotrix)

###########################
# Load All Treatment Data #
###########################

cAMP = read.csv(file="20140903_cAMPassay.csv", header=TRUE)
dim(cAMP)
head(cAMP)
table(cAMP$Treatment)

# Get B0 and NSB averages

NSB = mean(cAMP[which(cAMP$Treatment == "NSB"), "BlankAdj"])
B0 = mean(cAMP[which(cAMP$Treatment == "B0"), "BlankAdj"])

# Get standard curve info;


sc2250 = mean(cAMP[which(cAMP$std_dose == 2250), "BlankAdj"])
sc6750 = mean(cAMP[which(cAMP$std_dose == 6750), "BlankAdj"])

sc750 = mean(cAMP[which(cAMP$std_dose == 750), "BlankAdj"])
sc250 = mean(cAMP[which(cAMP$std_dose == 250), "BlankAdj"])
sc83.3 = mean(cAMP[which(cAMP$std_dose == 83.3), "BlankAdj"])
sc27.8 = mean(cAMP[which(cAMP$std_dose == 27.8), "BlankAdj"])
sc9.3 = mean(cAMP[which(cAMP$std_dose == 9.3), "BlankAdj"])
sc3.1 = mean(cAMP[which(cAMP$std_dose == 3.1), "BlankAdj"])
sc1 = mean(cAMP[which(cAMP$std_dose == 1), "BlankAdj"])
sc.3 = mean(cAMP[which(cAMP$std_dose == 0.3), "BlankAdj"])

standardCurve = as.matrix(rbind(sc6750, sc2250, sc750, sc250, sc83.3, sc27.8, sc9.3, sc3.1, sc1, sc.3))
dose = as.matrix(c(6750, 2250, 750, 250, 83.3, 27.8, 9.3, 3.1, 1, 0.3))
standardCurve = cbind(standardCurve, dose)
colnames(standardCurve) = c("A405means", "dose")

correctA = standardCurve[,"A405means"] - NSB
ratio = correctA/(B0 - NSB)
percentRatio = ratio*100

standardCurve = as.data.frame(cbind(standardCurve, correctA, ratio, percentRatio))


### Make Dose Response Curve Model
curve = rep(1, nrow(standardCurve))

drc.model.sc = drm(percentRatio~dose, curve, data=standardCurve, fct = LL.4())
summary(drc.model.sc)
save(drc.model.sc, file="20140903_cAMP_stdCurve.Rdata")

plot(drc.model.sc, main="cAMP assay (20140903): standard curve 4 parameter drc", ylab = "percent B/B0 ratio", 
	xlab = "Non-Acetylated Cyclic AMP (pmol/ml)", cex.lab = 1, cex.axis = 1,pch=19)

	
abline(h=7,lty=2, col="red")
abline(h=10,lty=2, col="blue")
abline(h=15,lty=2, col="green")


### Get Sample B/B0% And Use Standard Curve to Calculate [Non-Acetylated cAMP]

#get % ratio
cAMP$correctedA = cAMP[,"BlankAdj"] - NSB
cAMP$ratio = cAMP$correctedA/(B0 - NSB)
cAMP$percentRatio = cAMP$ratio*100

head(cAMP)
# remove samples that you don't want/need the [cAMP]
IBMXsamples = cAMP[-which(cAMP$Treatment == "NSB"),]
IBMXsamples = IBMXsamples[-which(IBMXsamples$Treatment == "B0"),]
IBMXsamples = IBMXsamples[-which(IBMXsamples$Treatment == "Standard"),]
IBMXsamples = IBMXsamples[-which(IBMXsamples$Treatment == "blank"),]
head(IBMXsamples)

table(IBMXsamples$Treatment)

getConcentration = IBMXsamples$percentRatio
IBMXsamples[,c("Treatment", "percentRatio", "BlankAdj")]

results =  ED(drc.model.sc, getConcentration, type="absolute")

IBMXresults = cbind(IBMXsamples, results)
IBMXresults[,c("Treatment", "BlankAdj", "percentRatio", "Estimate")]
head(IBMXresults)

write.csv(IBMXresults, file="cAMP20140903.results.csv", row.names=FALSE)








