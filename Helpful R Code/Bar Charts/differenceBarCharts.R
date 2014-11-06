rm(list=ls())
setwd("C:/Users/kiemelel/Documents/LXS/LXSpaper2014/DID/results")
options(stringsAsFactors=FALSE)


### Load phenotype Data into R;


#PCA of logE DID values;
info = read.csv(file="interestingMarkers_v2.csv")
head(info)
table(info$marker)

#Original Values;
orig = read.csv(file="DIDmeans_4peaks.csv")


#################################
# BarChart Showing Each SNP*Sex #
#################################

infoList = list()
markers = as.matrix(unique(info$rsNumber))
for(i in markers){
	infoList[[i]] = as.matrix(info[which(info$rsNumber==i), c("GenoL", "GenoS")])
	rownames(infoList[[i]])= c("F", "M")
}

stderrList = list()
markers = as.matrix(unique(info$rsNumber))
for(i in markers){
	stderrList[[i]] = as.matrix(info[which(info$rsNumber==i), c("GenoLStdErr", "GenoSStdErr")])
	rownames(stderrList[[i]])= c("F", "M")
}


####Get it to work in a loop;
for(i in 1:length(infoList)){
	fileName = paste(names(infoList[i]), ".SexSNPbarChart.v2.jpg", sep="")
	jpeg(file=fileName)

   BP <- barplot(infoList[[i]], 
	main=unique(paste(names(infoList[i]), " (Chr", info[which(info$rsNumber==names(infoList[i])),"chr"], ": ", round(info[which(info$rsNumber==names(infoList[i])),"Mb"], digits=1), " Mb)", sep="")),
	beside=TRUE,
	xlab="Genotype",
	ylab="1st Principal Component Score DID phenotypes",
	names.arg = c("L", "S"),
	args.legend=TRUE,
	col=c("red", "blue"),
	ylim = c(min(infoList[[i]]-stderrList[[i]]), max(infoList[[i]]+stderrList[[i]]))
   )

   arrows(x0=BP,y0=infoList[[i]],x1=BP,y1=infoList[[i]]+stderrList[[i]],angle=90, length=.15)
   arrows(x0=BP,y0=infoList[[i]],x1=BP,y1=infoList[[i]]-stderrList[[i]],angle=90, length=.15)
   abline(0,0)
   legend("topleft", c("Females","Males"),cex=0.85, bty="n", fill = c("red", "blue"))

dev.off()
}




#################################################################
# BarChart Showing Each SNP*Sex for original DID values stacked #
#################################################################

#Original Values;
orig = as.matrix(read.csv(file="DIDmeans_4peaks_difference.csv", row.names=1))

plot.new()
BP <- barplot(orig, main="Genotype Differences (L-S) \nin Alcohol Consumption (DID) for Days 2, 3, and 4",
	xlab="Alcohol Consumption Phenotype", 
   	ylab="Difference in AlcoholConsumption (g/kg)", 
	las=3, 
	ylim=c(-2,2),
	names.arg = rep("", ncol(orig)),
	cex.names=0.75,
	beside=TRUE,
	args.legend = TRUE,
	col=c("red", "blue"))

	legend("topright", c("female","male"),cex=0.85, bty="n", fill = c("red", "blue"))
	abline(h = 0, col = "black")
			
 	text(10.5, -1.5, "Chr 1 Peak               Chr 4 Peak             Ch 7 Peak               Chr 19 Peak", cex = .8)

text(4, 7, expression(bar(x) == sum(frac(x[i], n), i==1, n)))