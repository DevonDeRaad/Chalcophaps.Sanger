library(phangorn)
library(reshape2)
library(seqRFLP)
library(phylotools)
library(pegas)
library(geoR)
library(ape)
library(ggplot2)
library(dplyr)
library(ggtree)
library(tidytree)
library(gridExtra)
library(geosphere)

library(geosphere)
library(phangorn)
library(MASS)
library(ade4)


### Isolation by distance
chalc.mapping<- read.csv("~/Downloads/included.samples.csv")
#split by species
indica.spec<-chalc.mapping[chalc.mapping$Species == "indica",]
#longirostris
longi.spec<-chalc.mapping[chalc.mapping$Species == "longirostris",]
#stephani
stephani.spec<-chalc.mapping[chalc.mapping$Species == "stephani",]

indica.xy.coords.only<- subset(indica.spec, select=c("Long","Lat"))
rownames(indica.xy.coords.only)<-indica.spec$Tissue
longi.xy.coords.only<- subset(longi.spec, select=c("Long","Lat"))
rownames(longi.xy.coords.only)<-longi.spec$Tissue
stephani.xy.coords.only<- subset(stephani.spec, select=c("Long","Lat"))
rownames(stephani.xy.coords.only)<-stephani.spec$Tissue

#indica make dist matrix convert to km distance
indica.Dgeo <- as.dist(distm(indica.xy.coords.only, fun=distGeo))
indica.Dgeo<-indica.Dgeo/1000
indica.Dgeo
#
longi.Dgeo <- as.dist(distm(longi.xy.coords.only, fun=distGeo))
longi.Dgeo<-longi.Dgeo/1000
longi.Dgeo
#steph
stephani.Dgeo <- as.dist(distm(stephani.xy.coords.only, fun=distGeo))
stephani.Dgeo<-stephani.Dgeo/1000
stephani.Dgeo


#we now have Dgeo which is a distance matrix of euclidean distance btwn points (class dist)
#create a genetic distance matrix (class matrix) with same individual order as Dgeo
all_concat <- read.phyDat(file.path("~/Dropbox/Chalcophaps/allgenes_concatenated_nooutgroup_namedbysubspecies.fasta"), format = "fasta")
#create subset for each species matching the order of the geographic dataframe
stephani.concat <- subset(all_concat, subset = c(79:84,70,71,78,72:77))
indica.concat <- subset(all_concat, subset = c(14,36:60,62,63,1,67:69))
longi.concat <- subset(all_concat, subset = c(8,2:7,9:13,15:35,61,64:66))

#create distance matrix based on all genes concatenated
indica.gendist<-phangorn::dist.hamming(indica.concat, ratio = TRUE) #this computes raw, uncorrected 'p' distance
#create distance matrix based on all genes concatenated
longi.gendist<-dist.hamming(longi.concat, ratio = TRUE)
#create distance matrix based on all genes concatenated
stephani.gendist<-dist.hamming(stephani.concat, ratio = TRUE)


#test IBD indica
IBD <- mantel.randtest(indica.Dgeo,indica.gendist)
IBD
plot(indica.Dgeo,indica.gendist, pch=20,cex=1.5, ylab = "Genetic Distance", xlab = "Geographic Distance")
abline(lm(indica.gendist~indica.Dgeo), lty = 2)
#plot and check for denser areas in the plot indicating sub-groups
dens <- kde2d(indica.Dgeo,indica.gendist, n=300, lims=c(-5, 35, 0, 0.08))
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(indica.Dgeo, indica.gendist, pch=20,cex=1.5, ylab = "Nei's Genetic Distance", xlab = "Geographic Distance")
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(indica.gendist~indica.Dgeo), lty =2)
title("indica")


#test IBD longirostris
IBD <- mantel.randtest(longi.Dgeo,longi.gendist)
IBD
plot(longi.Dgeo,longi.gendist, pch=20,cex=1.5, ylab = "Genetic Distance", xlab = "Geographic Distance")
abline(lm(longi.gendist~longi.Dgeo), lty = 2)
#plot and check for denser areas in the plot indicating sub-groups
dens <- kde2d(longi.Dgeo,longi.gendist, n=300, lims=c(-5, 45, 0, 0.08))
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(longi.Dgeo, longi.gendist, pch=20,cex=1.5, ylab = "Nei's Genetic Distance", xlab = "Geographic Distance")
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(longi.gendist~longi.Dgeo), lty =2)
title("Correlation of Genetic and Geographical distances")

#test IBD stephani
IBD <- mantel.randtest(stephani.Dgeo,stephani.gendist)
IBD
plot(stephani.Dgeo,stephani.gendist, pch=20,cex=1.5, ylab = "Genetic Distance", xlab = "Geographic Distance")
abline(lm(stephani.gendist~stephani.Dgeo), lty = 2)
#plot and check for denser areas in the plot indicating sub-groups
dens <- kde2d(stephani.Dgeo,stephani.gendist, n=300, lims=c(-5, 45, 0, 0.08))
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(stephani.Dgeo, stephani.gendist, pch=20,cex=1.5, ylab = "Nei's Genetic Distance", xlab = "Geographic Distance")
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(stephani.gendist~stephani.Dgeo), lty =2)
title("Correlation of Genetic and Geographical distances")


pdf(file = "~/Dropbox/Chalcophaps/spring.2020/ibd.byspecies.pdf",width =5,height = 5)
#plot them all together and color points by species
plot(indica.Dgeo, indica.gendist, pch=20,cex=1.5, ylab = "Uncorrected 'p' Distance",
     xlab = "Geographic Distance (km)", ylim=c(0,.03), col = alpha("palegreen2", 0.6), frame.plot = FALSE)
points(longi.Dgeo, longi.gendist, pch=20,cex=1.5, col=alpha("skyblue1",.6))
points(stephani.Dgeo, stephani.gendist, pch=20,cex=1.5, col=alpha("sienna2",.6))
abline(lm(stephani.gendist~stephani.Dgeo), lty = 1, lwd=4, col="sienna2")
abline(lm(longi.gendist~longi.Dgeo), lty =1, lwd=4, col="skyblue1")
abline(lm(indica.gendist~indica.Dgeo), lty =1, lwd=4, col="palegreen2")
legend("topright",
       legend = c("indica", "longirostris","stephani"),
       text.font = 3,
       col = c("palegreen2","skyblue1","sienna2"), 
       pch = c(20), 
       bty = "n", 
       pt.cex = 3, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F)

dev.off()

