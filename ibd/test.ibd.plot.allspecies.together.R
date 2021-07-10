library(geosphere)
library(phangorn)
library(MASS)
library(ade4)
library(scales)
library(viridis)


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


#plot
pdf(file = "~/Dropbox/Chalcophaps/spring.2020/ibd.byspecies.pdf",width =6.5,height = 6)

layout.matrix <- matrix(c(1,2,1,3,1,4), nrow = 2, ncol = 3)

layout(mat = layout.matrix,
       heights = c(1,.75), # Heights of the two rows
       widths = c(1,1,1)) # Widths of the three columns

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

#plot residuals
steph.lm<-lm(stephani.gendist~stephani.Dgeo)
steph.res = resid(steph.lm)
plot(stephani.Dgeo, steph.res, ylab="Residuals", xlab="Geographic distance", 
     main=expression(italic('stephani')), col = alpha("sienna2", 0.6), pch=16) 
  abline(0, 0)
  
longi.lm<-lm(longi.gendist~longi.Dgeo)
longi.res = resid(longi.lm)
plot(longi.Dgeo, longi.res, ylab="Residuals", xlab="Geographic distance", 
     main=expression(italic('longirostris')), col = alpha("skyblue1", 0.6), pch=16) 
abline(0, 0)   

indi.lm<-lm(indica.gendist~indica.Dgeo)
indi.res = resid(indi.lm)
plot(indica.Dgeo, indi.res, ylab="Residuals", xlab="Geographic distance", 
     main=expression(italic('indica')), col = alpha("palegreen2", 0.6), pch=16) 
abline(0, 0)   

dev.off()









#alternatively, if you want to plot each cluster separately
#separately plot relationship within each pop cluster for stephani (start with solomons)
#test IBD stephani start with dgeo
sols.spec<-stephani.spec[stephani.spec$Country == "Solomon Islands",]
sols.xy.coords.only<- subset(sols.spec, select=c("Long","Lat"))
sols.Dgeo <- as.dist(distm(sols.xy.coords.only, fun=distGeo))
sols.Dgeo<-sols.Dgeo/1000
sols.Dgeo

#now genetic matrix
#create subset for each species matching the order of the geographic dataframe
sols.concat <- subset(all_concat, subset = c(79:84))
#create distance matrix based on all genes concatenated
sols.gendist<-dist.hamming(sols.concat, ratio = TRUE)
IBD <- mantel.randtest(sols.Dgeo,sols.gendist)
IBD
plot(sols.Dgeo,sols.gendist, pch=20,cex=1.5, ylab = "Genetic Distance", xlab = "Geographic Distance")
abline(lm(sols.gendist~sols.Dgeo), lty = 2)


#separately plot relationship within each pop cluster for stephani (now PNG without New Britain)
#test IBD stephani start with dgeo
png.spec<-stephani.spec[stephani.spec$Country == "Papua New Guinea",]
png.spec<-png.spec[-4,]
png.xy.coords.only<- subset(png.spec, select=c("Long","Lat"))
png.Dgeo <- as.dist(distm(png.xy.coords.only, fun=distGeo))
png.Dgeo<-png.Dgeo/1000
png.Dgeo

#now genetic matrix
#create subset for each species matching the order of the geographic dataframe
png.concat <- subset(all_concat, subset = c(70,71,78,73,74,75,76,77))
#create distance matrix based on all genes concatenated
png.gendist<-dist.hamming(png.concat, ratio = TRUE)
IBD <- mantel.randtest(png.Dgeo,png.gendist)
IBD
plot(png.Dgeo,png.gendist, pch=20,cex=1.5, ylab = "Genetic Distance", xlab = "Geographic Distance")
abline(lm(png.gendist~png.Dgeo), lty = 2)

#plot both together
plot(png.Dgeo, png.gendist, pch=20,cex=1.5, ylab = "Uncorrected 'p' Distance",
     xlab = "Geographic Distance (km)", ylim=c(0,.0012), col = "#D54941FF", frame.plot = FALSE,
     main=expression( italic(stephani)))
points(sols.Dgeo, sols.gendist, pch=20,cex=1.5, col="#FBBE22FF")
abline(lm(png.gendist~png.Dgeo), lty = 1, lwd=4, col="#D54941FF")
abline(lm(sols.gendist~sols.Dgeo), lty =1, lwd=4, col="#FBBE22FF")
legend("topleft",
       legend = c("PNG", "Solomon Islands"),
       #text.font = 3,
       col = c("#D54941FF","#FBBE22FF"), 
       pch = c(20), 
       bty = "n", 
       pt.cex = 3, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F)


#repeat the same for longirostris
#separately plot relationship within each pop cluster for longi (start with solomons)
#test IBD stephani start with dgeo
sols.spec<-stephani.spec[stephani.spec$Country == "Solomon Islands",]
sols.xy.coords.only<- subset(sols.spec, select=c("Long","Lat"))
sols.Dgeo <- as.dist(distm(sols.xy.coords.only, fun=distGeo))
sols.Dgeo<-sols.Dgeo/1000
sols.Dgeo

#now genetic matrix
#create subset for each species matching the order of the geographic dataframe
sols.concat <- subset(all_concat, subset = c(79:84))
#create distance matrix based on all genes concatenated
sols.gendist<-dist.hamming(sols.concat, ratio = TRUE)
IBD <- mantel.randtest(sols.Dgeo,sols.gendist)
IBD
plot(sols.Dgeo,sols.gendist, pch=20,cex=1.5, ylab = "Genetic Distance", xlab = "Geographic Distance")
abline(lm(sols.gendist~sols.Dgeo), lty = 2)




