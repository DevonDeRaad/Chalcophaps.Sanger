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

library(ggplot2)
library(viridis)
library(gridExtra)

#import csv with locality info
chalc.mapping<- read.csv("~/Downloads/included.samples.csv")
table(chalc.mapping$Subspecies)
table(chalc.mapping$Lat)

#split df by species
spec.dfs<-split(chalc.mapping, chalc.mapping$Species)

#use dplyr to make an indica table with unique lat longs only and the count
indica.sampling<-spec.dfs$indica %>% group_by(Lat, Long) %>% summarize(count=n())
#add species column
indica.sampling$species<-rep("indica", times=nrow(indica.sampling))

#use dplyr to make an indica table with unique lat longs only and the count
longi.sampling<-spec.dfs$longirostris %>% group_by(Lat, Long) %>% summarize(count=n())
#add species column
longi.sampling$species<-rep("longirostris", times=nrow(longi.sampling))

#use dplyr to make an indica table with unique lat longs only and the count
stephani.sampling<-spec.dfs$stephani %>% group_by(Lat, Long) %>% summarize(count=n())
#add species column
stephani.sampling$species<-rep("stephani", times=nrow(stephani.sampling))

#combine in tidy format
all.sampling<-as.data.frame(rbind(indica.sampling,longi.sampling,stephani.sampling))
all.sampling$species<-as.factor(all.sampling$species)

#make map 1: stephani
pac<-map_data("world")
#make full map
map1<-ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(140, 165), ylim = c(-14, 5)) + 
  geom_point(data = all.sampling[all.sampling$species == "stephani",], aes(x = Long, y = Lat, size=count,
             color=c("png","sol","png","sol","sol","sol","png","bis","png","png")), alpha =.8) +
  theme_classic()+
  scale_color_viridis(discrete = TRUE, option = "inferno", begin=.3, end=.85)+
  scale_size_continuous(range = c(3,5), breaks = c(1,2,3))+
  theme(legend.position = "none")+
  ggtitle("stephani")+
  theme(plot.title = element_text(face="italic"),axis.title.x=element_blank(),axis.title.y=element_blank())


#make map 2: longirostris
map2<-ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(120, 170), ylim = c(-35, 3)) + 
  geom_point(data = all.sampling[all.sampling$species == "longirostris",], aes(x = Long, y = Lat, size=count,
             color=c(rep("aus", times=6),"van","aus","aus","png","van",rep("png", times=3),"tim",rep("png", times=10))), alpha =.8) +
  theme_classic()+
  scale_color_viridis(discrete = TRUE)+
  scale_size_continuous(range = c(3,5), breaks = c(1,2,3))+
  theme(legend.position = "none")+
  ggtitle("longirostris")+
  theme(plot.title = element_text(face="italic"),axis.title.x=element_blank(),axis.title.y=element_blank())


#make map 3: indica
map3<-ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(95, 145), ylim = c(-13, 25)) + 
  geom_point(data = all.sampling[all.sampling$species == "indica",], aes(x = Long, y = Lat, size=count,
             color=c(rep("xmas", times=3),rep("malaysia", times=2),rep("phil", times=10),"viet",rep("phil", times=7),rep("china", times=2)), alpha =.9)) +
  theme_classic()+
  scale_color_viridis(discrete = TRUE, option = "plasma")+
  scale_size_continuous(range = c(3,5), breaks = c(1,2,3))+
  theme(legend.position = "none")+
  ggtitle("indica")+
  theme(plot.title = element_text(face="italic"),axis.title.x=element_blank(),axis.title.y=element_blank())

grid.arrange(map3,map2,map1, nrow=3)
#g <- arrangeGrob(map3,map2,map1, nrow=3) #generates g
#
ggsave("~/Dropbox/Chalcophaps/species.maps.pdf", g, width = 168, height = 239*.667, units = "mm")





#Didn't end up using these haplotype networks because they render so poorly.
#made them with popart instead for the paper

#make haplotype networks
all_renamed_nd2 <- read.phyDat(file.path("~/Dropbox/Chalcophaps/all_renamed_ND2.fasta"), format = "fasta")
all_renamed_ordered_nd2 <- subset(all_renamed_nd2, subset = c(2:14,16:36,67,62,65:66,40:61,63,37:39,15,64,1,68:70,71:73,80:85,74:79))
longi.fas<-subset(all_renamed_ordered_nd2, subset = c(1:34,36:38))
indi.fas<-subset(all_renamed_ordered_nd2, subset = c(35,39:70))
steph.fas<-subset(all_renamed_ordered_nd2, subset = c(71:85))

#plot longi
nd2.bin<-as.DNAbin(longi.fas)
nd2haps<-haplotype(nd2.bin)
nd2.net<- haploNet(nd2haps)
#build a dataframe that tells you which individuals correspond to which haplotypes
nd2.hapcalls<-stack(setNames(attr(nd2haps, "index"), rownames(nd2haps)))
individs<-attr(longi.fas, "names")
individs<-individs[nd2.hapcalls$values]
cbind(individs,nd2.hapcalls$ind)
#start calculating haplotype frequencies
specs<-c(rep("aus", times=34),rep("van", times=3))
nd2.ind.hap2<-table(hap=nd2.hapcalls$ind,pop=specs)
nd2.ind.hap2
par(bg=NA)
dev.off()
#plot colored by species
pdf(file = "~/Downloads/longi.hapnetwork.pdf",width =3.11 ,height = 3.11)
plot(nd2.net, size=attr(nd2.net, "freq"), labels = F, scale.ratio = 1.2, cex = 1, pie=nd2.ind.hap2, show.mutation = 2, threshold = 0,
     bg=c("#a6cee3","#1f78b4"))
dev.off()

#plot indi
nd2.bin<-as.DNAbin(indi.fas)
nd2haps<-haplotype(nd2.bin)
nd2.net<- haploNet(nd2haps)
#build a dataframe that tells you which individuals correspond to which haplotypes
nd2.hapcalls<-stack(setNames(attr(nd2haps, "index"), rownames(nd2haps)))
individs<-attr(indi.fas, "names")
individs<-individs[nd2.hapcalls$values]
cbind(individs,nd2.hapcalls$ind)
#start calculating haplotype frequencies
specs<-c("long",rep("indi", times=19),"xmas",rep("indi", times=9),rep("xmas", times=3))
nd2.ind.hap2<-table(hap=nd2.hapcalls$ind,pop=specs)
nd2.ind.hap2
par(bg=NA)
dev.off()
#plot colored by species
pdf(file = "~/Downloads/indi.hapnetwork.pdf",width =3.11 ,height = 3.11)
plot(nd2.net, size=attr(nd2.net, "freq"), labels = F,
                   scale.ratio = 1.2, cex = 1, pie=nd2.ind.hap2, show.mutation = 2, threshold = 0,
                    bg=c("forestgreen","blue","#a6d854"))
dev.off()

#plot stephani
nd2.bin<-as.DNAbin(steph.fas)
nd2haps<-haplotype(nd2.bin)
nd2.net<- haploNet(nd2haps)
#build a dataframe that tells you which individuals correspond to which haplotypes
nd2.hapcalls<-stack(setNames(attr(nd2haps, "index"), rownames(nd2haps)))
individs<-attr(steph.fas, "names")
individs<-individs[nd2.hapcalls$values]
cbind(individs,nd2.hapcalls$ind)
#start calculating haplotype frequencies
specs<-c(rep("png", times=6),"bis","png","png",rep("sol", times=6))
nd2.ind.hap2<-table(hap=nd2.hapcalls$ind,pop=specs)
nd2.ind.hap2
par(bg=NA)
dev.off()
#plot colored by species
pdf(file = "~/Downloads/steph.hapnetwork.pdf",width =3.11,height = 3.11)
plot(nd2.net, size=attr(nd2.net, "freq"), labels = F, scale.ratio = .5, pie=nd2.ind.hap2, show.mutation = 2, threshold = 0,
                   bg=c("darkorange","#fec44f","#d95f0e"), lwd =.3)
dev.off()




