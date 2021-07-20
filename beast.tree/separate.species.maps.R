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

library(ggplot2)
library(ggtree)


#Chalcophaps analysis

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

#make Figure 1#
pac<-map_data("world")
#"palegreen","palegreen4","powderblue","skyblue1","royalblue3","navyblue","lightsalmon","sienna3")
#make full map
map<-ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(92, 166), ylim = c(-40, 28)) + 
  geom_point(data = all.sampling, aes(x = Long, y = Lat, col=species, size=count), alpha =.8, show.legend=TRUE) +
  theme_classic()+
  scale_color_manual(values=c("palegreen2","skyblue1","sienna2"))+
  scale_size_continuous(range = c(3,5), breaks = c(1,2,3))+
  guides(colour = guide_legend(override.aes = list(size = 4), order=1, label.theme = element_text(face = "italic")),
         size = guide_legend(nrow = 1, order = 2))+
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0.01, 0.01))+
  labs(x = "longitude", y = "latitude")

map
#map complete#

#shadow map
ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="black", col="black", cex=0)+
  coord_sf(xlim = c(75, 166), ylim = c(-40, 30))+ 
  theme_void()

ggsave(filename = "~/Downloads/shadow.map.pdf", width = 6, height=5)


#indica map
ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(75, 136), ylim = c(-10, 30)) + 
  geom_point(data = all.sampling[all.sampling$species == "indica",], aes(x = Long, y = Lat, col=species, size=count), alpha =.8, show.legend=TRUE) +
  theme_void()+
  scale_color_manual(values=c("palegreen2"))+
  scale_size_continuous(range = c(3,5), breaks = c(1,2,3))+
  guides(colour = guide_legend(override.aes = list(size = 4), order=1, label.theme = element_text(face = "italic")),
         size = guide_legend(nrow = 1, order = 2))+
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0.01, 0.01), legend.title = element_blank())+
  labs(x = "longitude", y = "latitude")

ggsave(filename = "~/Downloads/indica.map.pdf", width = 6, height=5)

#longirostris map
ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(120, 166), ylim = c(-40, -3)) + 
  geom_point(data = all.sampling[all.sampling$species == "longirostris",], aes(x = Long, y = Lat, col=species, size=count), alpha =.8, show.legend=TRUE) +
  theme_void()+
  scale_color_manual(values=c("skyblue1"))+
  scale_size_continuous(range = c(3,5), breaks = c(1,2,3))+
  guides(colour = guide_legend(override.aes = list(size = 4), order=1, label.theme = element_text(face = "italic")),
         size = guide_legend(nrow = 1, order = 2))+
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0.01, 0.01), legend.title = element_blank())+
  labs(x = "longitude", y = "latitude")

ggsave(filename = "~/Downloads/longirostris.map.pdf", width = 6, height=5)


#stephani map
ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(119, 163), ylim = c(-12, 2)) + 
  geom_point(data = all.sampling[all.sampling$species == "stephani",], aes(x = Long, y = Lat, col=species, size=count), alpha =.8, show.legend=TRUE) +
  theme_void()+
  scale_color_manual(values=c("sienna2"))+
  scale_size_continuous(range = c(3,5), breaks = c(1,2,3))+
  guides(colour = guide_legend(override.aes = list(size = 4), order=1, label.theme = element_text(face = "italic")),
         size = guide_legend(nrow = 1, order = 2))+
  theme(legend.position = c(0.8, 0.55), legend.justification = c(0.01, 0.01), legend.title = element_blank())+
  labs(x = "longitude", y = "latitude")

ggsave(filename = "~/Downloads/stephani.map.pdf", width = 6, height=5)
