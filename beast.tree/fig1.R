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

#map complete#

#read in tree
beast_tree <- read.beast("~/Dropbox/Chalcophaps/beast.beauti/sy.sum.tree")

#tree no color coding
ggtree(beast_tree)+ geom_tiplab(cex=2)+
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior) == 1), cex=3, vjust=-.2, hjust=1)+
  xlim(0, .06)

#color code labels 
ggtree(beast_tree)+ geom_tiplab(cex=2, color=c(rep("red", times=40), rep("blue", times=46)))+
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.9), cex=3, vjust=-.2, hjust=1)+
  xlim(0, .06)

#make manual color vector
col.vec<-substr(beast_tree@phylo$tip.label, 1,3)
col.vec[col.vec == "Tur"]<-"black"
col.vec[col.vec == "ind"]<-"palegreen2"
col.vec[col.vec == "lon"]<-"skyblue1"
col.vec[col.vec == "ste"]<-"sienna2"

#tree color coding
ggtree(beast_tree)+
  geom_tiplab(cex=2)+
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.95), cex=3, vjust=-.2, hjust=1)+
  xlim(0, .06)+
  geom_tippoint(color=col.vec)

#show node numbers
ggtree(beast_tree)+
  geom_tiplab(cex=2)+
  geom_text2(aes(label=node), 
                 cex=3, vjust=-.2, hjust=1)+
  xlim(0, .06)+
  geom_tippoint(color=col.vec)

#only view Chalcophaps without the outgroup
viewClade(ggtree(beast_tree)+
            geom_tiplab(cex=2)+
            geom_text2(aes(label=round(as.numeric(posterior), 2), 
                           subset=as.numeric(posterior) == 1), cex=3, vjust=-.2, hjust=1)+
            xlim(0, .06)+
            geom_tippoint(color=col.vec, cex=2), 
          node=89)

#rotate tree around root node
viewClade(ggtree(beast_tree) %>% rotate(89)+
            geom_tiplab(cex=2)+
            geom_text2(aes(label=round(as.numeric(posterior), 2), 
                           subset=as.numeric(posterior) == 1), cex=3, vjust=-.2, hjust=1)+
            xlim(0, .06)+
            geom_tippoint(color=col.vec, cex=2), 
          node=89)

#remove subspecies from the tip labels
beast_tree@phylo$tip.label
beast_tree@phylo$tip.label<-gsub("-.*-", "-", beast_tree@phylo$tip.label)

#fix CAS labels
beast_tree@phylo$tip.label<-gsub("longirostris-DPM530", "longirostris-CAS-DPM530", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-DPM536", "longirostris-CAS-DPM536", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JF3019", "longirostris-CAS97849", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JF3031", "longirostris-CAS-JF3031", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JF3083", "longirostris-CAS97896", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JF3084", "longirostris-CAS97897", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JF3097", "longirostris-CAS97905", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JF3110", "longirostris-CAS97910", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JF3115", "longirostris-CAS97913", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JF3117", "longirostris-CAS97915", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JF3126", "longirostris-CAS97918", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JPD453", "longirostris-CAS-JPD453", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JPD522", "longirostris-CAS-JPD522", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JPD526", "longirostris-CAS99398", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JPD777", "longirostris-CAS-JPD777", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JPD782", "longirostris-CAS-JPD782", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-JF3048", "longirostris-CAS97883", beast_tree@phylo$tip.label)
beast_tree@phylo$tip.label<-gsub("longirostris-ZRH857", "longirostris-CAS98019", beast_tree@phylo$tip.label)

#show node #s
viewClade(ggtree(beast_tree) %>% rotate(90)+
            geom_tiplab(cex=2)+
            geom_text2(aes(label=node), 
                       cex=3, vjust=-.2, hjust=1)+
            xlim(0, .06)+
            geom_tippoint(color=col.vec, cex=2, alpha=.8), 
          node=89)

#add clade labels
tree<-viewClade(ggtree(beast_tree) %>% rotate(89)+
            geom_tiplab(cex=2)+
            geom_text2(aes(label=round(as.numeric(posterior), 2), 
                           subset=as.numeric(posterior)== 1), cex=4, vjust=-.25, hjust=1.3)+
            xlim(0, .0655)+
            geom_tippoint(color=col.vec, cex=2, alpha=.8)+
            geom_cladelabel(159, paste("Solomon Islands"), offset = .0048, align = TRUE, extend = c(0.15, 0.15), fontsize =3)+
            geom_cladelabel(81, paste("New Britain"), offset = .0048, align = TRUE, extend = c(0.15, 0.15), fontsize =3)+
            geom_cladelabel(165, paste("Papua New Guinea"), offset = .0048, align = TRUE, extend = c(0.15, 0.15), fontsize =3)+
            geom_cladelabel(71, paste("Timor Island"), offset = .0048, align = TRUE, extend = c(0.15, 0.15), fontsize =3)+
            geom_cladelabel(118, paste("China, Philippines"), offset = .0048, align = TRUE, extend = c(0.15, 0.15), fontsize =3)+
            geom_cladelabel(120, paste("Christmas Island"), offset = .0048, align = TRUE, extend = c(0.15, 0.15), fontsize =3)+
            geom_cladelabel(93, paste("China,\nPhilippines,\nVietnam,\nMalaysia"), offset = .0048, align = TRUE, extend = c(0.15, 0.15), fontsize =3)+
            geom_cladelabel(156, paste("Vanuatu,\nSanta Cruz Islands"), offset = .0048, align = TRUE, extend = c(0.15, 0.15), fontsize =3)+
            geom_cladelabel(124, paste("Australia,\nPapua New Guinea,\nLouisiade Archipelago"), offset = .0048, align = TRUE, extend = c(0.15, 0.15), fontsize =3)+
            geom_treescale(fontsize=3, linesize=2, width = 1, x=.02,y=.02)
          ,node=89)

tree

grid.arrange(map,tree, nrow=1)

map
ggsave(filename = "~/Downloads/chalc.map.pdf", width = 168, height = 168, units = "mm")

tree
ggsave(filename = "~/Downloads/chalc.tree.pdf", width = 280, height = 190, units = "mm")

#knit these together in 


