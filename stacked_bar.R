setwd("/media/data/R/treatment_fig/")
library(ggplot2)
library(RSvgDevice)
library(reshape)

theme_set(theme_bw(16))

b.phy <- read.table(file="../bac.phyla.env", header=T, row.names=1)
f.phy <- read.table(file="../fung.phyla.csv", header=T, row.names=1)
b.phy <- b.phy[b.phy$om!="99",]
f.phy <- f.phy[f.phy$om!="99",] 
b.phy <- droplevels(b.phy)
f.phy <- droplevels(f.phy)
b.phy$hor.om <- paste(b.phy$horizon, b.phy$om, sep=".")
f.phy$hor.om <- paste(f.phy$horizon, f.phy$om, sep=".")
# str(b.phy)


b.p<- melt(b.phy, id.vars=c("zone", "om", "horizon", "site", "erick", "hor.om"))
f.p<- melt(f.phy, id.vars=c("zone", "om", "horizon", "site", "hor.om"))
str(b.p)
str(f.p)

bac.phyla.color<-c("Alphaproteobacteria" = "#A6CEE3", "Betaproteobacteria" = "#7DB4D5", "Deltaproteobacteria" = "#5C9FC9", "Gammaproteobacteria" = "#3A89BD", "Proteobacteria" = "#1F78B4", "Actinobacteria" = "#B2DF8A", "Acidobacteria" = "#79C360", "Bacteroidetes" = "#33A02C", "Firmicutes" = "#FB9A99", "Chloroflexi" = "#E31A1C", "Cyanobacteria" = "#FDBF6F", "Gemmatimonadetes" = "#FF7F00", "Planctomycetes" = "#CAB2D6", "Verrucomicrobia" = "#6A3D9A", "Nitrospirae" = "#FFFF99", "other" = "#ededed", "unclassified" = "#000000")

b.p2 <- aggregate(value~variable+hor.om, data=b.p, FUN=sum) #ggplot aggregation makes a box for each sample=horizontal lines in bars

devSVG(file="b.phyla.svg", width=8, height=10)

b.plot <- ggplot(b.p2,  aes( y=value, x=factor(hor.om), color=NULL, fill=factor(variable), order=-as.numeric(variable)))+ #order=- makes the stack the same order as the legend
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=bac.phyla.color)
b.plot
dev.off()

ggsave(file="b.phyla.pdf", plot=b.plot, width=8, height=10)



### fung phyla
fung.phyla.color<-c("Asc_Dothideomycetes" = "#deebf7", "Asc_Eurotiomycetes" = "#c6dbef", "Asc_Leotiomycetes" =  "#6baed6", "Asc_Pezizomycetes" = "#4292c6", "Asc_Saccharomycetes" = "#2171b5", "Asc_Sordariomycetes" = "#08519c", "other_Ascomycota" = "#08306b", "Bas_Agaricomycetes" = "#fdae6b", "Bas_Microbotryomycetes" = "#fd8d3c", "Bas_Tremellomycetes" = "#f16913", "Bas_Wallemiomycetes" = "#d94801", "other_Basidiomycota" = "#8c2d04", "Zygomycota" = "#b2df8a", "Unclassified" =  "#EdEdEd", "Other_Fungi" = "#000000")

f.p2 <- aggregate(value~variable+hor.om, data=f.p, FUN=sum)

devSVG(file="f.phyla1.svg", width=8, height=10)

f.plot <- ggplot(f.p,  aes( y=value, x=factor(hor.om), fill=factor(variable), order=-as.numeric(variable)))+ #order=- makes the stack the same order as the legend
  geom_bar(position="fill", stat="identity", linetype=0)+
  scale_fill_manual(values=fung.phyla.color)
#   theme(axis.text.x=element_text(angle=45, hjust=1))
f.plot
dev.off()

ggsave(file="f.phyla.pdf", plot=f.plot, width=8, height=10)



################
#random plotting code, should delete
#colors for NMS
#BS   "#A6CEE3"
#IDF  "#B2DF8A"
#JP   "#1F78B4"
#MD   "#FB9A99"
#SBS  "#33A02C"
#TX   "#E31A1C"
#scatterplot 
# ggplot(b.disp)+geom_point(mapping=aes(x=om, y=zone.om.hor, colour=zone), size=3, position='jitter')+scale_colour_manual(values=c("#B2DF8A","#33A02C","#FB9A99","#A6CEE3","#1F78B4","#E31A1C"))

# ggplot(b.phyla.om.long, aes(y=value, x=factor(om), fill=factor(variable)))+
#   geom_bar(position="fill", stat="identity")+
#   scale_fill_manual(values=bac.phyla.color)

# scale_y_continuous(limits=c(0.2, 0.85))+

# ggplot(b.disp)+geom_point(mapping=aes(x=om, y=zone.om.hor, colour=zone), size=3)+scale_colour_manual(values=c("#B2DF8A","#33A02C","#FB9A99","#A6CEE3","#1F78B4","#E31A1C"))+ggtitle("Bacterial dispersion, zone+om+horizon as group")+facet_grid(~zone)
