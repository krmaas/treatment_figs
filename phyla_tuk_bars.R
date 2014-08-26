library(ggplot2)
library(RSvgDevice)
library(reshape2)
library(plyr)


#read in Tukey HSD sig treatment pairs
tuk <- read.table(file="tuk.zone.phyla.csv", header=T)
# 
# tuk$diff.b <- tuk$diff
# tuk$diff.b[tuk$diff.b>=0] <- 0
# tuk$diff.b[tuk$diff.b<0] <- "n"
# tuk$diff.b[tuk$diff.b=="0"] <- "p"
# 
# str(tuk)
# ###unnecessary to do this twice, should have just reordered
# # tuk$treat.dir <- paste(tuk$treatment, tuk$diff.b, sep=".")
# # tuk$dir.treat <- paste(tuk$diff.b, tuk$treatment, sep=".")
# tuk$dir.treat <- factor(paste(tuk$diff.b, tuk$treatment, sep="."), levels=c("p.0-1", "p.0-2", "p.1-2", "p.0-3", "p.1-3", "p.2-3", "n.0-1", "n.0-2", "n.1-2", "n.0-3", "n.1-3", "n.2-3"))
# 
# theme_set(theme_bw(16))
# 
# ggplot(tuk, aes(x=treatment, fill=zone))+
#   geom_histogram()+
#   facet_grid(domain~diff.b~horizon~phyla)+
#   scale_fill_manual(values=zone.color)
# 
# ggplot(tuk, aes(x=dir.treat, fill=zone))+
#   geom_histogram()+
#   facet_grid(horizon~domain)+
#   scale_fill_manual(values=zone.color)
# 
# ggplot(tuk[order(tuk$zone),], aes(x=dir.treat, y=diff, fill=zone))+
#   geom_bar(stat="identity")+
#   facet_grid(horizon~domain)+
#   scale_fill_manual(values=zone.color)
# 
# 
# ggplot(tuk[order(tuk$zone),], aes(x=factor(treatment, levels=c("0-1", "0-2", "1-2", "0-3", "1-3", "2-3")), y=abs(diff), fill=zone))+
#   geom_bar(stat="identity")+
#   facet_grid(horizon~domain)+
#   scale_fill_manual(values=zone.color)
# 
# tuk.pos <- tuk[tuk$diff>0,]
# tuk.neg <- tuk[tuk$diff<0.000000001,]
# 
# p <- ggplot(tuk[order(tuk$phyla),], aes(x=factor(dir.treat,levels=c("p.0-1", "p.0-2", "p.1-2", "p.0-3", "p.1-3", "p.2-3", "n.0-1", "n.0-2", "n.1-2", "n.0-3", "n.1-3", "n.2-3")), y=diff, fill=phyla))+
#   geom_bar(stat="identity")+
#   facet_grid(horizon~domain)+
#   scale_fill_manual(values=bac.fung.tuk.color)
# 
# p+guides(col=guide_legend(nrow=3))

p <- ggplot()+
  geom_hline(yintercept=0, colour="grey")+
  geom_bar(data= tuk.pos[order(tuk.pos$phyla),], aes(x=factor(treatment,levels=c("0-1", "0-2", "1-2", "0-3", "1-3", "2-3")), y=diff, fill=phyla), stat="identity")+
  geom_bar(data= tuk.neg[order(tuk.neg$phyla),], aes(x=factor(treatment,levels=c("0-1", "0-2", "1-2", "0-3", "1-3", "2-3")), y=diff, fill=phyla), stat="identity")+
  facet_grid(domain~horizon, scales="free_y")+
  scale_fill_manual(values=bac.fung.tuk.color, guide=guide_legend(ncol=2))
p  


bac.fung.tuk.color<-c("#A6CEE3","#7DB4D5", "#5C9FC9","#3A89BD", "#1F78B4", "#B2DF8A","#79C360", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#ededed", "#000000","#deebf7",   "#6baed6", "#4292c6", "#08519c", "#08306b",  "#fdae6b",  "#f16913", "#d94801", "#8c2d04",  "#b2df8a", "#EdEdEd","#000000")

zone.color <- c("#A6CEE3","#1F78B4", "#FB9A99","#B2DF8A", "#33A02C","#E31A1C")

# @rgb(166,206,227)  IDF  #A6CEE3
# @rgb(31,120,180)	SBS #1F78B4
# @rgb(251,154,153)	PP  #FB9A99
# @rgb(178,223,138)	BS  #B2DF8A
# @rgb(51,160,44)	JP  #33A02C
# @rgb(227,26,28)	LP  #E31A1C