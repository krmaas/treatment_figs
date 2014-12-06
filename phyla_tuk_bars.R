library(ggplot2)
library(RSvgDevice)
library(reshape2)
library(plyr)


#read in Tukey HSD sig treatment pairs
tuk <- read.table(file="tuk.zone.phyla.csv", header=T)
tuk$phyla <- factor(tuk$phyla, levels=c("Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria", 
                                                   "Deltaproteobacteria", "Proteobacteria", "Acidobacteria", "Actinobacteria", 
                                                   "Bacteroidetes", "Chloroflexi", "Cyanobacteria", "Firmicutes", "Gemmatimonadetes",
                                                   "Planctomycetes", "Verrucomicrobia", "unclassified", "other",
                                                   "Asc_Dothideomycetes", "Asc_Eurotiomycetes", "Asc_Leotiomycetes", "Asc_Pezizomycetes", 
                                                   "Asc_Saccharomycetes", "Asc_Sordariomycetes", "Asc_other", "Bas_Agaricomycetes" , 
                                                   "Bas_Microbotryomycetes" , "Bas_Tremellomycetes", "Bas_Wallemiomycetes", 
                                                   "Bas_other", "Zygomycota", "Unclassified", "Other_Fungi" ), ordered=TRUE)
# str(tuk)

tuk.pos <- tuk[tuk$diff>0,]
tuk.neg <- tuk[tuk$diff<0,]


theme_set(theme_bw(16))
# 
phyla.color <- c("Asc_Dothideomycetes" = "#deebf7", "Asc_Eurotiomycetes" = "#c6dbef", "Asc_Leotiomycetes" =  "#6baed6", 
                 "Asc_Pezizomycetes" = "#4292c6", "Asc_Saccharomycetes" = "#2171b5", "Asc_Sordariomycetes" = "#08519c", 
                 "Asc_other" = "#08306b", "Bas_Agaricomycetes" = "#fdae6b", "Bas_Microbotryomycetes" = "#fd8d3c", 
                 "Bas_Tremellomycetes" = "#f16913", "Bas_Wallemiomycetes" = "#d94801", "Bas_other" = "#8c2d04", 
                 "Zygomycota" = "#b2df8a", "Unclassified" =  "#EdEdEd", "Other_Fungi" = "#000000", 
                 "Alphaproteobacteria" = "#A6CEE3", "Betaproteobacteria" = "#7DB4D5", "Deltaproteobacteria" = "#5C9FC9", 
                 "Gammaproteobacteria" = "#3A89BD", "Proteobacteria" = "#1F78B4", "Actinobacteria" = "#B2DF8A", 
                 "Acidobacteria" = "#79C360", "Chloroflexi" = "#33A02C", "Cyanobacteria" = "#FB9A99", "Bacteroidetes" = "#E31A1C", 
                 "Firmicutes" = "#FDBF6F", "Gemmatimonadetes" = "#FF7F00", "Planctomycetes" = "#CAB2D6", 
                 "Verrucomicrobia" = "#6A3D9A", "Nitrospirae" = "#FFFF99", "unclassified" = "#ededed", "other" = "#000000")

p <- ggplot()+
  geom_hline(yintercept=0, colour="grey")+
  geom_bar(data = tuk.pos, aes (x= factor (treatment, levels=c("0-1", "0-2", "1-2", "0-3", "1-3", "2-3")),
                                y=diff, fill=phyla), stat="identity")+
  geom_bar(data = tuk.neg, aes (x= factor (treatment, levels=c("0-1", "0-2", "1-2", "0-3", "1-3", "2-3")), 
                                y=diff, fill=phyla), stat="identity") +
  facet_grid(domain~horizon, scales="free_y")+
  scale_fill_manual(values=phyla.color, breaks=levels(tuk$phyla))+
  labs(x="Treatment Pair", y="Significant Differences")
# +
#   coord_flip()
p  


devSVG(file="phyla.tuk.horizontal.svg", width=10, height=8)

  
# p <- ggplot()+
#   geom_hline(yintercept=0, colour="grey")+
#   geom_bar(data= tuk.pos[order(tuk.pos$phyla),], 
#            aes(x=factor(treatment,levels=c("0-1", "0-2", "1-2", "0-3", "1-3", "2-3")), y=diff, fill=phyla), stat="identity")+
#   geom_bar(data= tuk.neg[order(tuk.neg$phyla),], 
#            aes(x=factor(treatment,levels=c("0-1", "0-2", "1-2", "0-3", "1-3", "2-3")), y=diff, fill=phyla), stat="identity") +
#   facet_grid(domain~horizon, scales="free_y")+
#   scale_fill_manual(values=phyla.color)+
#   labs(x="Treatment Pair", y="Significant Differences")+
#   coord_flip()
# p  

dev.off()

ggsave(file="phyla.tuk.horizontal.pdf", plot=p, width=10, height=8)



bac.fung.tuk.color<-c("#A6CEE3","#7DB4D5", "#5C9FC9","#3A89BD", "#1F78B4", "#B2DF8A","#79C360", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#ededed", "#000000","#deebf7",   "#6baed6", "#4292c6", "#08519c", "#08306b",  "#fdae6b",  "#f16913", "#d94801", "#8c2d04",  "#b2df8a", "#EdEdEd","#000000")

zone.color <- c("bIDF"= "#A6CEE3", "bSBS"= "#1F78B4", "MD"= "#FB9A99","oBS"= "#B2DF8A", "oJP" = "#33A02C", "TX" ="#E31A1C")

# @rgb(166,206,227)  IDF  #A6CEE3
# @rgb(31,120,180)	SBS #1F78B4
# @rgb(251,154,153)	PP  #FB9A99
# @rgb(178,223,138)	BS  #B2DF8A
# @rgb(51,160,44)	JP  #33A02C
# @rgb(227,26,28)	LP  #E31A1C