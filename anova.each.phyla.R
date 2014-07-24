library(plyr)
library(reshape2)
library(ggplot2)
library(RSvgDevice)
 

##########################
#Anova on phyla differences relativized to reference plots
b.phy <- read.table(file="../bac.phyla.env", header=T, row.names=1)
f.phy <- read.table(file="../fung.phyla.csv", header=T, row.names=1)
b.phy <- b.phy[b.phy$om!="99",]
f.phy <- f.phy[f.phy$om!="99",] 
b.phy <- droplevels(b.phy)
f.phy <- droplevels(f.phy)
b.phy$site.hor <- paste(b.phy$site, b.phy$horizon, sep=".")
str(b.phy)
f.phy$site.hor <- paste(f.phy$site, f.phy$horizon, sep=".")
str(f.phy)
f.p <- f.phy
f.p$zone <- NULL
f.p$om <- NULL
f.p$horizon <- NULL
f.p$site.hor <- NULL
f.p$site <- NULL
str(f.p)

b.p <- b.phy
b.p$zone <- NULL
b.p$om <- NULL
b.p$horizon <- NULL
b.p$site.hor <- NULL
b.p$erick <- NULL
b.p$site <- NULL
str(b.p)

#average
means <- aggregate(b.p[b.phy$om==0,], by=list(b.phy$site.hor[b.phy$om==0]), FUN=mean)
means.dtf <- data.frame(means, row.names=1)
means.ext <- means.dtf[b.phy$site.hor,]
b.phy.rel <- (b.p-means.ext)/means.ext

means <- aggregate(f.p[f.phy$om==0,], by=list(f.phy$site.hor[f.phy$om==0]), FUN=mean)
means.dtf <- data.frame(means, row.names=1)
means.ext <- means.dtf[f.phy$site.hor,]
f.phy.rel <- (f.p-means.ext)/means.ext

b.p <- cbind(b.phy.rel, as.factor(b.phy$zone), as.factor(b.phy$om), as.factor(b.phy$horizon))
f.p <- cbind(f.phy.rel, as.factor(f.phy$zone), as.factor(f.phy$om), as.factor(f.phy$horizon))
str(b.p)
str(f.p)
names(b.p)[18] <- ("zone")
names(b.p)[19] <- ("om")
names(b.p)[20] <- ("horizon")
names(f.p)[16] <- ("zone")
names(f.p)[17] <- ("om")
names(f.p)[18] <- ("horizon")
str(b.p)
str(f.p)
b.p[is.na(b.p)] <- 0
b.p[b.p==Inf] <- 0
f.p[is.na(f.p)] <- 0
f.p[f.p==Inf] <- 0

b.p<- melt(b.p, id.vars=c("zone", "om", "horizon"))
f.p<- melt(f.p, id.vars=c("zone", "om", "horizon"))

#stacked bar plot of major taxa

bac.phyla.color<-c("#A6CEE3","#7DB4D5", "#5C9FC9","#3A89BD", "#1F78B4", "#B2DF8A","#79C360", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99","#ededed", "#000000")

devSVG(file="b.phyla.svg", width=8, height=10)

b.plot <- ggplot(b.p,  aes( y=value, x=factor(hor.om), fill=factor(variable), order=-as.numeric(variable)))+ #order=- makes the stack the same order as the legend
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=bac.phyla.color)+
  theme_bw()+ 
  theme(text=element_text(family="Times"))
b.plot
dev.off()

ggsave(file="b.phyla.pdf", plot=b.plot, width=8, height=10)


fung.phyla.color<-c("#deebf7", "#c6dbef",  "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b",  "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#8c2d04", "#b2df8a",  "#EdEdEd","#000000")

devSVG(file="f.phyla.svg", width=8, height=10)

f.plot <- ggplot(f.p,  aes( y=value, x=factor(hor.om), fill=factor(variable), order=-as.numeric(variable)))+ #order=- makes the stack the same order as the legend
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=fung.phyla.color)+
  theme_bw()+ 
  theme(text=element_text(family="Times"))
f.plot
dev.off()

ggsave(file="f.phyla.pdf", plot=f.plot, width=8, height=10)


# b.x <- aov(Firmicutes~om, data=b.p[b.p$horizon==1 & b.p$om==c("1", "2"),])
# summary(b.x)  #this worked



#anova on each phyla and horizon individually, have to split horizone manually because of OM3 imbalance
#organic
anova.out <- dlply(b.p[b.p$horizon==1 & b.p$om!=3,], .(variable), function(b.p) aov(value~om, data=b.p[b.p$horizon==1 & b.p$om!=3,]))
juicy.bits <- function(x)
{c(anova(x)[1,4], anova(x)[1,5])} #don't understand why this works, got it from http://oardc.osu.edu/culman/wp-content/uploads/2014/05/Reshape-and-Plyr.txt
ldply(anova.out, juicy.bits)

tukey.out <- dlply(b.p[b.p$horizon==1 & b.p$om!=3,], .(variable), function(b.p) TukeyHSD(aov(value~om, data=b.p[b.p$horizon==1 & b.p$om!=3,])))
# tukey.out
# str(tukey.out)
tuk.juice <- function(x)
{(x)$om[,c(1,4)]}
ldply(tukey.out, tuk.juice)

#mineral
anova.out <- dlply(b.p[b.p$horizon==2 ,], .(variable), function(b.p) aov(value~om, data=b.p[b.p$horizon==2,]))
juicy.bits <- function(x)
{c(anova(x)[1,4], anova(x)[1,5])}
ldply(anova.out, juicy.bits)

tukey.out <- dlply(b.p[b.p$horizon==2 ,], .(variable), function(b.p) TukeyHSD(aov(value~om, data=b.p[b.p$horizon==2 ,])))
# tukey.out
# str(tukey.out)
tuk.juice <- function(x)
{(x)$om[,c(1,4)]}
ldply(tukey.out, tuk.juice)

#organic
anova.out <- dlply(f.p[f.p$horizon==1 & f.p$om!=3,], .(variable), function(f.p) aov(value~om, data=f.p[f.p$horizon==1 & f.p$om!=3,]))
juicy.bits <- function(x)
{c(anova(x)[1,4], anova(x)[1,5])} #don't understand why this works, got it from http://oardc.osu.edu/culman/wp-content/uploads/2014/05/Reshape-and-Plyr.txt
ldply(anova.out, juicy.bits)

tukey.out <- dlply(f.p[f.p$horizon==1 & f.p$om!=3,], .(variable), function(f.p) TukeyHSD(aov(value~om, data=f.p[f.p$horizon==1 & f.p$om!=3,])))
# tukey.out
# str(tukey.out)
tuk.juice <- function(x)
{(x)$om[,c(1,4)]}
ldply(tukey.out, tuk.juice)

#mineral
anova.out <- dlply(f.p[f.p$horizon==2 ,], .(variable), function(f.p) aov(value~om, data=f.p[f.p$horizon==2,]))
juicy.bits <- function(x)
{c(anova(x)[1,4], anova(x)[1,5])}
ldply(anova.out, juicy.bits)

tukey.out <- dlply(f.p[f.p$horizon==2 ,], .(variable), function(f.p) TukeyHSD(aov(value~om, data=f.p[f.p$horizon==2 ,])))
# tukey.out
# str(tukey.out)
tuk.juice <- function(x)
{(x)$om[,c(1,4)]}
ldply(tukey.out, tuk.juice)

#attempting to graph these data-this isn't very useful.  too many frames
b <- ggplot(data=b.p, aes(x=om, y=value))+
  geom_hline(yintercept=0, aes(col="grey", linetype=2))+
  geom_boxplot()
b+  facet_grid(horizon~variable)+
#   scale_y_log10()
#   scale_fill_manual(values=c("grey", "white"))+
  theme_bw()+
  labs(title="Bacterial Phyla Relative to Reference")
ylim1 <- boxplot.stats(b.a.rel.env$invsimpson)$stats[c(1,5)]
b1 <- invsimp + coord_cartesian(ylim=ylim1*1.5)


devSVG(file="b.invsimp.svg", width=10, height=8)
invsimp <- ggplot(data=subset(b.a.rel.env, b.a.rel.env$om<98), aes(x=factor(om), y=invsimpson, fill=factor(horizon)))+
  geom_hline(yintercept=0, aes(col="grey", linetype=b))+
  geom_boxplot(outlier.size=0)+
  scale_fill_manual(values=c("grey","white"))+
  theme_bw()+
  labs(title="Bacterial Relative to Reference")
invsimp
#calculate ylim to include the wiskers but not the outliers, outliers are still used to calculate the boxes
ylim1 <- boxplot.stats(b.a.rel.env$invsimpson)$stats[c(1,5)]
b.i <- invsimp + coord_cartesian(ylim=ylim1*1.5)
dev.off()













# b.cov <- aov(coverage~factor(zone)+factor(om), data=b.a.env)
# summary(b.cov)
# # TukeyHSD(b.cov)
# 
# b.invs <- aov(invsimpson~factor(zone)*factor(om)*factor(horizon), data=b.a.env)
# summary(b.invs)
# # TukeyHSD(b.invs)
# b.invs <- aov(invsimpson~factor(zone)*factor(om)*factor(horizon), data=b.a.rel.env)
# summary(b.invs)
# 
# f.invs <- aov(invsimpson~factor(om)+factor(zone), data=subset(f.a.env, horizon==1& om>=1 &om<3))
# summary(f.invs)
# # TukeyHSD(f.invs)
# f.invs <- aov(invsimpson~factor(om)+factor(zone), data=subset(f.a.rel.env, horizon==1& om>=1 &om<3))
# summary(f.invs)
# # TukeyHSD(f.invs)
# 
# f.invs <- aov(invsimpson~factor(om)+factor(zone), data=subset(f.a.env, horizon==2 & om >=1))
# summary(f.invs)
# TukeyHSD(f.invs)
# f.invs <- aov(invsimpson~factor(zone)*factor(om), data=subset(f.a.rel.env, horizon==2 & om>=1))
# summary(f.invs)
# TukeyHSD(f.invs)
# 
# f.invs <- aov(invsimpson~factor(om), data=subset(f.a.rel.env, horizon==2 & om>=1))
# summary(f.invs)
# TukeyHSD(f.invs)
# 
# f.invs <- aov(invsimpson~factor(om), data=subset(f.a.env, horizon==2 & om>=1))
# summary(f.invs)
# TukeyHSD(f.invs)
# 
# b.even <- aov(shannoneven~factor(zone)+factor(om), data=subset(b.a.env, horizon==1))
# summary(b.even)
# # TukeyHSD(b.even)
# 
# b.sobs <- aov(sobs~factor(zone)+factor(om), data=b.a.env)
# summary(b.sobs)
# # TukeyHSD(b.sobs)
# 
# 
# f.cov <- aov(coverage~factor(zone)+factor(om), data=f.a.env)
# summary(f.cov)
# # TukeyHSD(f.cov)
# 
# f.invs <- aov(invsimpson~factor(zone)+factor(om), data=f.a.env)
# summary(f.invs)
# # TukeyHSD(f.invs)
# 
# f.even <- aov(shannoneven~factor(zone)+factor(om), data=f.a.env)
# summary(f.even)
# # TukeyHSD(f.even)
# 
# f.sobs <- aov(sobs~factor(zone)+factor(om), data=f.a.env)
# summary(f.sobs)
# # TukeyHSD(f.sobs)
# 
# b.invs <- aov(invsimpson~factor(zone)*factor(om), data=subset(b.a.env, horizon==1& om>=1 &om<3))
# summary(b.invs)
# # TukeyHSD(b.invs)
# b.invs <- aov(invsimpson~factor(zone)*factor(om), data=subset(b.a.rel.env, horizon==1& om>=1 &om<3))
# summary(b.invs)
# # TukeyHSD(b.invs)
# 
# b.invs <- aov(invsimpson~factor(zone)*factor(om), data=subset(b.a.env, horizon==2 & om >=1))
# summary(b.invs)
# # TukeyHSD(b.invs)
# b.invs <- aov(invsimpson~factor(zone)*factor(om), data=subset(b.a.rel.env, horizon==2 & om>=1))
# summary(b.invs)
# # TukeyHSD(b.invs)
# 
# 