####regress relativized shift in phyla abundance from ref agains standardized C, N, CN, pH


setwd('/media/data/R/treatment_fig/')

library(ggplot2)
library(RSvgDevice)
library(reshape2)
library(plyr)
library(dplyr)

zone.color <- c("bIDF"= "#A6CEE3", "bSBS"= "#1F78B4", "MD"= "#FB9A99","oBS"= "#B2DF8A", "oJP" = "#33A02C", "TX" ="#E31A1C")
theme_set(theme_bw())

#read in and sort env data
f.env <- read.table(file="../f.all.env.csv", header=T)


#create site.hor variable to melt on
f.env$site.hor <- factor(paste(f.env$site, f.env$horizon, sep="."))
f.toz <- melt(f.env, id.vars=c("sample","zone", "plot", "site", "om", "horizon", "site.hor"))

#calculate z scores for each parameter
f.z <- ddply(f.toz, .(site.hor, variable), function(x) {x$zscore  = scale(x$value) 
                                                        return(x)})
#turn back into wide matrix with zscore rather than raw data
f.z.w <- dcast(f.z, sample+zone+plot+site+om+horizon+site.hor~variable, value.var="zscore")

#read in the phyla abundances
f.phy <- read.table(file="../fung.phyla.csv", header=T, row.names=1)
#must read in sample names as row names for aggregate means proceedure to work
f.phy <- f.phy[order(row.names(f.phy)),]

f.z.w <- f.z.w[order(f.z.w$sample),]
f.phy$zone <- NULL
f.phy$om <- NULL
f.phy$horizon <- NULL
f.phy$site <- NULL

str(f.phy)
str(f.z.w)

#####   Relativize taxa abundances as shifts from ref for each site.hor

means <- aggregate(f.phy[f.z.w$om==0,], by=list(f.z.w$site.hor[f.z.w$om==0]), FUN=mean)
warnings()
means.dtf <- data.frame(means, row.names=1)
means.ext <- means.dtf[f.z.w$site.hor,]
f.phy.rel <- (f.phy-means.ext)/means.ext
f.phy.rel[f.phy.rel==Inf] <- 1


str(f.phy.rel)

f.p.rel <- cbind(factor(f.z.w$zone), factor(f.z.w$horizon), factor(f.z.w$om), f.z.w$C, f.z.w$N, f.z.w$CN, f.z.w$pH_H2O, f.phy.rel)
names(f.p.rel)[1] <- ("zone")
names(f.p.rel)[2] <- ("horizon")
names(f.p.rel)[3] <- ("om")
names(f.p.rel)[4] <- ("C")
names(f.p.rel)[5] <- ("N")
names(f.p.rel)[6] <- ("CN")
names(f.p.rel)[7] <- ("pH")

str(f.p.rel)

# remove all samples that have no env data (this takes care of fire as well)
f.p.na <- f.p.rel[!is.na(f.p.rel$C),]


#write function to use with plyr to return a dataframe of intercept, slope, and R2 for each ecozone, each horizon, each om
f.p.na$z.h.om <- paste(f.p.na$zone, f.p.na$horizon, f.p.na$om, sep=".")


#need to remove the one weird JP org OM3 sample because it's throwing NAs
f.p.na <- droplevels(f.p.na[f.p.na$z.h.om!="oJP.1.3",])

f.p <- melt(f.p.na, id.vars=c("zone", "horizon", "om", "C", "N", "CN", "pH", "z.h.om"))
str(f.p)
#remove all taxa instances that have NA
f.p <- f.p[!is.na(f.p$value),]

lmfun <- function(x)   {
  lmfit <- lm(value~C, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

f.C.p <- ddply(f.p, ~zone+horizon+om+variable, lmfun)

C <- ggplot(data = f.C.p[f.C.p$rsq>0.25 & f.C.p$horizon==2,], aes(x=om, y=slope, color=factor(zone)))+
  facet_wrap( ~ variable, ncol=8)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Fungal taxa ~ C")+ labs(title="Fungal relativized abundances, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)+
#   ylim(-10,10)+
  geom_hline(y=0)

lmfun <- function(x)   {
  lmfit <- lm(value~N, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

f.N.p <- ddply(f.p, ~zone+horizon+om+variable, lmfun)

N <- ggplot(data = f.N.p[f.N.p$rsq>0.25 & f.N.p$horizon==2,], aes(x=om, y=slope, color=factor(zone)))+
  facet_wrap( ~ variable, ncol=8)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Fungal taxa ~ N")+ labs(title="Fungal relativized abundances, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)+
  #   ylim(-10,10)+
  geom_hline(y=0)


lmfun <- function(x)   {
  lmfit <- lm(value~CN, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

f.CN.p <- ddply(f.p, ~zone+horizon+om+variable, lmfun)

CN <- ggplot(data = f.CN.p[f.CN.p$rsq>0.25 & f.CN.p$horizon==2,], aes(x=om, y=slope, color=factor(zone)))+
  facet_wrap( ~ variable, ncol=8)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Fungal taxa ~ N")+ labs(title="Fungal relativized abundances, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)+
  #   ylim(-10,10)+
  geom_hline(y=0)


lmfun <- function(x)   {
  lmfit <- lm(value~pH, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

f.pH.p <- ddply(f.p, ~zone+horizon+om+variable, lmfun)

pH <- ggplot(data = f.pH.p[f.pH.p$rsq>0.25 & f.pH.p$horizon==2,], aes(x=om, y=slope, color=factor(zone)))+
  facet_wrap( ~ variable, ncol=8)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Fungal taxa ~ pH")+ labs(title="Fungal relativized abundances, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)+
#     ylim(-5,5)+
  geom_hline(y=0)

C 
ggsave(file="f.taxa.C.slopes.svg", height=12, width=16)

CN
ggsave(file="f.taxa.CN.slopes.svg", height=12, width=16)

N
ggsave(file="f.taxa.N.slopes.svg", height=12, width=16)

pH
ggsave(file="f.taxa.pH.slopes.svg", height=12, width=16)


#read in and sort env data
b.env <- read.table(file="../b.all.env.csv", header=T)


#create site.hor variable to melt on
b.env$site.hor <- factor(paste(b.env$site, b.env$horizon, sep="."))
b.toz <- melt(b.env, id.vars=c("sample","zone", "plot", "site", "om", "horizon", "site.hor"))

#calculate z scores for each parameter
b.z <- ddply(b.toz, .(site.hor, variable), function(x) {x$zscore  = scale(x$value) 
                                                        return(x)})
#turn back into wide matrix with zscore rather than raw data
b.z.w <- dcast(b.z, sample+zone+plot+site+om+horizon+site.hor~variable, value.var="zscore")

#read in the phyla abundances
b.phy <- read.table(file="../bac.phyla.csv", header=T, row.names=1, sep=",")
b.phy <- b.phy[order(row.names(b.phy)),]

b.z.w <- b.z.w[order(b.z.w$sample),]
b.phy$zone <- NULL
b.phy$om <- NULL
b.phy$horizon <- NULL
b.phy$site <- NULL

str(b.phy)
str(b.z.w)

#####   Relativize taxa abundances as shifts from ref for each site.hor

means <- aggregate(b.phy[b.z.w$om==0,], by=list(b.z.w$site.hor[b.z.w$om==0]), FUN=mean)
means.dtf <- data.frame(means, row.names=1)
means.ext <- means.dtf[b.z.w$site.hor,]
b.phy.rel <- (b.phy-means.ext)/means.ext
b.phy.rel[b.phy.rel==Inf] <- 1


str(b.phy.rel)

b.p.rel <- cbind(factor(b.z.w$zone), factor(b.z.w$horizon), factor(b.z.w$om), b.z.w$C, b.z.w$N, b.z.w$CN, b.z.w$pH_H2O, b.phy.rel)
names(b.p.rel)[1] <- ("zone")
names(b.p.rel)[2] <- ("horizon")
names(b.p.rel)[3] <- ("om")
names(b.p.rel)[4] <- ("C")
names(b.p.rel)[5] <- ("N")
names(b.p.rel)[6] <- ("CN")
names(b.p.rel)[7] <- ("pH")

str(b.p.rel)

# remove all samples that have no env data (this takes care of fire as well)
b.p.na <- b.p.rel[!is.na(b.p.rel$C),]


#write function to use with plyr to return a dataframe of intercept, slope, and R2 for each ecozone, each horizon, each om
b.p.na$z.h.om <- paste(b.p.na$zone, b.p.na$horizon, b.p.na$om, sep=".")


#need to remove the one weird JP org OM3 sample because it's throwing NAs
b.p.na <- droplevels(b.p.na[b.p.na$z.h.om!="oJP.1.3",])

b.p <- melt(b.p.na, id.vars=c("zone", "horizon", "om", "C", "N", "CN", "pH", "z.h.om"))
str(b.p)
#remove all taxa instances that have NA
b.p <- b.p[!is.na(b.p$value),]

lmfun <- function(x)   {
  lmfit <- lm(value~C, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

b.C.p <- ddply(b.p, ~zone+horizon+om+variable, lmfun)

C <- ggplot(data = b.C.p[b.C.p$rsq>0.25 & b.C.p$horizon==2,], aes(x=om, y=slope, color=factor(zone)))+
  facet_wrap( ~ variable, ncol=9)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Bacterial taxa ~ C")+ labs(title="Bacterial relativized abundances, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)+
  #   ylim(-10,10)+
  geom_hline(y=0)
C

lmfun <- function(x)   {
  lmfit <- lm(value~N, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

b.N.p <- ddply(b.p, ~zone+horizon+om+variable, lmfun)

N <- ggplot(data = b.N.p[b.N.p$rsq>0.25 & b.N.p$horizon==2,], aes(x=om, y=slope, color=factor(zone)))+
  facet_wrap( ~ variable, ncol=9)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Bacterial taxa ~ N")+ labs(title="Bacterial relativized abundances, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)+
  #   ylim(-10,10)+
  geom_hline(y=0)


lmfun <- function(x)   {
  lmfit <- lm(value~CN, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

b.CN.p <- ddply(b.p, ~zone+horizon+om+variable, lmfun)

CN <- ggplot(data = b.CN.p[b.CN.p$rsq>0.25 & b.CN.p$horizon==2,], aes(x=om, y=slope, color=factor(zone)))+
  facet_wrap( ~ variable, ncol=9)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Bacterial taxa ~ CN")+ labs(title="Bacterial relativized abundances, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)+
  #   ylim(-10,10)+
  geom_hline(y=0)


lmfun <- function(x)   {
  lmfit <- lm(value~pH, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

b.pH.p <- ddply(b.p, ~zone+horizon+om+variable, lmfun)

pH <- ggplot(data = b.pH.p[b.pH.p$rsq>0.25 & b.pH.p$horizon==2,], aes(x=om, y=slope, color=factor(zone)))+
  facet_wrap( ~ variable, ncol=9)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Bacterial taxa ~ pH")+ labs(title="Bacterial relativized abundances, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)+
  #     ylim(-5,5)+
  geom_hline(y=0)

C
ggsave(file="b.taxa.C.slopes.svg", height=12, width=16)

CN
ggsave(file="b.taxa.CN.slopes.svg", height=12, width=16)

N
ggsave(file="b.taxa.N.slopes.svg", height=12, width=16)

pH
ggsave(file="b.taxa.pH.slopes.svg", height=12, width=16)
