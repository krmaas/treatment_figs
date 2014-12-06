##### I want to compare changes in diversity with changes in environmental variables.  I'm standardizing C, N, CN, and pH 
#     around the mean and relativatizing shifts in alpha diversity from the reference baseline.  
#####
##I'm not totally happly with the plots-I like the thick line around the points showing R2 but can't figure out how to do
# and jitter the points


setwd('/media/data/R')

library(ggplot2)
library(RSvgDevice)
library(reshape2)
library(plyr)
library(dplyr)

zone.color <- c("bIDF"= "#A6CEE3", "bSBS"= "#1F78B4", "MD"= "#FB9A99","oBS"= "#B2DF8A", "oJP" = "#33A02C", "TX" ="#E31A1C")
theme_set(theme_bw())

#read in and sort env data
f.env <- read.table(file="f.all.env.csv", header=T)
f.env <- f.env[order(f.env$sample),]
#create site.hor variable to melt on
f.env$site.hor <- paste(f.env$site, f.env$horizon, sep=".")
f.toz <- melt(f.env, id.vars=c("sample","zone", "plot", "site", "om", "horizon", "site.hor"))
#calculate z scores for each parameter
f.z <- ddply(f.toz, .(site.hor, variable), function(x) {x$zscore  = scale(x$value) 
                                                        return(x)})
#turn back into wide matrix with zscore rather than raw data
f.z.w <- dcast(f.z, sample+zone+plot+site+om+horizon+site.hor~variable, value.var="zscore")

#read in alpha numbers
f.a <- read.table(file="fung.alpha.cvs", header=T, row.names=1)
f.a <- f.a[order(row.names(f.a)),]
f.a.env <- cbind(f.a, f.env$zone, f.env$horizon, f.env$om)
names(f.a.env)[7] <- ("zone")
names(f.a.env)[8] <- ("horizon")
names(f.a.env)[9] <- ("om")

#cal mean ref values
means <- aggregate(f.a[f.env$om==0,], by=list(f.env$site.hor[f.env$om==0]), FUN=mean)
means.dtf <- data.frame(means, row.names=1)
means.ext <- means.dtf[f.env$site.hor,]
f.a.rel.sub <- (f.a-means.ext)/means.ext

f.a.rel.sub$sample <- row.names(f.a.rel.sub)

f.a.z <- merge(f.a.rel.sub, f.z.w, by="sample")
str(f.a.z)

f.a.z <- f.a.z[f.a.z$om<98,]


#read in and sort env data
b.env <- read.table(file="b.all.env.csv", header=T)
b.env <- b.env[order(b.env$sample),]
#create site.hor variable to melt on
b.env$site.hor <- paste(b.env$site, b.env$horizon, sep=".")
b.toz <- melt(b.env, id.vars=c("sample","zone", "plot", "site", "om", "horizon", "site.hor"))
#calculate z scores for each parameter
b.z <- ddply(b.toz, .(site.hor, variable), function(x) {x$zscore  = scale(x$value) 
                                                        return(x)})
#turn back into wide matrix with zscore rather than raw data
b.z.w <- dcast(b.z, sample+zone+plot+site+om+horizon+site.hor~variable, value.var="zscore")

#read in alpha numbers
b.a <- read.table(file="bac.alpha.csv", header=T, row.names=1)
b.a <- b.a[order(row.names(b.a)),]
b.a.env <- cbind(b.a, b.env$zone, b.env$horizon, b.env$om)
names(b.a.env)[7] <- ("zone")
names(b.a.env)[8] <- ("horizon")
names(b.a.env)[9] <- ("om")

#cal mean ref values
means <- aggregate(b.a[b.env$om==0,], by=list(b.env$site.hor[b.env$om==0]), FUN=mean)
means.dtf <- data.frame(means, row.names=1)
means.ext <- means.dtf[b.env$site.hor,]
b.a.rel.sub <- (b.a-means.ext)/means.ext

b.a.rel.sub$sample <- row.names(b.a.rel.sub)

b.a.z <- merge(b.a.rel.sub, b.z.w, by="sample")
str(b.a.z)

b.a.z <- b.a.z[b.a.z$om<98,]


###unfortuantly shows nothing other than abstract art
# p <- ggplot(data = f.a.z, aes(x=invsimpson, y=C, color=factor(om), shape=factor(zone)))+
#   geom_point(size=2)+
#   facet_grid(horizon ~ .)+
#   xlim(-1, 2.5)+
#   ylim(-2.5, 2.5)+
#   scale_color_manual(values=c("green4", "yellowgreen", "yellow3", "tan4"))
# p
# p +  stat_smooth(method="lm", size=2)
# ggsave(file="f.invsimp_C_zone.pdf")
# 
# p <- ggplot(data = b.a.z, aes(x=invsimpson, y=CN, color=factor(om), shape=zone))+
#   geom_point(size=2)+
#   facet_grid(horizon ~ .)+
#   xlim(-1, 2)+
#   ylim(-2.5, 2.5)+
#   scale_color_manual(values=c("green4", "yellowgreen", "yellow3", "tan4"))
# p +  stat_smooth(method="lm", size=2)
# ggsave(file="b.invsimp_C.zone.pdf", height=10, width=8)

#write function to use with plyr to return a dataframe of intercept, slope, and R2 for each ecozone, each horizon, each om
f.a.z$z.h.om <- paste(f.a.z$zone, f.a.z$horizon, f.a.z$om, sep=".")

str(f.a.z)


#need to remove the one weird JP org OM3 sample because it's throwing NAs
f.a.z <- droplevels(f.a.z[f.a.z$z.h.om!="oJP.1.3",])

#remove all samples that have no env data                        
f.a.z.na <- f.a.z[!is.na(f.a.z$C),]



#can't figure out how to specify x and y in the ddply command, have to replace in function each call

lmfun <- function(x)   {
  lmfit <- lm(invsimpson~C, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

f.C.invsimps <- ddply(f.a.z.na, ~zone+horizon+om, lmfun)


C <- ggplot(data = f.C.invsimps, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Inv. Simpson ~ C")+ labs(title="Fungal relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)
    

lmfun <- function(x)   {
  lmfit <- lm(invsimpson~CN, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

f.CN.invsimps <- ddply(f.a.z.na, ~zone+horizon+om, lmfun)

CN <- ggplot(data = f.CN.invsimps, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Inv. Simpson ~ CN")+ labs(title="Fungal relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)

lmfun <- function(x)   {
  lmfit <- lm(invsimpson~N, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

f.N.invsimps <- ddply(f.a.z.na, ~zone+horizon+om, lmfun)


N <- ggplot(data = f.N.invsimps, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Inv. Simpson ~ N")+ labs(title="Fungal relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)

lmfun <- function(x)   {
  lmfit <- lm(invsimpson~pH_H2O, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

f.pH.invsimps <- ddply(f.a.z.na, ~zone+horizon+om, lmfun)

pH <- ggplot(data = f.pH.invsimps, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Inv. Simpson ~ pH")+ labs(title="Fungal relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)
pH

########not sure why this isn't working
# 
# test <- pH <- ggplot(data = f.pH.invsimps, aes(x=om, y=slope, fill=factor(zone), color=rsq))+
#   facet_grid(horizon ~ .)+
#   scale_fill_manual(values=zone.color)+
#   xlab("OM")+ ylab("Slope Inv. Simpson ~ pH")+ labs(title="Fungal relativized diversity, standardized env")+
#   geom_point2(aes(size=4.5, thickness=2), shape=21, position=position_jitter(w=0.1, h=0))+
#   scale_color_gradient(high="black", low="grey50")
# test
# devSVG(file="test.svg", height=10, width=8)
# test
# dev.off()
# #source multiplot
# multiplot(C, N, CN, pH, cols=2)

#looks great in RStudio doesn't save well

devSVG(file="f.C.slopes.svg", height=10, width=8)
C
dev.off()
devSVG(file="f.CN.slopes.svg", height=10, width=8)
CN
dev.off()
devSVG(file="f.N.slopes.svg", height=10, width=8)
N
dev.off()
devSVG(file="f.pH.slopes.svg", height=10, width=8)
pH
dev.off()



###again in the key of Bacteria

#write function to use with plyr to return a dataframe of intercept, slope, and R2 for each ecozone, each horizon, each om
b.a.z$z.h.om <- paste(b.a.z$zone, b.a.z$horizon, b.a.z$om, sep=".")

str(b.a.z)


#need to remove the one weird JP org OM3 sample because it's throwing NAs
b.a.z <- droplevels(b.a.z[b.a.z$z.h.om!="oJP.1.3",])

#remove all samples that have no env data                        
b.a.z.na <- b.a.z[!is.na(b.a.z$C),]



#can't figure out how to specify x and y in the ddply command, have to replace in function each call

lmfun <- function(x)   {
  lmfit <- lm(invsimpson~C, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

b.C.invsimps <- ddply(b.a.z.na, ~zone+horizon+om, lmfun)


C <- ggplot(data = b.C.invsimps, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Inv. Simpson ~ C")+ labs(title="Bacterial relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)


lmfun <- function(x)   {
  lmfit <- lm(invsimpson~CN, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

b.CN.invsimps <- ddply(b.a.z.na, ~zone+horizon+om, lmfun)

CN <- ggplot(data = b.CN.invsimps, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Inv. Simpson ~ CN")+ labs(title="Bacterial relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)

lmfun <- function(x)   {
  lmfit <- lm(invsimpson~N, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

b.N.invsimps <- ddply(b.a.z.na, ~zone+horizon+om, lmfun)


N <- ggplot(data = b.N.invsimps, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Inv. Simpson ~ N")+ labs(title="Bacterial relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)

lmfun <- function(x)   {
  lmfit <- lm(invsimpson~pH_H2O, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

b.pH.invsimps <- ddply(b.a.z.na, ~zone+horizon+om, lmfun)


pH <- ggplot(data = b.pH.invsimps, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Inv. Simpson ~ pH")+ labs(title="Bacterial relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)

#source multiplot
# devSVG(file="b.slopes.svg", width=8, height=10)
# multiplot(C, N, CN, pH, cols=2)
# dev.off()


devSVG(file="b.C.slopes.svg", height=10, width=8)
C
dev.off()
devSVG(file="b.CN.slopes.svg", height=10, width=8)
CN
dev.off()
devSVG(file="b.N.slopes.svg", height=10, width=8)
N
dev.off()
devSVG(file="b.pH.slopes.svg", height=10, width=8)
pH
dev.off()




#########

### Run everything with Richness



lmfun <- function(x)   {
  lmfit <- lm(sobs~C, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

f.C.sobs <- ddply(f.a.z.na, ~zone+horizon+om, lmfun)


C <- ggplot(data = f.C.sobs, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Richness ~ C")+ labs(title="Fungal relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)


lmfun <- function(x)   {
  lmfit <- lm(sobs~CN, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

f.CN.sobs <- ddply(f.a.z.na, ~zone+horizon+om, lmfun)

CN <- ggplot(data = f.CN.sobs, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Richness ~ CN")+ labs(title="Fungal relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)

lmfun <- function(x)   {
  lmfit <- lm(sobs~N, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

f.N.sobs <- ddply(f.a.z.na, ~zone+horizon+om, lmfun)


N <- ggplot(data = f.N.sobs, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Richness ~ N")+ labs(title="Fungal relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)

lmfun <- function(x)   {
  lmfit <- lm(sobs~pH_H2O, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

f.pH.sobs <- ddply(f.a.z.na, ~zone+horizon+om, lmfun)

pH <- ggplot(data = f.pH.sobs, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Richness ~ pH")+ labs(title="Fungal relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)
pH

devSVG(file="f.C.sobs.slopes.svg", height=10, width=8)
C
dev.off()
devSVG(file="f.CN.sobs.slopes.svg", height=10, width=8)
CN
dev.off()
devSVG(file="f.N.sobs.slopes.svg", height=10, width=8)
N
dev.off()
devSVG(file="f.pH.sobs.slopes.svg", height=10, width=8)
pH
dev.off()


#can't figure out how to specify x and y in the ddply command, have to replace in function each call

lmfun <- function(x)   {
  lmfit <- lm(sobs~C, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

b.C.sobs <- ddply(b.a.z.na, ~zone+horizon+om, lmfun)


C <- ggplot(data = b.C.sobs, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Richness ~ C")+ labs(title="Bacterial relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)


lmfun <- function(x)   {
  lmfit <- lm(sobs~CN, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

b.CN.sobs <- ddply(b.a.z.na, ~zone+horizon+om, lmfun)

CN <- ggplot(data = b.CN.sobs, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Richness ~ CN")+ labs(title="Bacterial relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)

lmfun <- function(x)   {
  lmfit <- lm(sobs~N, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

b.N.sobs <- ddply(b.a.z.na, ~zone+horizon+om, lmfun)


N <- ggplot(data = b.N.sobs, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Richness ~ N")+ labs(title="Bacterial relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)

lmfun <- function(x)   {
  lmfit <- lm(sobs~pH_H2O, x)
  output <- data.frame(t(coef(lmfit)), summary(lmfit)$r.square)
  names(output) <- c("intercept", "slope", "rsq")
  return(output)
}

b.pH.sobs <- ddply(b.a.z.na, ~zone+horizon+om, lmfun)


pH <- ggplot(data = b.pH.sobs, aes(x=om, y=slope, color=factor(zone)))+
  facet_grid(horizon ~ .)+
  scale_color_manual(values=zone.color)+
  xlab("OM")+ ylab("Slope Richness ~ pH")+ labs(title="Bacterial relativized diversity, standardized env")+
  geom_point(aes(alpha=rsq), color="black", size=4.5)+
  geom_point(size=3)

#source multiplot
# devSVG(file="b.slopes.svg", width=8, height=10)
# multiplot(C, N, CN, pH, cols=2)
# dev.off()


devSVG(file="b.C.sobs.slopes.svg", height=10, width=8)
C
dev.off()
devSVG(file="b.CN.sobs.slopes.svg", height=10, width=8)
CN
dev.off()
devSVG(file="b.N.sobs.slopes.svg", height=10, width=8)
N
dev.off()
devSVG(file="b.pH.sobs.slopes.svg", height=10, width=8)
pH
dev.off()