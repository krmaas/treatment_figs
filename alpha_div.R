library(ggplot2)
library(RSvgDevice)




#read in data
b.a <- read.table(file="../bac.alpha.csv", header=T, row.names=1)
f.a <- read.table(file="../fung.alpha.cvs", header=T, row.names=1)

b.env <- read.table(file="../bac_env.txt", header=T, row.names=1)
f.env <- read.table(file="../Fung_env.csv", header=T, row.names=1)

#sort both to make sure they match
b.a <- b.a[order(row.names(b.a)),]
f.a <- f.a[order(row.names(f.a)),]
b.env <- b.env[order(row.names(b.env)),]
f.env <- f.env[order(row.names(f.env)),]

#will relativize to Ref plots, need to add site/hor variable
b.env$site.hor <- paste(b.env$site, b.env$horizon, sep=".")
f.env$site.hor <- paste(f.env$site, f.env$horizon, sep=".")

b.a.env <- cbind(b.a, b.env$zone, b.env$horizon, b.env$om)
f.a.env <- cbind(f.a, f.env$zone, f.env$horizon, f.env$om)
names(b.a.env)[7] <- ("zone")
names(b.a.env)[8] <- ("horizon")
names(b.a.env)[9] <- ("om")
names(f.a.env)[7] <- ("zone")
names(f.a.env)[8] <- ("horizon")
names(f.a.env)[9] <- ("om")

b.a.env <- subset(b.a.env, om<"98")
f.a.env <- subset(f.a.env, om<"98")


#cal mean ref values
means <- aggregate(b.a[b.env$om==0,], by=list(b.env$site.hor[b.env$om==0]), FUN=mean)
means.dtf <- data.frame(means, row.names=1)
means.ext <- means.dtf[b.env$site.hor,]
b.a.rel.sub <- (b.a-means.ext)#/means.ext

means <- aggregate(f.a[f.env$om==0,], by=list(f.env$site.hor[f.env$om==0]), FUN=mean)
means.dtf <- data.frame(means, row.names=1)
means.ext <- means.dtf[f.env$site.hor,]
f.a.rel.sub <- (f.a-means.ext)#/means.ext

b.a.sub.env <- cbind(b.a.rel.sub, b.env$zone, b.env$om, b.env$horizon)
f.a.sub.env <- cbind(f.a.rel.sub, f.env$zone, f.env$om, f.env$horizon)
names(b.a.sub.env)[7] <- ("zone")
names(b.a.sub.env)[8] <- ("om")
names(b.a.sub.env)[9] <- ("horizon")
names(f.a.sub.env)[7] <- ("zone")
names(f.a.sub.env)[8] <- ("om")
names(f.a.sub.env)[9] <- ("horizon")


# b.a.rel.env.wfire <- b.a.rel.env
# f.a.rel.env.wfire <- f.a.rel.env
# b.a.rel.env1 <- b.a.rel.env
# f.a.rel.env1 <- f.a.rel.env
b.a.rel.env  <- subset(b.a.sub.env, b.env$om<98)

f.a.rel.env  <- subset(f.a.sub.env, f.env$om<98)



#quick plots to explore shifts from ref


#can't see anything with just points
# cov <- ggplot(data=b.a.rel.env, aes(x=factor(b.env$om), y=coverage, color=b.env$zone))+
#   geom_point()
# cov <- ggplot(data=subset(b.a.rel.env, b.a.rel.env$om<98), aes(x=factor(om), y=coverage, fill=factor(horizon)))+
#   geom_hline(yintercept=0, aes(col="grey", linetype=b))+
#   geom_boxplot()+
#   scale_fill_manual(values=c("grey","white"))+
#   theme_bw()+
#   ggtitle("Bacterial Relative to Reference")
# cov

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

ggsave(file="b.invsimp.pdf", plot=b.i, width=10, height=8)


# even <- ggplot(data=subset(b.a.rel.env, b.a.rel.env$om<98), aes(x=factor(om), y=shannoneven, fill=factor(horizon)))+
#   geom_hline(yintercept=0, aes(col="grey", linetype=b))+
#   geom_boxplot()+
#   scale_fill_manual(values=c("grey","white"))+
#   theme_bw()+
#   ggtitle("Bacterial Relative to Reference")
# even
# 
# sobs <- ggplot(data=subset(b.a.rel.env, b.a.rel.env$om<98), aes(x=factor(om), y=sobs, fill=factor(horizon)))+
#   geom_hline(yintercept=0, aes(col="grey", linetype=b))+
#   geom_boxplot()+
#   scale_fill_manual(values=c("grey","white"))+
#   theme_bw()+
#   ggtitle("Bacterial Relative to Reference")
# sobs
# 
# cov <- ggplot(data=subset(f.a.rel.env, f.a.rel.env$om<98), aes(x=factor(om), y=coverage, fill=factor(horizon)))+
#   geom_hline(yintercept=0, aes(col="grey", linetype=b))+
#   geom_boxplot()+
#   scale_fill_manual(values=c("grey","white"))+
#   theme_bw()+
#   ggtitle("Fungal Relative to Reference")
# cov

devSVG(file="f.invsimp.svg", width=10, height=8)

invsimp <- ggplot(data=subset(f.a.rel.env, f.a.rel.env$om<98), aes(x=factor(om), y=invsimpson, fill=factor(horizon)))+
  geom_hline(yintercept=0, aes(col="grey", linetype=b))+
  geom_boxplot(outlier.size=0)+
  scale_fill_manual(values=c("grey","white"))+
  theme_bw()+
  ggtitle("Fungal Relative to Reference")
invsimp
ylim1 <- boxplot.stats(f.a.rel.env$invsimpson)$stats[c(1,5)]
f.i <- invsimp + coord_cartesian(ylim=ylim1*1.5)
f.i
dev.off()

ggsave(file="f.invsimp.pdf", plot=f.i, width=10, height=8)
# 
# even <- ggplot(data=subset(f.a.rel.env, f.a.rel.env$om<98), aes(x=factor(om), y=shannoneven, fill=factor(horizon)))+
#   geom_hline(yintercept=0, aes(col="grey", linetype=b))+
#   geom_boxplot()+
#   scale_fill_manual(values=c("grey","white"))+
#   #   scale_y_continuous(limits=c(-1,1))+
#   theme_bw()+
#   ggtitle("Fungal Relative to Reference")
# even
# 
# sobs <- ggplot(data=subset(f.a.rel.env, f.a.rel.env$om<98), aes(x=factor(om), y=sobs, fill=factor(horizon)))+
#   geom_hline(yintercept=0, aes(col="grey", linetype=b))+
#   geom_boxplot()+
#   scale_fill_manual(values=c("grey","white"))+
#   #   scale_y_continuous(limits=c(-1,1))+
#   theme_bw()+
#   ggtitle("Fungal Relative to Reference")
# sobs
# sobs <- ggplot(data=f.a.rel.env, aes(x=factor(horizon), y=sobs))+
#   geom_boxplot()+
#   ggtitle("Fungal Relative to Reference")
# even
# 
# data <- f.a.rel.env[f.a.rel.env$om!="99",]
# 

#not seeing tons of trends with the box plots
b.cov <- aov(coverage~factor(zone)+factor(om), data=b.a.rel.env)
summary(b.cov)
# TukeyHSD(b.cov)

b.invs <- aov(invsimpson~factor(zone)+factor(om), data=b.a.rel.env)
summary(b.invs)
# TukeyHSD(b.invs)

b.even <- aov(shannoneven~factor(zone)+factor(om), data=b.a.rel.env)
summary(b.even)
# TukeyHSD(b.even)

b.sobs <- aov(sobs~factor(zone)+factor(om), data=b.a.rel.env)
summary(b.sobs)
# TukeyHSD(b.sobs)


f.cov <- aov(coverage~factor(zone)+factor(om), data=f.a.rel.env)
summary(f.cov)
# TukeyHSD(f.cov)

f.invs <- aov(invsimpson~factor(zone)+factor(om), data=f.a.rel.env)
summary(f.invs)
# TukeyHSD(f.invs)

f.even <- aov(shannoneven~factor(zone)+factor(om), data=f.a.rel.env)
summary(f.even)
# TukeyHSD(f.even)

f.sobs <- aov(sobs~factor(zone)+factor(om), data=f.a.rel.env)
summary(f.sobs)
# TukeyHSD(f.sobs)


#alpha's vs zone and hor
b.cov.hor <- aov(coverage~factor(horizon), data=b.a.rel.env)
summary(b.cov.hor)
# TukeyHSD(b.cov.hor)

b.invs.hor <- aov(invsimpson~factor(horizon), data=b.a.rel.env)
summary(b.invs.hor)
TukeyHSD(b.invs.hor)

b.even.hor <- aov(shannoneven~factor(horizon), data=b.a.rel.env)
summary(b.even.hor)
TukeyHSD(b.even.hor)

b.sobs.hor <- aov(sobs~factor(horizon), data=b.a.rel.env)
summary(b.sobs.hor)
TukeyHSD(b.sobs.hor)


f.cov.hor <- aov(coverage~factor(horizon), data=f.a.rel.env)
summary(f.cov.hor)
TukeyHSD(f.cov.hor)

f.invs.hor <- aov(invsimpson~factor(horizon), data=f.a.rel.env)
summary(f.invs.hor)
TukeyHSD(f.invs.hor)

f.even.hor <- aov(shannoneven~factor(horizon), data=f.a.rel.env)
summary(f.even.hor)
TukeyHSD(f.even.hor)

f.sobs.hor <- aov(sobs~factor(horizon), data=f.a.rel.env)
summary(f.sobs.hor)
TukeyHSD(f.sobs.hor)


# #run plots one at a time
# plot(b.cov)
# plot(b.invs)
# plot(b.even)
# plot(b.sobs)


######really basic stats on the data
#create function that calculates standard error of the mean
sem <- function(x)
   {
       sqrt(var(x)/length(x))
    }

mean(f.a$coverage)
sem(f.a$coverage)
mean(b.a$coverage)
sem(b.a$coverage)
mean(f.a$shannoneven)
sem(f.a$shannoneven)
mean(b.a$shannoneven)
sem(b.a$shannoneven)
mean(f.a$invsimpson)
sem(f.a$invsimpson)
mean(b.a$invsimpson)
sem(b.a$invsimpson)

sem(f.a$coverage)/mean(f.a$coverage)
sem(b.a$coverage)/mean(b.a$coverage)
sem(f.a$shannoneven)/mean(f.a$shannoneven)
sem(b.a$shannoneven)/mean(b.a$shannoneven)
sem(f.a$invsimpson)/mean(f.a$invsimpson)
sem(b.a$invsimpson)/mean(b.a$invsimpson)


