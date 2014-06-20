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

#cal mean ref values
means <- aggregate(b.a[b.env$om==0,], by=list(b.env$site.hor[b.env$om==0]), FUN=mean)
means.dtf <- data.frame(means, row.names=1)
means.ext <- means.dtf[b.env$site.hor,]
b.a.rel <- b.a/means.ext

means <- aggregate(f.a[f.env$om==0,], by=list(f.env$site.hor[f.env$om==0]), FUN=mean)
means.dtf <- data.frame(means, row.names=1)
means.ext <- means.dtf[f.env$site.hor,]
f.a.rel <- f.a/means.ext

b.a.rel.env <- cbind(b.a.rel, b.env$zone, b.env$om, b.env$horizon)
f.a.rel.env <- cbind(f.a.rel, f.env$zone, f.env$om, f.env$horizon)
names(b.a.rel.env)[7] <- ("zone")
names(b.a.rel.env)[8] <- ("om")
names(b.a.rel.env)[9] <- ("horizon")
names(f.a.rel.env)[7] <- ("zone")
names(f.a.rel.env)[8] <- ("om")
names(f.a.rel.env)[9] <- ("horizon")


b.a.rel.env.wfire <- b.a.rel.env
f.a.rel.env.wfire <- f.a.rel.env

b.a.rel.env  <- subset(b.a.rel.env, b.env$om<98)

f.a.rel.env  <- subset(f.a.rel.env, f.env$om<98)



#quick plots to explore shifts from ref
library(ggplot2)
#can't see anything with just points
# cov <- ggplot(data=b.a.rel.env, aes(x=factor(b.env$om), y=coverage, color=b.env$zone))+
#   geom_point()
cov <- ggplot(data=subset(b.a.rel.env, b.a.rel.env$om<98), aes(x=factor(om), y=coverage))+
  geom_boxplot()+
  ggtitle("Bacterial Relative to Reference")
cov

invsimp <- ggplot(data=subset(b.a.rel.env, b.a.rel.env$om<98), aes(x=factor(om), y=invsimpson))+
  geom_boxplot()+
  ggtitle("Bacterial Relative to Reference")
invsimp

even <- ggplot(data=subset(b.a.rel.env, b.a.rel.env$om<98), aes(x=factor(om), y=shannoneven))+
  geom_boxplot()+
  ggtitle("Bacterial Relative to Reference")
even


cov <- ggplot(data=subset(f.a.rel.env, f.a.rel.env$om<98), aes(x=factor(om), y=coverage))+
  geom_boxplot()+
  ggtitle("Fungal Relative to Reference")
cov

invsimp <- ggplot(data=subset(f.a.rel.env, f.a.rel.env$om<98), aes(x=factor(om), y=invsimpson))+
  geom_boxplot()+
  ggtitle("Fungal Relative to Reference")
invsimp

even <- ggplot(data=subset(f.a.rel.env, f.a.rel.env$om<98), aes(x=factor(om), y=shannoneven))+
  geom_boxplot()+
  ggtitle("Fungal Relative to Reference")
even






