#import the whole datasets
f.otu<-t(read.table(file="../fung.all.mar14.hil.0.945.subsample.shared1.t.csv", header=T, row.names=1))
b.otu<-t(read.table(file="../bac_final.an.0.03.subsample.shared", header=T, row.names=1))

b.env <- read.table(file="../bac_env.txt", header=T, row.names=1)
f.env <- read.table(file="../Fung_env.csv", header=T, row.names=1)
library(labdsv)

#plot frequency abundance, quick and dirty can't control anything about the plots.  prettied up in inscape
#abuocc(f.otu)
#abuocc(b.otu)

#better frequency abundance plots
b.pres<-apply(b.otu>0,2,sum)
plot(sort(b.pres))
f.pres<-apply(f.otu>0,2,sum)
plot(sort(f.pres))
#otu mean where present
b.mean <- (apply(b.otu, 2, sum))/b.pres
f.mean <- (apply(f.otu, 2, sum))/f.pres

library(ggplot2)
b <- as.data.frame(cbind(b.mean, b.pres))
f <- as.data.frame(cbind(f.mean, f.pres))
# b1 <- ggplot(b, aes(x=b.pres, y=b.mean)) +
#   geom_point(alpha=1/10, size=3)+
#   scale_y_log10(name="Mean OTU abundance")+
#   scale_x_continuous(name="Sample occurance")+
#   theme_bw()+
#   ggtitle("Bacterial OTU abundance and occurance")
# 
# 
# f1 <- ggplot(f, aes(x=f.pres, y=f.mean)) +
#   geom_point(alpha=1/10, size=3)+
#   scale_y_log10(name="Mean OTU abundance")+
#   scale_x_continuous(name="Sample occurance")+
#   theme_bw()+
#   ggtitle("Fungal OTU abundance and occurance")

#add glm selected otus
b.sel<-t(read.table(file="../b.big.glm.otu.t", header=T, row.names=1))
f.sel <- t(read.table(file="../f.glm.otu.t", header=T, row.names=1))

b.pres<-apply(b.sel>0,2,sum)
f.pres<-apply(f.sel>0,2,sum)

#otu mean where present
b.mean <- (apply(b.sel, 2, sum))/b.pres
f.mean <- (apply(f.sel, 2, sum))/f.pres

b.s <- as.data.frame(cbind(b.mean, b.pres))
f.s <- as.data.frame(cbind(f.mean, f.pres))


# didn't work
# bs1 <- ggplot(b.s, aes(x=b.sel.pres, y=b.sel.mean))+
#  geom_point(size=4, color="red")+
#   scale_shape(solid=FALSE)
# b1+ggplot(b.s, aes(x=b.sel.pres, y=b.sel.mean))+
#  geom_point(aes(size=4, color="red")+
#   scale_shape(solid=FALSE)

f.s$group <- 2
f$group <- 1
f.s.d <- rbind(f,f.s)
f.s.d$group <- factor(f.s.d$group, labels=c("All OTUs", "OTUs responding to Harvesting Treatment"))
str(f.s.d)


library(RSvgDevice)#correctly saves text as text rather than converting to lines
devSVG(file="f.freq.abund.svg", width=10, height=8)
fs.plot <- ggplot(f.s.d, aes(x=f.pres, y=f.mean, alpha=group, color=group, size=group))+
  geom_point()+
  scale_size_manual(values=c(2,4))+
  scale_color_manual(values=c("black", "red"))+
  scale_y_log10(name="Mean OTU abundance")+
  scale_x_continuous(name="Sample occurance")+
  theme_bw()+
  ggtitle("Fungal OTU abundance and occurance")+
  theme(legend.position=c(.75,.9) ,legend.background=element_rect(colour="black", fill="white", size=0.3, linetype="solid"),
        legend.key=element_rect(colour="white"))+
  theme(panel.border=element_rect(fill=NA, color="black", size=.5))

fs.plot
dev.off()







# b.s$group <- 2
# b$group <- 1
# b.s.d <- rbind(b,b.s)
# b.s.d$group <- factor(b.s.d$group, labels=c("All OTUs", "OTUs responding to Harvesting Treatment"))
# #check that the groups copied over ok
# min(b.s.d$group)
# max(b.s.d$group)

devSVG(file="b.freq.abund.svg", width=10, height=8)
bs.plot <- ggplot(b.s.d, aes(x=b.pres, y=b.mean, alpha=group, color=group, size=group))+
  geom_point()+
  scale_color_manual(values=c("black", "red"))+
  scale_size_manual(values=c(2,4))+
  scale_y_log10(name="Mean OTU abundance")+
  scale_x_continuous(name="Sample occurance")+
  theme_bw()+
  ggtitle("Bacterial OTU abundance and occurance")+
  theme(legend.position=c(.75,.9) ,legend.background=element_rect(colour="black", fill="white", size=0.3, linetype="solid"),
        legend.key=element_rect(colour="white"))+
  theme(panel.border=element_rect(fill=NA, color="black", size=.5))

bs.plot
dev.off()

#doesn't preserve text in svg's
# ggsave(file="b.freq.abund.svg", plot=bs.plot, width=10, height=8)
# ggsave(file="f.freq.abund.svg", plot=fs.plot, width=10, height=8)
ggsave(file="b.freq.abund.pdf", plot=bs.plot, width=10, height=8)
ggsave(file="f.freq.abund.pdf", plot=fs.plot, width=10, height=8)



