setwd("/media/data/final_bac/full")
library(vegan)
library(ecodist)
library(ade4)

b_bc_d <- "bac_final.an.braycurtis.0.09.square.ave.dist"
b_bc09_ave <-data.matrix(read.table(b_bc_d, row.names=1, header=T))
b_bc09ave <- lower(b_bc09_ave)
scree09 <- nmds(b_bc09ave, mindim=1, maxdim=6, nits=100)
stress09 <- scree09$stress
plot(stress09)

#read in thetayc
b_bc_d <- "bac_final.an.thetayc.0.09.square.ave.nomouse.dist.csv"
b_th09_ave <-data.matrix(read.table(b_bc_d, row.names=1, header=T))
#determine dim for nms
b_th09ave <- lower(b_th09_ave)
scree09th <- nmds(b_th09ave, mindim=1, maxdim=6, nits=10)
stress09th <- scree09th$stress
plot(stress09th)


#nms 3 axes
nms09 <-metaMDS(b_th09_ave, k=3, trymin=50, trymax=250, wascores=FALSE)


#read lt distance 2
b_bc_d <- "bac_final.an.03.filter2.braycurtis.0.03.lt.ave.dist"
b_bc2_mat <-data.matrix(read.table(b_bc_d, fill=T, row.names=1, header=T))
View(b_bc2_mat)

#read lt distance 10
b_bc_d <- "bac_final.an.03.filter10.braycurtis.0.03.lt.ave.dist"
b_bc10_mat <-data.matrix(read.table(b_bc_d, fill=T, row.names=1, header=T))

#read lt distance 100
)b_bc_d <- "bac_final.an.03.filter100.braycurtis.0.03.lt.ave.dist"
b_bc100_mat <-data.matrix(read.table(b_bc_d, fill=T, row.names=1, header=T))

#read lt distance 1000
b_bc_d <- "bac_final.an.03.filter1000.braycurtis.0.03.lt.ave.dist"
b_bc1000_mat <-data.matrix(read.table(b_bc_d, fill=T, row.names=1, header=T))


#read 3% bac bc for comparing effect of removing rares
b_bc_d <- "bac_final.an.braycurtis.0.03.square.nomouse.csv"
b_bc03 <-data.matrix(read.table(b_bc_d, row.names=1, header=T))

b_bc_d <- "bac_final.an.braycurtis.0.05.square.ave.dist"
b_bc05_ave <-data.matrix(read.table(b_bc_d, row.names=1, header=T))

b_bc_d <- "bac_final.an.braycurtis.0.09.square.ave.dist"
b_bc09_ave <-data.matrix(read.table(b_bc_d, row.names=1, header=T))

#mantel on pairs of matrix
mantel.rtest(b_bc2_mat, b_bc10_mat, nrepet=9999)


#determine dim for nms
b_bc03ave <- lower(b_bc03_ave)
scree03 <- nmds(b_bc03ave, mindim=1, maxdim=6, nits=100)
stress03 <- scree03$stress
plot(stress03)

#read in thetayc
b_bc_d <- "bac_final.an.thetayc.0.03.square.ave.nomouse.dist.csv"
b_th03_ave <-data.matrix(read.table(b_bc_d, row.names=1, header=T))
#determine dim for nms
b_th03ave <- lower(b_th03_ave)
scree03th <- nmds(b_th03ave, mindim=1, maxdim=6, nits=10)
stress03th <- scree03th$stress
plot(stress03th)

#nms 3 axes
nms03 <-metaMDS(b_th03_ave, k=3, trymin=50, trymax=250, wascores=FALSE)

#read in env
bac_env <- read.table("bac_seq_names1.csv", row.names=1, header=T)
bac_env_noplot <- read.table("bac_seq_names1_noplot.csv", row.names=1, header=T)

envfit(nms03, bac_env, perm=999, choices=1:3)

b_th03ave <- as.dist(b_th03_ave)
disp<-betadisper(b_th03ave, bac_env$Zone)

 zone_om<-paste(bac_env$Zone, bac_env$om, sep="_")
fix(zone_om)
disp_zoneom <-betadisper(b_th03ave, zone_om, choices=1:3)
Error in betadisper(b_th03ave, zone_om, choices = 1:3) : 
  unused argument(s) (choices = 1:3)
disp_zoneom <-betadisper(b_th03ave, zone_om)
plot(disp_zoneom)
 boxplot(disp_zoneom)
 om_hor<-paste(bac_env$om, bac_env$Horizon, sep="_")
 disp_omhor<-betadisper(b_th03ave, om_hor)
 boxplot(disp_omhor)
 View(bac_env)
 om_hor<-paste(bac_env$om, bac_env$horizon, sep="_")
 disp_omhor<-betadisper(b_th03ave, om_hor)
 boxplot(disp_omhor)
 zone_hor_om<-paste(bac_env$Zone, bac_env$horizon, bac_env$om)
 disp_zonehorom<-betadisper(b_th03ave, zone_hor_om)
 boxplot(disp_zonehorom, las="2")
#########Fung

#read in thetayc
f_th03_d <- "fung_all.unique.chop.precluster.pick.an.thetayc.0.03.lt.ave.dist"
f_th03_ave <-(read.table(f_th03_d, row.names=1, header=T, fill=TRUE))
#determine dim for nms
f_th03_ave <-(read.table(f_th03_d, row.names=1, header=T, fill=TRUE))
f_th03_ave<-lower(f_th03_ave)
f_scree<-nmds(f_th03_ave, mindim=1, maxdim=6, nits=100)
fstress03th <- f_scree$stress
plot(fstress03th)


bSBS_nms<-metaMDS(bcSBSb03f100_otu, distance="bray", mindim=1, maxdim=6, nits=1000)
bIDF_nms<-nmds(bcIDFb03f100_otu, distance="bray", mindim=1, maxdim=6, nits=1000)
oBS_nms<-nmds(oBSb03f100_otu, distance="bray", mindim=1, maxdim=6, nits=1000)
oJP_nms<-nmds(oJPb03f100_otu, distance="bray", mindim=1, maxdim=6, nits=1000)
tx_nms<-nmds(txb03f100_otu, distance="bray", mindim=1, maxdim=6, nits=1000)
MD_nms<-nmds(MDb03f100_otu, distance="bray", mindim=1, maxdim=6, nits=1000)
