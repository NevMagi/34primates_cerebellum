############################################
### Analysis carried out in response
### to Reviewer #1. Specifically, these
### analyses seek to explore whether the main
### result holds when we remove specimens that
### may have suspiciously low or high cerebellar
### volumes.
### 
### We scanned the literature specifically for
### cerebellar volumes, and found volumes for 
### most of the species in our analysis.
###
### Sources (see manuscript):
### Stephan et al. 1981
### Maseko et al. 2012
### Rilling et al. 1998
### Navarrete et al. 2018
### Smaers et al. 2011
### MacLeod et al. 2003

## NB outlier rhesus and crab-eating macaque
### specimens were removed before calculating
### species median. Species medians for these
### species were not suspicious versus the
### literature.

### NB2 species we did not find in the liter-
### ature were assumed to have normal cereb-
### ellar volumes.
############################################

library(tidyverse)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(ape)
library(phytools)
library(Rphylopars)
library(dispRity)
library(nlme)

### SET UP ###
getwd()
setwd("~/Documents/PhD/EvolutionPrimates/RevisionCode")

outputdir <- "5.Suspicious"
dir.create(outputdir)

suspicious <- read_csv2("./input/PhenotypesMedian_split_SUSPICIOUS.csv") %>% dplyr::rename(species = WilsonReederName) %>% dplyr::rename(Mass = `Body mass species mean`) %>% dplyr::select(MedianCerebellum, MedianCerebrum, MedianROC, MedianCrus, Suspicious) %>% as.data.frame()
suspicious$Suspicious <- as.factor(suspicious$Suspicious)
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
row.names(suspicious) <- tree[["tip.label"]] # Depends on alphabet sorting of your data.
tree<-treedata(tree,suspicious, sort=T,warnings=T)$phy # Match tree to the data


#One could manually assign groupings and filter based on that.
SUSPICIOUS <- suspicious %>% filter(Suspicious == "Yes")
tree_suspicious <- treedata(tree,SUSPICIOUS, sort=T,warnings=T)$phy # Match tree to the suspicious data

#One could manually assign groupings and filter based on that.
NOSUSPICIOUS <- suspicious %>% filter(Suspicious == "No")
tree_nosuspicious <- treedata(tree,NOSUSPICIOUS, sort=T,warnings=T)$phy # Match tree to the non-suspicious data.


#####################################
### CALCULATE PICs
#####################################
Csuspicious <- SUSPICIOUS$MedianCerebellum
Vsuspicious <- SUSPICIOUS$MedianCerebrum

Cnosuspicious <- NOSUSPICIOUS$MedianCerebellum
Vnosuspicious <- NOSUSPICIOUS$MedianCerebrum

# Compute contrasts
pic.C.suspicious <- pic(Csuspicious,tree_suspicious)
pic.V.suspicious <- pic(Vsuspicious, tree_suspicious)

pic.C.nosuspicious <- pic(Cnosuspicious,tree_nosuspicious)
pic.V.nosuspicious <- pic(Vnosuspicious, tree_nosuspicious)


#------------------------------------
# PIC-regressions suspicious
#------------------------------------
# Compute a linear regression of C on V
result <- lm(pic.C.suspicious ~ pic.V.suspicious)
sink("5.Suspicious/5.1.lm-C-V.suspicious.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("5.Suspicious/5.2.lm-C-V-mean-nointercept-plot.suspicious.pdf")
plot(pic.C.suspicious ~ pic.V.suspicious)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.
## They are the same; check! ##

## Check if PIC/PGLS models are the same.
bm.suspicious <-corBrownian(1,tree_suspicious)
pglsModel<-gls(MedianCerebellum~MedianCerebrum,data=SUSPICIOUS,correlation=bm.suspicious)
summary(pglsModel)
PICmodel <- lm(pic.C.suspicious ~pic.V.suspicious + 0)
sink("5.Suspicious/5.3.PGLS-PIC-C-V.suspicious.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


#------------------------------------
# PIC-regressions non-suspicious
#------------------------------------
# Compute a linear regression of C on V
result <- lm(pic.C.nosuspicious ~ pic.V.nosuspicious)
sink("5.Suspicious/5.4.lm-C-V.nosuspicious.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("5.Suspicious/5.5.lm-C-V-mean-nointercept-plot.nosuspicious.pdf")
plot(pic.C.nosuspicious ~ pic.V.nosuspicious)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.
## They are the same; check! ##

## Check if PIC/PGLS models are the same.
bm.nosuspicious <-corBrownian(1,tree_nosuspicious)
pglsModel<-gls(MedianCerebellum~MedianCerebrum,data=NOSUSPICIOUS,correlation=bm.nosuspicious)
summary(pglsModel)
PICmodel <- lm(pic.C.nosuspicious ~pic.V.nosuspicious + 0)
sink("5.Suspicious/5.6.PGLS-PIC-C-V.nosuspicious.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


###########################################################
### Use pANCOVA to test deviations in intercept or slope
### and with that deviation from the general allometric 
### relationship for cerebellar-to-cerebral scaling depending
### on wheter they are suspicious or not.
###########################################################


library(evomap)
library(geiger)

# This is an alternative way to subset the trees :). Very convenient
Suspicious <-getTips(tree_suspicious,findMRCA(tree_suspicious,tips = tree_suspicious[["tip.label"]]))
Nosuspicious <-getTips(tree_nosuspicious,findMRCA(tree_nosuspicious,tips = tree_nosuspicious[["tip.label"]]))


#------------------------------------
# pANCOVA - suspicious
#------------------------------------
#For differences in slope: 
grpS<-rep("A",length(rownames(suspicious))) 
grpS[Nosuspicious]<-"B" 
grpS<-as.factor(grpS) 
names(grpS)<-rownames(suspicious)

#For differences in intercept: 
grpI<-rep("A",length(rownames(suspicious))) 
grpI[Nosuspicious]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(suspicious)


## Data subsetting
## Cerebellum vs. Cerebrum
data<-primateData; tree<-primateTree
Y<-"Brain"; X<-"Body"
data<-data[,c(which(colnames(data)==Y),which(colnames(data)==X)),drop=F]
data<-na.omit(data); data<-log(data)
tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)   #match the data to the tree
colnames(data)<-c("Dependent","Independent")  
## Very hacky but bypasses weird error.
data <- data[1:34,]
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
tree <-treedata(tree,suspicious, sort=T,warnings=T)$phy
rownames(data) <- tree[["tip.label"]]

## Plot separate PGLS regressions for suspicious and non-suspicious data
## Cerebellum vs. Cerebrum

## Again, extremely hacky, but works.
data$Dependent <- suspicious$MedianCerebellum
data$Independent <- suspicious$MedianCerebrum

pdf("5.Suspicious/5.7.PGLS.suspicious.CC.pdf")
plot(data$Dependent~data$Independent,col="black",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=Suspicious,col="green",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=Nosuspicious,col="grey",lwd=5,cex=1,pch=19)
graphics.off()

## Model construction
Model.CC <-model.matrix(as.formula(Dependent~Independent),data)
Model.CC.s <-model.matrix(as.formula(Dependent~grpS:Independent),data) #slopes differ
Model.CC.i <-model.matrix(as.formula(Dependent~grpI + Independent),data) #intercepts differ
Model.CC.si <-model.matrix(as.formula(Dependent~grpI + grpS:Independent),data) #slopes & intercepts differ

sink("5.Suspicious/5.8.ModelTesting.SlopesIntercepts.CC.txt")
## Testing versus baseline model (common intercept and slope)
print("(1) Differences in slopes, holding intercept constant:")
gls.ancova(Dependent~Independent,vcv(tree),Model.CC,Model.CC.s)

print("(2) Differences in intercept, holding slopes constant:")
gls.ancova(Dependent~Independent,vcv(tree),Model.CC,Model.CC.i)

print("(3) Differences in slopes and differences in intercept:")
gls.ancova(Dependent~Independent,vcv(tree),Model.CC,Model.CC.si)
sink()



###########################################################
### Repeat for ansiform area regressions
###########################################################

suspiciousAA <- suspicious %>% drop_na()
suspiciousAA <- suspiciousAA[-c(9),] # We manually drop the crab-eating macaque, which was 
treeAA<-treedata(tree,suspiciousAA, sort=T,warnings=T)$phy # Match tree to the data


#One could manually assign groupings and filter based on that.
SUSPICIOUSAA <- suspiciousAA %>% filter(Suspicious == "Yes")
tree_suspiciousAA <- treedata(treeAA,SUSPICIOUSAA, sort=T,warnings=T)$phy # Match tree to the suspicious data

#One could manually assign groupings and filter based on that.
NOSUSPICIOUSAA <- suspiciousAA %>% filter(Suspicious == "No")
tree_nosuspiciousAA <- treedata(treeAA,NOSUSPICIOUSAA, sort=T,warnings=T)$phy # Match tree to the non-suspicious data.


#####################################
### CALCULATE PICs
#####################################
AAsuspicious <- SUSPICIOUSAA$MedianCrus
ROCsuspicious <- SUSPICIOUSAA$MedianROC
Vsuspicious <- SUSPICIOUSAA$MedianCerebrum

AAnosuspicious <- NOSUSPICIOUSAA$MedianCrus
ROCnosuspicious <- NOSUSPICIOUSAA$MedianROC
Vnosuspicious <- NOSUSPICIOUSAA$MedianCerebrum

# Compute contrasts
pic.AA.suspiciousAA <- pic(AAsuspicious,tree_suspiciousAA)
pic.ROC.suspiciousAA <- pic(ROCsuspicious, tree_suspiciousAA)
pic.V.suspiciousAA <- pic(VsuspiciousAA, tree_suspiciousAA)

pic.AA.nosuspiciousAA <- pic(AAnosuspicious, tree_nosuspiciousAA)
pic.ROC.nosuspiciousAA <- pic(ROCnosuspicious, tree_nosuspiciousAA)
pic.V.nosuspiciousAA <- pic(VnosuspiciousAA, tree_nosuspiciousAA)

#------------------------------------
# PIC-regressions suspicious AA-ROC
#------------------------------------
# Compute a linear regression of AA on ROC
result <- lm(pic.AA.suspiciousAA ~ pic.ROC.suspiciousAA)
sink("5.Suspicious/5.9.lm-AA-ROC.suspicious.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("5.Suspicious/5.10.lm-AA-ROC-mean-nointercept-plot.suspicious.pdf")
plot(pic.AA.suspiciousAA ~ pic.ROC.suspiciousAA)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.
## They are the same; check! ##

## Check if PIC/PGLS models are the same.
bm.suspiciousAA <-corBrownian(1,tree_suspiciousAA)
pglsModel<-gls(MedianCrus~MedianROC,data=SUSPICIOUSAA,correlation=bm.suspiciousAA)
summary(pglsModel)
PICmodel <- lm(pic.AA.suspiciousAA ~pic.ROC.suspiciousAA + 0)
sink("5.Suspicious/5.11.PGLS-PIC-AA-ROC.suspicious.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


#------------------------------------
# PIC-regressions non-suspicious AA-ROC
#------------------------------------
# Compute a linear regression of AA on ROC
result <- lm(pic.AA.nosuspiciousAA ~ pic.ROC.nosuspiciousAA)
sink("5.Suspicious/5.12.lm-AA-ROC.nosuspicious.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("5.Suspicious/5.13.lm-AA-ROC-mean-nointercept-plot.nosuspicious.pdf")
plot(pic.AA.nosuspiciousAA ~ pic.ROC.nosuspiciousAA)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.
## They are the same; check! ##

## Check if PIC/PGLS models are the same.
bm.nosuspiciousAA <-corBrownian(1,tree_nosuspiciousAA)
pglsModel<-gls(MedianCrus~MedianROC,data=NOSUSPICIOUSAA,correlation=bm.nosuspiciousAA)
summary(pglsModel)
PICmodel <- lm(pic.AA.nosuspiciousAA ~pic.ROC.nosuspiciousAA + 0)
sink("5.Suspicious/5.14.PGLS-PIC-AA-ROC.nosuspicious.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

#------------------------------------
# PIC-regressions non-suspicious AA-Cerebllum
#------------------------------------
# Compute a linear regression of AA on Cerebellum
result <- lm(pic.AA.nosuspiciousAA ~ pic.C.nosuspiciousAA)
sink("5.Suspicious/5.12a.lm-AA-C.nosuspicious.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("5.Suspicious/5.13a.lm-AA-C-mean-nointercept-plot.nosuspicious.pdf")
plot(pic.AA.nosuspiciousAA ~ pic.C.nosuspiciousAA)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.
## They are the same; check! ##

## Check if PIC/PGLS models are the same.
bm.nosuspiciousAA <-corBrownian(1,tree_nosuspiciousAA)
pglsModel<-gls(MedianCrus~MedianCerebellum,data=NOSUSPICIOUSAA,correlation=bm.nosuspiciousAA)
summary(pglsModel)
PICmodel <- lm(pic.AA.nosuspiciousAA ~pic.C.nosuspiciousAA + 0)
sink("5.Suspicious/5.14a.PGLS-PIC-AA-C.nosuspicious.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

#------------------------------------
# PIC-regressions suspicious AA-V
#------------------------------------
# Compute a linear regression of AA on V
result <- lm(pic.AA.suspiciousAA ~ pic.V.suspiciousAA)
sink("5.Suspicious/5.15.lm-AA-V.suspicious.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("5.Suspicious/5.16.lm-AA-V-mean-nointercept-plot.suspicious.pdf")
plot(pic.AA.suspiciousAA ~ pic.V.suspiciousAA)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.
## They are the same; check! ##

## Check if PIC/PGLS models are the same.
bm.suspiciousAA <-corBrownian(1,tree_suspiciousAA)
pglsModel<-gls(MedianCrus~MedianCerebrum,data=SUSPICIOUSAA,correlation=bm.suspiciousAA)
summary(pglsModel)
PICmodel <- lm(pic.AA.suspiciousAA ~pic.V.suspiciousAA + 0)
sink("5.Suspicious/5.17.PGLS-PIC-AA-V.suspicious.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

#------------------------------------
# PIC-regressions non-suspicious AA-V
#------------------------------------
# Compute a linear regression of AA on V
result <- lm(pic.AA.nosuspiciousAA ~ pic.V.nosuspiciousAA)
sink("5.Suspicious/5.18.lm-AA-V.nosuspicious.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("5.Suspicious/5.19.lm-AA-V-mean-nointercept-plot.nosuspicious.pdf")
plot(pic.AA.nosuspiciousAA ~ pic.V.nosuspiciousAA)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.
## They are the same; check! ##

## Check if PIC/PGLS models are the same.
bm.nosuspiciousAA <-corBrownian(1,tree_nosuspiciousAA)
pglsModel<-gls(MedianCrus~MedianCerebrum,data=NOSUSPICIOUSAA,correlation=bm.nosuspiciousAA)
summary(pglsModel)
PICmodel <- lm(pic.AA.nosuspiciousAA ~pic.V.nosuspiciousAA + 0)
sink("5.Suspicious/5.20.PGLS-PIC-AA-V.nosuspicious.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


###-----------------------------------------------###
# Plot figures supplement - non-suspicious data
####----------------------------------------------###

### -------------------------- ###
### PGLS Plotting with Caper
### -------------------------- ###   
library(caper)
library(evomap)

# Model 1: Cerebellar volume regressed on cerebral volume.
tiff("./5.Suspicious/5.21.BM-PGLS-cerebellumCerebrum-nonsuspicious-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod1 <- gls(MedianCerebellum ~ MedianCerebrum, data=NOSUSPICIOUS, correlation=bm.nosuspicious)
plot(mod1)
plot(NOSUSPICIOUS$MedianCerebellum ~ NOSUSPICIOUS$MedianCerebrum, xlab = "Cerebral volume (log10-transformed)", ylab = "Cerebellar volume (log10-transformed)")
abline(mod1)
abline(a=-0.5656694, b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod1.ci <-gls.ci(NOSUSPICIOUS$MedianCerebellum, NOSUSPICIOUS$MedianCerebrum,vcv(tree_nosuspicious))
lines(pGLS.mod1.ci$CI.plot$X,pGLS.mod1.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod1.ci$CI.plot$X,pGLS.mod1.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod1.pi <-gls.pi(NOSUSPICIOUS$MedianCerebellum, NOSUSPICIOUS$MedianCerebrum,vcv(tree_nosuspicious), 1)
lines(pGLS.mod1.pi$PI.plot$X,pGLS.mod1.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod1.pi$PI.plot$X,pGLS.mod1.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()

# Model 2: Ansiform area volume regressed on cerebellar volume.
tiff("./5.Suspicious/5.22.BM-PGLS-ansiformareaCerebellum-nonsuspicious-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod2 <- gls(MedianCrus ~ MedianCerebellum, data=NOSUSPICIOUSAA, correlation=bm.nosuspiciousAA)
plot(mod2)
plot(NOSUSPICIOUSAA$MedianCrus ~ NOSUSPICIOUSAA$MedianCerebellum, xlab = "Cerebellar volume (log10-transformed)", ylab = "Ansiform area volume (log10-transformed)")
abline(mod2)
abline(a=-2.090477, b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod2.ci <-gls.ci(NOSUSPICIOUSAA$MedianCrus, NOSUSPICIOUSAA$MedianCerebellum,vcv(tree_nosuspiciousAA))
lines(pGLS.mod2.ci$CI.plot$X,pGLS.mod2.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod2.ci$CI.plot$X,pGLS.mod2.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod2.pi <-gls.pi(NOSUSPICIOUSAA$MedianCrus, NOSUSPICIOUSAA$MedianCerebellum,vcv(tree_nosuspiciousAA), 1)
lines(pGLS.mod2.pi$PI.plot$X,pGLS.mod2.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod2.pi$PI.plot$X,pGLS.mod2.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()

# Model 3-better: Ansiform area volume regressed on rest of cerebellar volume.
tiff("./5.Suspicious/5.22a.BM-PGLS-ansiformareaROC-nonsuspicious-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod3 <- gls(MedianCrus ~ MedianROC, data=NOSUSPICIOUSAA, correlation=bm.nosuspiciousAA)
plot(mod3)
plot(NOSUSPICIOUSAA$MedianCrus ~ NOSUSPICIOUSAA$MedianROC, xlab = "Rest of cerebellar volume (log10-transformed)", ylab = "Ansiform area volume (log10-transformed)")
abline(mod3)
abline(a=-1.915490, b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod3.ci <-gls.ci(NOSUSPICIOUSAA$MedianCrus, NOSUSPICIOUSAA$MedianCerebellum,vcv(tree_nosuspiciousAA))
lines(pGLS.mod3.ci$CI.plot$X,pGLS.mod3.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod3.ci$CI.plot$X,pGLS.mod3.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod3.pi <-gls.pi(NOSUSPICIOUSAA$MedianCrus, NOSUSPICIOUSAA$MedianCerebellum,vcv(tree_nosuspiciousAA), 1)
lines(pGLS.mod3.pi$PI.plot$X,pGLS.mod3.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod3.pi$PI.plot$X,pGLS.mod3.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()

# Model 4: Ansiform area volume regressed on cerebral volume.
tiff("./5.Suspicious/5.23.BM-PGLS-ansiformareaCerebrum-nonsuspicious-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod4 <- gls(MedianCrus ~ MedianCerebrum, data=NOSUSPICIOUSAA, correlation=bm.nosuspiciousAA)
plot(mod4)
plot(NOSUSPICIOUSAA$MedianCrus ~ NOSUSPICIOUSAA$MedianCerebrum,  xlab = "Cerebral volume (log10-transformed)", ylab = "Ansiform area volume (log10-transformed)")
abline(mod4)
abline(a=-2.740146, b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod4.ci <-gls.ci(NOSUSPICIOUSAA$MedianCrus, NOSUSPICIOUSAA$MedianCerebrum,vcv(tree_nosuspiciousAA))
lines(pGLS.mod4.ci$CI.plot$X,pGLS.mod4.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod4.ci$CI.plot$X,pGLS.mod4.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod4.pi <-gls.pi(NOSUSPICIOUSAA$MedianCrus, NOSUSPICIOUSAA$MedianCerebrum,vcv(tree_nosuspiciousAA), 1)
lines(pGLS.mod4.pi$PI.plot$X,pGLS.mod4.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod4.pi$PI.plot$X,pGLS.mod4.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()


### -------------------------- ###
### Calculate R2s for PGLS models
### -------------------------- ###   

### LEFT OF HERE ###

sink("5.phylogenetic_regressions/5.14.r2_phylogenetic_regressions.txt")
print("Cerebellar volume regressed on cerebral volume")
mod1.r2 <- R2(mod = mod1, phy = tree, sigma2_d = "s2w") #check out the meanings of the different R2s
mod1.r2
print("Ansiform area volume regressed on cerebellar volume")
mod3.r2 <- R2(mod = mod2, phy = tree13, sigma2_d = "s2w")
mod3.r2
print("Ansiform area volume regressed on cerebellar volume")
mod3a.r2 <- R2(mod = mod3a, phy = tree13, sigma2_d = "s2w")
mod3a.r2
print("Ansiform area volume regressed on cerebral volume")
mod4.r2 <- R2(mod = mod3, phy = tree13, sigma2_d = "s2w")
mod4.r2
sink()










#------------------------------------
# END OF SCRIPT
#------------------------------------
