
############################################
### Analysis carried out in response
### to Reviewer #3. Specifically, these
### analyses seek to explore strepsirrhine
### versus haplorrhine scaling and how these
### groups interact with the main PGLS by
### dummy coding these groups.
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
library(geiger)

### SET UP ###
getwd()
setwd("~/Documents/PhD/EvolutionPrimates/RevisionCode")

outputdir <- "1.strepsirrhine"
dir.create(outputdir)

strepsirr <- read_csv2("./input/PhenotypesMedian_split_strepsirr.csv") %>% dplyr::rename(species = WilsonReederName) %>% dplyr::rename(Mass = `Body mass species mean` ) %>% dplyr::select(MedianCerebellum, MedianCerebrum, MedianCrus, MedianROC, MedianCerebroCerebellar, MedianCerebellarCrus, Mass, Strepsirrhine, BrainVol) %>% as.data.frame()
strepsirr$Strepsirrhine <- as.factor(strepsirr$Strepsirrhine)
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
row.names(strepsirr) <- tree[["tip.label"]] # Depends on alphabet sorting of your data.
strepsirr <- strepsirr %>% drop_na(Mass)
tree<-treedata(tree,strepsirr, sort=T,warnings=T)$phy # Match tree to the data

#One could manually assign groupings and filter based on that.
STREP <-strepsirr %>% filter(Strepsirrhine== "Yes")
tree_strep <- treedata(tree,STREP, sort=T,warnings=T)$phy # Match tree to the strepsirrhine data

# The haplorhine data contains a species that had no body mass data, that we will remove from the tree as well.
HAPL <-strepsirr %>% filter(Strepsirrhine== "No") %>% drop_na(Mass)
tree_hapl<-treedata(tree,HAPL, sort=T,warnings=T)$phy # Match tree to the data

## Remove the macaca fasciularis and mulatta rows for the  ansiform analysis
HAPLAA <- HAPL[-c(16:17),] %>% drop_na(MedianCrus) # manually removed the macaques.
tree_haplaa <- treedata(tree, HAPLAA, sort=T, warnings=T)$phy

# This is an alternative way to subset the trees :). Very convenient
Strepsirrhines <-getTips(tree,findMRCA(tree,c("Galago_demidoff","Microcebus_murinus")))
Haplorhines <-getTips(tree,findMRCA(tree,c("Saimiri_sciureus","Cercopithecus_cephus_cephus")))


########################################################
### REGULAR REGRESSION - body weight Isler et al. 2008
########################################################

## Cerebellum

strepsirr %>%
  group_by(Strepsirrhine) %>%
  summarize(cor=cor(MedianCerebellum, Mass))

pdf("1.strepsirrhine/1.0a.cerebellar.pdf")
strep_cerebellar <- ggplot(strepsirr, aes(x = Mass, y = MedianCerebellum, color = Strepsirrhine, grp.label = Strepsirrhine)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("grp", "eq", "R2"))) +
  geom_point()
strep_cerebellar 
graphics.off()

pdf("1.strepsirrhine/1.0b.cerebellar_scatter.pdf")
strep_cerebellar_scatter <- ggscatter(
  strepsirr, y = "MedianCerebellum", x = "Mass",
  color = "Strepsirrhine", palette = "jco",
  add = "reg.line"
) +
  facet_wrap(~Strepsirrhine) +
  stat_cor(label.y = 5.4) +
  stat_regline_equation(label.y = 5.3)
strep_cerebellar_scatter
graphics.off()

## Cerebrum
strepsirr %>%
  group_by(Strepsirrhine) %>%
  summarize(cor=cor(MedianCerebrum, `Mass`))

pdf("1.strepsirrhine/1.0c.cerebral.pdf")
strep_cerebral <- ggplot(strepsirr, aes(x = Mass, y = MedianCerebrum, color = Strepsirrhine)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("grp", "eq", "R2"))) +
  geom_point()
strep_cerebral
graphics.off()

pdf("1.strepsirrhine/1.0d.cerebral_scatter.pdf")
strep_cerebral_scatter <- ggscatter(
  strepsirr, y = "MedianCerebrum", x = "Mass",
  color = "Strepsirrhine", palette = "jco",
  add = "reg.line"
) +
  facet_wrap(~Strepsirrhine) +
  stat_cor(label.y = 6.2) +
  stat_regline_equation(label.y = 6.1)
strep_cerebral_scatter
graphics.off()

## Ansiform area
strepsirr %>%
  group_by(Strepsirrhine) %>%
  summarize(cor=cor(MedianCrus, `Mass`))

pdf("1.strepsirrhine/1.0e.ansiform.pdf")
strep_AA <- ggplot(strepsirr, aes(x = Mass, y = MedianCrus, color = Strepsirrhine)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("grp", "eq", "R2"))) +
  geom_point()
strep_AA
graphics.off()

pdf("1.strepsirrhine/1.0f.ansiform_scatter.pdf")
strep_AA_scatter <- ggscatter(
  strepsirr, y = "MedianCrus", x = "Mass",
  color = "Strepsirrhine", palette = "jco",
  add = "reg.line"
) +
  facet_wrap(~Strepsirrhine) +
  stat_cor(label.y = 4.5) +
  stat_regline_equation(label.y = 4.4)
strep_AA_scatter
graphics.off()

########################################################
### REGULAR REGRESSION - Brain Volume
########################################################

# Cerebellum
strep_cerebellar_bv <- ggplot(strepsirr, aes(x = BrainVol, y = MedianCerebellum, color = Strepsirrhine, grp.label = Strepsirrhine)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("grp", "eq", "R2"))) +
  geom_point()
strep_cerebellar_bv

strep_cerebellar_bv_scatter <- ggscatter(
  strepsirr, y = "MedianCerebellum", x = "BrainVol",
  color = "Strepsirrhine", palette = "jco",
  add = "reg.line"
) +
  facet_wrap(~Strepsirrhine) +
  stat_cor(label.y = 5.4) +
  stat_regline_equation(label.y = 5.3)
strep_cerebellar_bv_scatter

## Cerebrum
strep_cerebral_bv <- ggplot(strepsirr, aes(x = BrainVol, y = MedianCerebrum, color = Strepsirrhine, grp.label = Strepsirrhine)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("grp", "eq", "R2"))) +
  geom_point()
strep_cerebellar_bv

strep_cerebral_bv_scatter <- ggscatter(
  strepsirr, y = "MedianCerebrum", x = "BrainVol",
  color = "Strepsirrhine", palette = "jco",
  add = "reg.line"
) +
  facet_wrap(~Strepsirrhine) +
  stat_cor(label.y = 6.2) +
  stat_regline_equation(label.y = 6.1)
strep_cerebral_bv_scatter

##################################
### PIC/ PGLS REGRESSION
#################################

# Plot trees
plotTree(tree)
plotTree(tree_strep)
plotTree(tree_hapl)

# Check if node ages are appropriate
tree.age(tree)
tree.age(tree_strep)
tree.age(tree_hapl)

#####################################
### FIRST: strepsirrhines
#####################################

Cstrep <- STREP$MedianCerebellum
Vstrep <- STREP$MedianCerebrum
Mstrep <- STREP$Mass 

# Compute contrasts
# Full
pic.C.strep <- pic(Cstrep,tree_strep)
pic.V.strep <- pic(Vstrep, tree_strep)
pic.M.strep <- pic(Mstrep, tree_strep)


#------------------------------------
# PIC-regressions
#------------------------------------
# Compute a linear regression of C on V
result <- lm(pic.C.strep ~ pic.V.strep)
sink("1.strepsirrhine/1.1.lm-C-V.strepsirrhine.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("1.strepsirrhine/1.2.lm-C-V-mean-nointercept-plot.strepsirrhine.pdf")
plot(pic.C.strep ~ pic.V.strep)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.
## They are the same; check! ##

## Check if PIC/PGLS models are the same.
bm.strep <-corBrownian(1,tree_strep)
pglsModel<-gls(MedianCerebellum~MedianCerebrum,data=STREP,correlation=bm.strep)
summary(pglsModel)
PICmodel <- lm(pic.C.strep ~pic.V.strep + 0)
sink("1.strepsirrhine/1.3.PGLS-PIC-C-V.strepsirrhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


## Let's also examine strepsirrhine cerebellar and cerebral scaling versus body mass. This might help to illustrate the relative small cerebra in this lineage.
## Cerebellum
pglsModel<-gls(MedianCerebellum~Mass,data=STREP,correlation=bm.strep)
summary(pglsModel)
PICmodel <- lm(pic.C.strep~pic.M.strep + 0)
sink("1.strepsirrhine/1.4.PGLS-PIC-C-M.strepsirrhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

## Cerebrum
pglsModel<-gls(MedianCerebrum~Mass,data=STREP,correlation=bm.strep)
summary(pglsModel)
PICmodel <- lm(pic.V.strep ~pic.M.strep + 0)
sink("1.strepsirrhine/1.5.PGLS-PIC-V-M.strepsirrhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


####################################
### NEXT: haplorhines
####################################
Chapl <- HAPL$MedianCerebellum
Vhapl <- HAPL$MedianCerebrum
AAhapl <- HAPLAA$MedianCrus
ROChapl <- HAPLAA$MedianROC
Mhapl <- HAPL$Mass
BVhapl <- HAPL$BrainVol # is the number of observations the same?

# Compute contrasts
# Full
pic.C.hapl <- pic(Chapl,tree_hapl)
pic.V.hapl <- pic(Vhapl, tree_hapl)
pic.AA.hapl <- pic(AAhapl, tree_haplaa)
pic.ROC.hapl <- pic(ROChapl, tree_haplaa)
pic.M.hapl <- pic(Mhapl, tree_hapl)
pic.BV.hapl <- pic(BVhapl, tree_hapl) # does this work?

#------------------------------------
# PIC-regression
#------------------------------------
# Compute a linear regression of C on V
result <- lm(pic.C.hapl ~ pic.V.hapl)
sink("1.strepsirrhine/1.6.lm-C-V.haplorhine.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("1.strepsirrhine/1.7.lm-C-V-mean-nointercept-plot.haplorhine.pdf")
plot(pic.C.hapl ~ pic.V.hapl)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.

## Compare PIC regressions and PGLS.
bm.hapl <-corBrownian(1,tree_hapl)
pglsModel<-gls(MedianCerebellum~MedianCerebrum,data=HAPL,correlation=bm.hapl)
summary(pglsModel)
PICmodel <- lm(pic.C.hapl ~pic.V.hapl + 0)
sink("1.strepsirrhine/1.8.PGLS-PIC-C-V.haplorhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

## Let's also examine haplorhine cerebellar and cerebral scaling versus body mass. This might help to illustrate the relative small cerebra in this lineage.
## Cerebellum
pglsModel<-gls(MedianCerebellum~Mass,data=HAPL,correlation=bm.hapl)
summary(pglsModel)
PICmodel <- lm(pic.C.hapl~pic.M.hapl + 0)
sink("1.strepsirrhine/1.9.PGLS-PIC-C-M.haplorhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

## We also measured brain volume in most our specimens. For several reason this is more illustrative for relative scaling of brain areas in both groups.
## Cerebellum
pglsModel<-gls(MedianCerebellum~BrainVol,data=HAPL,correlation=bm.hapl)
summary(pglsModel)
PICmodel <- lm(pic.C.hapl~pic.BV.hapl + 0)
sink("1.strepsirrhine/1.9a.PGLS-PIC-C-BV.haplorhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


## Cerebrum
pglsModel<-gls(MedianCerebrum~Mass,data=HAPL,correlation=bm.hapl)
summary(pglsModel)
PICmodel <- lm(pic.V.hapl ~pic.M.hapl + 0)
sink("1.strepsirrhine/1.10.PGLS-PIC-V-M.haplorhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

## Brain volumes cerebrum
pglsModel<-gls(MedianCerebrum~BrainVol,data=HAPL,correlation=bm.hapl)
summary(pglsModel)
PICmodel <- lm(pic.V.hapl ~pic.BV.hapl + 0)
sink("1.strepsirrhine/1.10a.PGLS-PIC-V-BV.haplorhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

#######################################################################
## Ansiform area regression for the haplorhines only (12/13 species)
#######################################################################

# Compute a linear regression of AA on ROC
result <- lm(pic.AA.hapl ~ pic.ROC.hapl)
sink("1.strepsirrhine/1.11.lm-AA-ROC.haplorhine.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("1.strepsirrhine/1.12.lm-AA-ROC-mean-nointercept-plot.haplorhine.pdf")
plot(pic.AA.hapl ~ pic.ROC.hapl)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.

## Compare PIC regressions and PGLS.
bm.haplaa <-corBrownian(1,tree_haplaa)
pglsModel<-gls(MedianCrus~MedianROC,data=HAPLAA,correlation=bm.haplaa)
summary(pglsModel)
PICmodel <- lm(pic.AA.hapl ~pic.ROC.hapl + 0)
sink("1.strepsirrhine/1.13.PGLS-PIC-AA-ROC.haplorhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

###########################################################
### Use pANCOVA to test deviations in intercept or slope
### and with that deviation from the general allometric 
### relationship for cerebellar-to-cerebral and cerebellar
### & cerebral-to-body mass scaling
###########################################################

library(evomap)
library(geiger)

# This is an alternative way to subset the trees :). Very convenient
Strepsirrhines <-getTips(tree,findMRCA(tree,c("Galago_demidoff","Microcebus_murinus")))
Haplorhines <-getTips(tree,findMRCA(tree,c("Saimiri_sciureus","Cercopithecus_cephus_cephus")))

#For differences in slope: 
grpS<-rep("A",length(rownames(strepsirr))) 
grpS[Strepsirrhines]<-"B" 
grpS<-as.factor(grpS) 
names(grpS)<-rownames(strepsirr)

#For differences in intercept: 
grpI<-rep("A",length(rownames(strepsirr))) 
grpI[Strepsirrhines]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(strepsirr)


## Data subsetting
## Cerebellum vs. Cerebrum
data<-primateData; tree<-primateTree
Y<-"Brain"; X<-"Body"
data<-data[,c(which(colnames(data)==Y),which(colnames(data)==X)),drop=F]
data<-na.omit(data); data<-log(data)
tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)   #match the data to the tree
colnames(data)<-c("Dependent1", "Dependent2")  
## Very hacky but bypasses weird error.
data <- data[1:33,]
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
tree <-treedata(tree,strepsirr, sort=T,warnings=T)$ph
rownames(data) <- tree[["tip.label"]]

## Plot separate PGLS regressions for Strepsirrhines and Haplorhines
## Again, extremely hacky, but works.
data$Dependent1 <- strepsirr$MedianCerebellum
data$Dependent2 <- strepsirr$MedianCerebrum
data$Independent <- strepsirr$Mass

## Cerebellum vs. Cerebrum
pdf("1.strepsirrhine/1.14.PGLS.Strepsirrhine.CC.pdf")
plot(data$Dependent1~data$Dependent2,col="black",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent1","Dependent2",data,tree,model="BM",group=Strepsirrhines,col="plum3",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent1","Dependent2",data,tree,model="BM",group=Haplorhines,col="grey15",lwd=5,cex=1,pch=19)
graphics.off()

## Model construction
Model.CC <-model.matrix(as.formula(Dependent1~Dependent2),data)
Model.CC.s <-model.matrix(as.formula(Dependent1~grpS:Dependent2),data) #slopes differ
Model.CC.i <-model.matrix(as.formula(Dependent1~grpI + Dependent2),data) #intercepts differ
Model.CC.si <-model.matrix(as.formula(Dependent1~grpI + grpS:Dependent2),data) #slopes & intercepts differ

sink("1.strepsirrhine/1.15.ModelTesting.SlopesIntercepts.CC.txt")
## Testing versus baseline model (common intercept and slope)
print("(1) Differences in slopes, holding intercept constant:")
gls.ancova(Dependent1~Dependent2,vcv(tree),Model.CC,Model.CC.s)

print("(2) Differences in intercept, holding slopes constant:")
gls.ancova(Dependent1~Dependent2,vcv(tree),Model.CC,Model.CC.i)

print("(3) Differences in slopes and differences in intercept:")
gls.ancova(Dependent1~Dependent2,vcv(tree),Model.CC,Model.CC.si)
sink()


## Cerebellum vs. Body Mass
pdf("1.strepsirrhine/1.16.PGLS.Strepsirrhine.CM.pdf")
plot(data$Dependent1~data$Independent,col="black",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent1","Independent",data,tree,model="BM",group=Strepsirrhines,col="pink",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent1","Independent",data,tree,model="BM",group=Haplorhines,col="azure4",lwd=5,cex=1,pch=19)
graphics.off()

## Model construction
Model.CM <-model.matrix(as.formula(Dependent1~Independent),data)
Model.CM.s <-model.matrix(as.formula(Dependent1~grpS:Independent),data) #slopes differ
Model.CM.i <-model.matrix(as.formula(Dependent1~grpI + Independent),data) #intercepts differ
Model.CM.si <-model.matrix(as.formula(Dependent1~grpI + grpS:Independent),data) #slopes & intercepts differ

sink("1.strepsirrhine/1.17.ModelTesting.SlopesIntercepts.CM.txt")
## Testing versus baseline model (common intercept and slope)
print("(1) Differences in slopes, holding intercept constant:")
gls.ancova(Dependent1~Independent,vcv(tree),Model.CM,Model.CM.s)

print("(2) Differences in intercept, holding slopes constant:")
gls.ancova(Dependent1~Independent,vcv(tree),Model.CM,Model.CM.i)

print("(3) Differences in slopes and differences in intercept:")
gls.ancova(Dependent1~Independent,vcv(tree),Model.CM,Model.CM.si)
sink()

## Cerebrum vs. Body Mass
## The important question is if the intercept is different between the two groups in this PGLS.
pdf("1.strepsirrhine/1.18.PGLS.Strepsirrhine.VM.pdf")
plot(data$Dependent2~data$Independent,col="black",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent2","Independent",data,tree,model="BM",group=Strepsirrhines,col="pink",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent2","Independent",data,tree,model="BM",group=Haplorhines,col="azure4",lwd=5,cex=1,pch=19)
graphics.off()

## Model construction
Model.VM <-model.matrix(as.formula(Dependent2~Independent),data)
Model.VM.s <-model.matrix(as.formula(Dependent2~grpS:Independent),data) #slopes differ
Model.VM.i <-model.matrix(as.formula(Dependent2~grpI + Independent),data) #intercepts differ
Model.VM.si <-model.matrix(as.formula(Dependent2~grpI + grpS:Independent),data) #slopes & intercepts differ

sink("1.strepsirrhine/1.19.ModelTesting.SlopesIntercepts.VM.txt")


## Testing versus baseline model (common intercept and slope)
print("(1) Differences in slopes, holding intercept constant:")
gls.ancova(Dependent2~Independent,vcv(tree),Model.VM,Model.VM.s)

print("(2) Differences in intercept, holding slopes constant:")
gls.ancova(Dependent2~Independent,vcv(tree),Model.VM,Model.VM.i)

print("(3) Differences in slopes and differences in intercept:")
gls.ancova(Dependent2~Independent,vcv(tree),Model.VM,Model.VM.si)
sink()



## Plot them together
pdf("1.strepsirrhine/1.20.PGLS.Strepsirrhine.composite.pdf", width = 16, height = 8)
plot(data$Dependent1~data$Dependent2,col="black",pch=10,xlab="", asp=0.4, ylab="",cex.lab=2, ylim =c(2.0, 6.2))
pGLS.plotGrade("Dependent1","Dependent2",data,tree,model="BM",group=Strepsirrhines,col="plum3",lwd=5,cex=1,pch=10)
pGLS.plotGrade("Dependent1","Dependent2",data,tree,model="BM",group=Haplorhines,col="grey15",lwd=5,cex=1,pch=10)
## Add cerebrum - brainvol
pGLS.plotGrade("Dependent2","Independent",data,tree,model="BM",group=Strepsirrhines,col="pink",lwd=5,cex=1,pch=10, lty="dotted")
pGLS.plotGrade("Dependent2","Independent",data,tree,model="BM",group=Haplorhines,col="azure4",lwd=5,cex=1,pch=10,  lty="dotted")
## Add cerebellum - brainvol
pGLS.plotGrade("Dependent1","Independent",data,tree,model="BM",group=Strepsirrhines,col="pink",lwd=5,cex=1,pch=10,  lty="dotted")
pGLS.plotGrade("Dependent1","Independent",data,tree,model="BM",group=Haplorhines,col="azure4",lwd=5,cex=1,pch=10,  lty="dotted")
graphics.off()


### -------------------------- ###
### END of script
### -------------------------- ### 




