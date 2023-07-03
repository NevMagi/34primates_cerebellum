### -------------------------- ###
# Phylogenetic regressions 
# (PGLSs & PIC Regressions)
# N. Magielse, 2022
### -------------------------- ###     

## Load packages
library(ape)
#library(geiger)
library(nlme)
library(phytools)
library(corrplot)
library(tidyverse)
library(magrittr)
library(evomap)
#library(car) #maybe not necessary (trying some things)
library(rr2) # For PGLS R2, also see R2s for Correlated Data: Phylogenetic Models, LMMs, and GLMMs from Anthony Ives.
#library(phylolm)
#library(relaimpo)
library(ggpmisc)

## Set working directory
getwd()
setwd("~/Documents/PhD/EvolutionPrimates/Final_code/") ## Provide path to folder with scripts and input data.

## Set output directory
outputdir <- "5.phylogenetic_regressions"
dir.create(outputdir)

## Load data
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
tree13 <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_13.nex"), method="extend") # To be used for those species with full values (no NAs for Ansiform Area Volume)
pheno <- read.csv2("./0.cleanInput/Phenotypes_outliersRM.csv", sep=';')  %>% dplyr::rename(species = WilsonReederName) %>% dplyr::select(species, Hominoidea, CerebellarVol, SurfaceArea, CerebralVol, AbsGI, FoldingLength, FoldingNumber, Lambda, Delta, CrusVol_1_3, ROC)
phenoAA <- pheno %>% drop_na()
phenoMedian <- read.csv2("./0.cleanInput/PhenotypesMedian_outliersRM.csv", sep=';')  %>% dplyr::rename(species = WilsonReederName) %>% dplyr::select(species, Hominoidea, MedianCerebellum, MedianCerebrum, MedianCrus, MedianROC)
phenoMedianAA <- phenoMedian %>% drop_na()

# The function `pic.ortho` accepts multiple observations per individual.
# They are grouped by species name, as the tips of the tree
phenolist <- split(pheno,pheno$species)
phenolist13 <- split(phenoAA, phenoAA$species)
phenolist.Median <- split(phenoMedian, phenoMedian$species)
phenolist.Median13 <- split(phenoMedianAA, phenoMedianAA$species)

# Reorder to match tip label order
phenolist <- phenolist[tree$tip.label];
phenolist13 <- phenolist13[tree13$tip.label]
phenolist.Median <- phenolist.Median[tree$tip.label]
phenolist.Median13 <- phenolist.Median13[tree13$tip.label]

# Get phenotypes
C <- lapply(phenolist,"[[","CerebellarVol")
SA <- lapply(phenolist,"[[","SurfaceArea")
V <- lapply(phenolist,"[[","CerebralVol")
G <- lapply(phenolist,"[[","AbsGI")
L <- lapply(phenolist,"[[","FoldingLength")
N <- lapply(phenolist,"[[","FoldingNumber")
W <- lapply(phenolist,"[[","Lambda")
D <- lapply(phenolist,"[[","Delta")
AA <- lapply(phenolist,"[[","CrusVol_1_3")
ROC <- lapply(phenolist,"[[","ROC")

# Get phenotypes for species without NAs anywhere.
C13 <- lapply(phenolist13,"[[","CerebellarVol")
SA13 <- lapply(phenolist13,"[[","SurfaceArea")
V13 <- lapply(phenolist13,"[[","CerebralVol")
G13 <- lapply(phenolist13,"[[","AbsGI")
L13 <- lapply(phenolist13,"[[","FoldingLength")
N13 <- lapply(phenolist13,"[[","FoldingNumber")
W13 <- lapply(phenolist13,"[[","Lambda")
D13 <- lapply(phenolist13,"[[","Delta")
AA13 <- lapply(phenolist13,"[[","CrusVol_1_3")
ROC13 <- lapply(phenolist13,"[[","ROC")

# Get phenotypes for species median values (only Cerebellar and Cerebral volumes of interest)
C.med <- lapply(phenolist.Median,"[[","MedianCerebellum")
V.med <- lapply(phenolist.Median,"[[","MedianCerebrum")
AA.med <- lapply(phenolist.Median,"[[","MedianCrus")
ROC.med <- lapply(phenolist.Median,"[[","MedianROC")

C.med13 <- lapply(phenolist.Median13,"[[","MedianCerebellum")
V.med13 <- lapply(phenolist.Median13,"[[","MedianCerebrum")
AA.med13 <- lapply(phenolist.Median13,"[[","MedianCrus")
ROC.med13 <- lapply(phenolist.Median13,"[[","MedianROC")

#---------------
# Plot the tree
#---------------

# Plot tree annotated with branch length
pdf('5.phylogenetic_regressions/5.1.tree.pdf')
plot(tree, cex=0.5)
axisPhylo()
nodelabels(cex=0.5)
title("Branch length (Mya)")
edgelabels(round(tree$edge.length,1), frame="n", adj=c(0.5,-0.2), cex=0.5)
graphics.off()

# Plot 13-species tree annotated with branch length
pdf('5.phylogenetic_regressions/5.1a.tree13.pdf')
plot(tree13, cex=0.5)
axisPhylo()
nodelabels(cex=0.5)
title("Branch length (Mya)")
edgelabels(round(tree$edge.length,1), frame="n", adj=c(0.5,-0.2), cex=0.5)
graphics.off()

# Plot tree as variance-covariance matrix
pdf('5.phylogenetic_regressions/5.2.varcovmatrix.pdf')
heatmap(vcv.phylo(tree))
graphics.off()

# Plot 13-species tree as variance-covariance matrix
pdf('5.phylogenetic_regressions/5.2a.varcovmatrix13.pdf')
heatmap(vcv.phylo(tree13))
graphics.off()


#------------------------------------
# Phylogenetic independent contrasts
#------------------------------------

# Compute contrasts
pic.C <- pic.ortho(C,tree,intra=TRUE)
pic.SA <- pic.ortho(SA,tree,intra=TRUE)
pic.V <- pic.ortho(V,tree,intra=TRUE)
pic.G <- pic.ortho(G,tree,intra=TRUE)
pic.L <- pic.ortho(L,tree,intra=TRUE)
pic.N <- pic.ortho(N,tree,intra=TRUE)
pic.W <- pic.ortho(W,tree,intra=TRUE)
pic.D <- pic.ortho(D,tree,intra=TRUE)
pic.AA <- pic.ortho(AA, tree, intra=TRUE) ## This one does not really make sense. We are only able to calculate PIC for nodes 51 and 52. PICs are more sensibly calculated from the 13-species tree, that has complete data.
pic.ROC <- pic.ortho(ROC, tree, intra=TRUE) ## Ibidem.

# Compute contrasts for 13-species data
pic.C13 <- pic.ortho(C13,tree13,intra=TRUE)
pic.SA13 <- pic.ortho(SA13,tree13,intra=TRUE)
pic.V13 <- pic.ortho(V13,tree13,intra=TRUE)
pic.G13 <- pic.ortho(G13,tree13,intra=TRUE)
pic.L13 <- pic.ortho(L13,tree13,intra=TRUE)
pic.N13 <- pic.ortho(N13,tree13,intra=TRUE)
pic.W13 <- pic.ortho(W13,tree13,intra=TRUE)
pic.D13 <- pic.ortho(D13,tree13,intra=TRUE)
pic.AA13 <- pic.ortho(AA13, tree13,  intra=TRUE)
pic.ROC13 <- pic.ortho(ROC13, tree13, intra=TRUE)

# Compute contrast for species median values (to compare PIC-model forced through 0 and the PGLS-model)
pic.C.med <- pic(C.med, tree)
pic.V.med <- pic(V.med, tree)
#pic.AA.med <- pic(AA.med, tree)
#pic.ROC.med <- pic(ROC.med, tree)

pic.C13.med <- pic(C.med13, tree13)
pic.V13.med <- pic(V.med13, tree13)
pic.AA13.med <- pic(AA.med13, tree13)
pic.ROC13.med <- pic(ROC.med13, tree13)

# Plot tree annotated with C and V contrasts
pdf('5.phylogenetic_regressions/5.3.tree-C-V.pdf')
plot(tree, cex=0.5)
axisPhylo()
title("Cerebellar and cerebral volume phylogenetic independent contrasts")
nodelabels(round(pic.C,2),adj=c(-0.05,-0.25), cex=0.35, frame="n")
nodelabels(round(pic.V,2),adj=c(-0.05,1.25),cex=0.35, frame="n")
graphics.off()

# Plot tree annotated with C and V contrasts, for 13 species
pdf('5.phylogenetic_regressions/5.3a.tree13-C-V.pdf')
plot(tree13, cex=0.5)
axisPhylo()
title("Cerebellar and cerebral volume phylogenetic independent contrasts (13 species)")
nodelabels(round(pic.C13,2),adj=c(-0.05,-0.25), cex=0.35, frame="n")
nodelabels(round(pic.V13,2),adj=c(-0.05,1.25),cex=0.35, frame="n")
graphics.off()

# Plot scatterplots for all PICs.
pdf('5.phylogenetic_regressions/5.4.scatterplots-PICs.pdf')
df=data.frame(pic.C,pic.V,pic.SA,pic.G,pic.L,pic.N,pic.W,pic.D)
colnames(df) <- c("Cerebellar \n Volume","Cerebral \n Volume", "Cerebral \n Surface \n Area","AbsGI","Folding \n Length","Folding \n Number","Fold \n Wavelength","Fold Depth")
pairs(~pic.C+pic.V+pic.SA+pic.G+pic.L+pic.N+pic.W+pic.D,
      data=df,
      labels=colnames(df),
      panel=panel.smooth)
graphics.off()

# Plot scatterplots for 13-species PICs (full (non-NAs) data, has ansiform area as well).
pdf('5.phylogenetic_regressions/5.4a.scatterplots-PICs13.pdf')
df=data.frame(pic.C13,pic.V13,pic.AA13, pic.SA13,pic.G13,pic.L13,pic.N13,pic.W13,pic.D13)
colnames(df) <- c("Cerebellar \n Volume","Cerebral \n Volume", "Ansiform Area \n Volume", "Cerebral \n Surface \n Area","AbsGI","Folding \n Length","Folding \n Number","Fold \n Wavelength","Fold Depth")
pairs(~pic.C13+pic.V13+pic.AA13+pic.SA13+pic.G13+pic.L13+pic.N13+pic.W13+pic.D13,
      data=df,
      labels=colnames(df),
      panel=panel.smooth)
graphics.off()

# Plot correlation matrix for all PICs.
M <- cbind(pic.C,pic.V,pic.SA,pic.G,pic.L,pic.N,pic.W,pic.D)
R <- cor(M)
rownames(R) <- colnames(R) <- c("Cerebellar Volume", "Cerebral Surface Area","Cerebral Volume","AbsGI","Folding Length","Folding Number","Fold Wavelength","Fold Depth")
pdf('5.phylogenetic_regressions/5.5.correlations.pdf')
corrplot(R, method="ellipse", order="AOE")
graphics.off()

# Plot correlation matrix for 13-species PICs.
M <- cbind(pic.C13,pic.V13,pic.AA13, pic.SA13,pic.G13,pic.L13,pic.N13,pic.W13,pic.D13)
R <- cor(M)
rownames(R) <- colnames(R) <- c("Cerebellar Volume", "Ansiform Area Volume","Cerebral Surface Area","Cerebral Volume","AbsGI","Folding Length","Folding Number","Fold Wavelength","Fold Depth")
pdf('5.phylogenetic_regressions/5.5a.correlations13.pdf')
corrplot(R, method="ellipse", order="AOE")
graphics.off()

# Compute a linear regression of C on V
result <- lm(pic.C ~ pic.V)
sink("5.phylogenetic_regressions/5.6.lm-C-V.txt")
summary(result)
sink()

# Compute a linear regression of C on V, 13 species
result <- lm(pic.C13 ~ pic.V13)
sink("5.phylogenetic_regressions/5.6a.lm-C13-V13.txt")
summary(result)
sink()

# Compute a linear regression of AA on V, 13 species
result <- lm(pic.AA13 ~ pic.V13)
sink("5.phylogenetic_regressions/5.6b.lm-AA13-V13.txt")
summary(result)
sink()

# Compute a linear regression of AA on C, 13 species
result <- lm(pic.AA13 ~ pic.C13)
sink("5.phylogenetic_regressions/5.6c.lm-AA13-C13.txt")
summary(result)
sink()

# Compute a linear regression of AA on ROC, 13 species
result <- lm(pic.AA13 ~ pic.ROC13)
sink("5.phylogenetic_regressions/5.6d.lm-AA13-ROC13.txt")
summary(result)
sink()

# PICs have expected mean zero: fit a regression with intercept fixed at 0
# Linear regression of C on V
result <- lm(pic.C ~ 0 +  pic.V)
sink("5.phylogenetic_regressions/5.7.lm-C-V-mean-nointercept.txt")
summary(result)
sink()

# Linear regression of C on V, 13 species
result <- lm(pic.C13 ~ 0 + pic.V13)
sink("5.phylogenetic_regressions/5.7a.lm-C-V-mean-nointercept13.txt")
summary(result)
sink()

# Linear regression of AA on V, 13 species
result <- lm(pic.AA13 ~ 0 + pic.V13)
sink("5.phylogenetic_regressions/5.7b.lm-AA-V-mean-nointercept13.txt")
summary(result)
sink()

# Linear regression of AA on C, 13 species
result <- lm(pic.AA13 ~ 0 +  pic.C13)
sink("5.phylogenetic_regressions/5.7c.lm-AA-C-mean-nointercept13.txt")
summary(result)
sink()

# Linear regression of AA on ROC, 13 species
result <- lm(pic.AA13 ~ 0 +  pic.ROC13)
sink("5.phylogenetic_regressions/5.7d.lm-AA-ROC-mean-nointercept13.txt")
summary(result)
sink()


# Plot PIC regressions
# Regression of C on V
pdf("5.phylogenetic_regressions/5.8.lm-C-V-mean-nointercept-plot.pdf")
plot(pic.C ~ pic.V)
abline(a = 0, b = coef(result))
graphics.off()

# Regression of C on V, 13 species
pdf("5.phylogenetic_regressions/5.8a.lm-C-V-mean-nointercept-plot13.pdf")
plot(pic.C13 ~ pic.V13)
abline(a = 0, b = coef(result))
graphics.off()

# Regression of AA on V, 13 species
pdf("5.phylogenetic_regressions/5.8b.lm-AA-V-mean-nointercept-plot13.pdf")
plot(pic.AA13 ~ pic.V13)
abline(a = 0, b = coef(result))
graphics.off()

# Regression of AA on C, 13 species
pdf("5.phylogenetic_regressions/5.8c.lm-AA-C-mean-nointercept-plot13.pdf")
plot(pic.AA13 ~ pic.C13)
abline(a = 0, b = coef(result))
graphics.off()

# Regression of AA on ROC, 13 species
pdf("5.phylogenetic_regressions/5.8d.lm-AA-ROC-mean-nointercept-plot13.pdf")
plot(pic.AA13 ~ pic.ROC13)
abline(a = 0, b = coef(result))
graphics.off()


## PGLS under Brownian Motion should be the same as PIC-regression forced through 0.
# We use species median data here, so that neither the PGLS and PIC-model incorporate intraspecific variance.
bm<-corBrownian(1,tree)
pglsModel<-gls(MedianCerebellum~MedianCerebrum,data=phenoMedian,correlation=bm)
summary(pglsModel)
PICmodel <- lm(pic.C.med ~pic.V.med + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("5.phylogenetic_regressions/5.9.PGLS-PIC-C-V.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

## Allow Pagel's lambda to vary, using corPagel.
## We test three different models, starting at lambda = 0.8, and comparing it with ANOVAs to models starting at lamda = 0.0 and 1.0, respectively.
sink("5.phylogenetic_regressions/5.9a.PGLS-PIC-C-V-variable-lambda.txt")
PGLS.CV <-gls(MedianCerebellum~MedianCerebrum,correlation=corPagel(value=0.8,phy=tree),data=phenoMedian) 
intervals(PGLS.CV,which="var-cov")
PGLS.CV.l0 <- gls(MedianCerebellum~MedianCerebrum,correlation=corPagel(value=0.0,phy=tree, fixed = T),data=phenoMedian) 
PGLS.CV.l1 <- gls(MedianCerebellum~MedianCerebrum,correlation=corPagel(value=1.0,phy=tree, fixed = T),data=phenoMedian) 
anova(PGLS.CV, PGLS.CV.l0)
anova(PGLS.CV, PGLS.CV.l1)
sink()

pdf("5.phylogenetic_regressions/5.9b.PGLS-PIC-C-V-lambda-distribution.pdf")
lambda<-seq(0,1,length.out=500)
lik<-sapply(lambda,function(lambda)logLik(gls(MedianCerebellum~MedianCerebrum, correlation=corPagel(value=lambda,phy=tree,fixed=TRUE), data=phenoMedian)))
plot(lik~lambda,type="l",main=expression(paste("LikelihoodPlotfor", lambda)),ylab="LogLikelihood",xlab=expression(lambda))
abline(v=PGLS.CV$modelStruct,col="black")
graphics.off()

# Also for only 13 species with complete data.
# We take species median values.
sink("5.phylogenetic_regressions/5.9c.PGLS-PIC-C13-V13.txt")
bm13<-corBrownian(1,tree13)
pglsModel<-gls(MedianCerebellum~MedianCerebrum,data=phenoMedianAA,correlation=bm13)
PICmodel <- lm(pic.C13.med ~pic.V13.med + 0)
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

## PGLS and PIC-models for Ansiform Area regressed on Cerebellar Volume
# Both the PGLS and PIC-model incorporate intraspecific variation (PGLS at this step, PIC at the pic.ortho step)
pglsModel <- gls(MedianCrus ~ MedianCerebellum, correlation = corBrownian(phy = tree13),
                 data = phenoMedianAA, method = "ML")
PICmodel <- lm(pic.AA13.med~ 0 + pic.C13.med) ## NOT YET THE SAME AS PGLS, check what is up.
sink("5.phylogenetic_regressions/5.10.PGLS-PIC-AA-C.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

sink("5.phylogenetic_regressions/5.10a.PGLS-PIC-AA-C-variable-lambda.txt")
PGLS.AA <-gls(MedianCrus~MedianCerebellum, correlation=corPagel(value=0.8,phy=tree13),data=phenoMedianAA) 
intervals(PGLS.AA,which="var-cov")
sink()

pdf("5.phylogenetic_regressions/5.10b.PGLS-PIC-AA-C-lambda-distribution.pdf")
lambda<-seq(0,1,length.out=500)
lik<-sapply(lambda,function(lambda)logLik(gls(MedianCrus~MedianCerebellum, correlation=corPagel(value=lambda,phy=tree13,fixed=TRUE), data=phenoMedianAA)))
plot(lik~lambda,type="l",main=expression(paste("LikelihoodPlotfor", lambda)),ylab="LogLikelihood",xlab=expression(lambda))
abline(v=PGLS.AA$modelStruct,col="black")
graphics.off()

## PGLS and PIC-models for Ansiform Area regressed on Rest of Cerebellar Volume (ROC)
# Both the PGLS and PIC-model incorporate intraspecific variation (PGLS at this step, PIC at the pic.ortho step)
pglsModel <- gls(MedianCrus ~ MedianROC, correlation = corBrownian(phy = tree13),
                 data = phenoMedianAA, method = "ML")
PICmodel <- lm(pic.AA13.med~ 0 + pic.ROC13.med)
sink("5.phylogenetic_regressions/5.10c.PGLS-PIC-AA-ROC.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

sink("5.phylogenetic_regressions/5.10d.PGLS-PIC-AA-ROC-variable-lambda.txt")
PGLS.AA <-gls(MedianCrus~MedianROC, correlation=corPagel(value=0.8,phy=tree13),data=phenoMedianAA) 
intervals(PGLS.AA,which="var-cov")
sink()

pdf("5.phylogenetic_regressions/5.10e.PGLS-PIC-AA-ROC-lambda-distribution.pdf")
lambda<-seq(0,1,length.out=500)
lik<-sapply(lambda,function(lambda)logLik(gls(MedianCrus~MedianROC, correlation=corPagel(value=lambda,phy=tree13,fixed=TRUE), data=phenoMedianAA)))
plot(lik~lambda,type="l",main=expression(paste("LikelihoodPlotfor", lambda)),ylab="LogLikelihood",xlab=expression(lambda))
abline(v=PGLS.AA$modelStruct,col="black")
graphics.off()

## PGLS and PIC-models for Ansiform Area regressed on Cerebral Volume
# Both the PGLS and PIC-model incorporate intraspecific variation (PGLS at this step, PIC at the pic.ortho step)
pglsModel <- gls(MedianCrus ~ MedianCerebrum, correlation = corBrownian(phy = tree13),
                 data = phenoMedianAA, method = "ML")
PICmodel <- lm(pic.AA13.med ~ 0 + pic.V13.med)
sink("5.phylogenetic_regressions/5.10f.PGLS-PIC-AA-V.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

sink("5.phylogenetic_regressions/5.10g.PGLS-PIC-AA-V-variable-lambda.txt")
PGLS.AA <-gls(MedianCrus~MedianCerebrum,correlation=corPagel(value=0.8,phy=tree13),data=phenoMedianAA) 
intervals(PGLS.AA,which="var-cov")
sink()

pdf("5.phylogenetic_regressions/5.10h.PGLS-PIC-AA-V-lambda-distribution.pdf")
lambda<-seq(0,1,length.out=500)
lik<-sapply(lambda,function(lambda)logLik(gls(MedianCrus~MedianCerebellum, correlation=corPagel(value=lambda,phy=tree13,fixed=TRUE), data=phenoMedianAA)))
plot(lik~lambda,type="l",main=expression(paste("LikelihoodPlotfor", lambda)),ylab="LogLikelihood",xlab=expression(lambda))
abline(v=PGLS.AA$modelStruct,col="black")
graphics.off()

### -------------------------- ###
### PGLS Plotting with Caper
### -------------------------- ###   
# phenoMedian.cd <- comparative.data(tree, phenoMedian, species, vcv=TRUE, vcv.dim=3)
bm<-corBrownian(1,tree)
bm13 <-corBrownian(1,tree13)

# Model 1: Cerebellar volume regressed on cerebral volume.
tiff("./5.phylogenetic_regressions/5.11.BM-PGLS-cerebellumCerebrum-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod1 <- gls(MedianCerebellum ~ MedianCerebrum, data=phenoMedian, correlation=bm)
plot(mod1)
plot(phenoMedian$MedianCerebellum ~ phenoMedian$MedianCerebrum, xlab = "Cerebral volume (log10-transformed)", ylab = "Cerebellar volume (log10-transformed)")
abline(mod1)
abline(a=-0.6270014, b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod1.ci <-gls.ci(phenoMedian$MedianCerebellum, phenoMedian$MedianCerebrum,vcv(tree))
lines(pGLS.mod1.ci$CI.plot$X,pGLS.mod1.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod1.ci$CI.plot$X,pGLS.mod1.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod1.pi <-gls.pi(phenoMedian$MedianCerebellum, phenoMedian$MedianCerebrum,vcv(tree), 1)
lines(pGLS.mod1.pi$PI.plot$X,pGLS.mod1.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod1.pi$PI.plot$X,pGLS.mod1.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()

# Model 2: Cerebellar volume regressed on cerebral volume, only for the 13 species with complete data.
tiff("./5.phylogenetic_regressions/5.11a.BM-PGLS-cerebellumCerebrum13-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod2 <- gls(MedianCerebellum ~ MedianCerebrum, data=phenoMedianAA, correlation=bm13)
plot(mod2)
plot(phenoMedianAA$MedianCerebellum ~ phenoMedianAA$MedianCerebrum, xlab = "Cerebral volume (log10-transformed)", ylab = "Cerebellar volume (log10-transformed)")
abline(mod2)
abline(a=-0.5641225, b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod2.ci <-gls.ci(phenoMedianAA$MedianCerebellum, phenoMedianAA$MedianCerebrum,vcv(tree13))
lines(pGLS.mod2.ci$CI.plot$X,pGLS.mod2.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod2.ci$CI.plot$X,pGLS.mod2.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod2.pi <-gls.pi(phenoMedianAA$MedianCerebellum, phenoMedianAA$MedianCerebrum,vcv(tree13), 1)
lines(pGLS.mod2.pi$PI.plot$X,pGLS.mod2.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod2.pi$PI.plot$X,pGLS.mod2.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()

# Model 3: Ansiform area volume regressed on cerebellar volume.
tiff("./5.phylogenetic_regressions/5.12.BM-PGLS-ansiformareaCerebellum-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod3 <- gls(MedianCrus ~ MedianCerebellum, data=phenoMedianAA, correlation=bm13)
plot(mod3)
plot(phenoMedianAA$MedianCrus ~ phenoMedianAA$MedianCerebellum, xlab = "Cerebellar volume (log10-transformed)", ylab = "Ansiform area volume (log10-transformed)")
abline(mod3)
abline(a=-1.976440, b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod3.ci <-gls.ci(phenoMedianAA$MedianCrus, phenoMedianAA$MedianCerebellum,vcv(tree13))
lines(pGLS.mod3.ci$CI.plot$X,pGLS.mod3.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod3.ci$CI.plot$X,pGLS.mod3.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod3.pi <-gls.pi(phenoMedianAA$MedianCrus, phenoMedianAA$MedianCerebellum,vcv(tree13), 1)
lines(pGLS.mod3.pi$PI.plot$X,pGLS.mod3.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod3.pi$PI.plot$X,pGLS.mod3.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()

# Model 3-better: Ansiform area volume regressed on rest of cerebellar volume.
tiff("./5.phylogenetic_regressions/5.12a.BM-PGLS-ansiformareaROC-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod3a <- gls(MedianCrus ~ MedianROC, data=phenoMedianAA, correlation=bm13)
plot(mod3a)
plot(phenoMedianAA$MedianCrus ~ phenoMedianAA$MedianROC, xlab = "Rest of cerebellar volume (log10-transformed)", ylab = "Ansiform area volume (log10-transformed)")
abline(mod3a)
abline(a=-1.833182, b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod3a.ci <-gls.ci(phenoMedianAA$MedianCrus, phenoMedianAA$MedianCerebellum,vcv(tree13))
lines(pGLS.mod3a.ci$CI.plot$X,pGLS.mod3a.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod3a.ci$CI.plot$X,pGLS.mod3a.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod3a.pi <-gls.pi(phenoMedianAA$MedianCrus, phenoMedianAA$MedianCerebellum,vcv(tree13), 1)
lines(pGLS.mod3a.pi$PI.plot$X,pGLS.mod3a.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod3a.pi$PI.plot$X,pGLS.mod3a.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()

# Model 4: Ansiform area volume regressed on cerebral volume.
tiff("./5.phylogenetic_regressions/5.13.BM-PGLS-ansiformareaCerebrum-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod4 <- gls(MedianCrus ~ MedianCerebrum, data=phenoMedianAA, correlation=bm13)
plot(mod4)
plot(phenoMedianAA$MedianCrus ~ phenoMedianAA$MedianCerebrum,  xlab = "Cerebral volume (log10-transformed)", ylab = "Ansiform area volume (log10-transformed)")
abline(mod4)
abline(a=-2.783793, b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod4.ci <-gls.ci(phenoMedianAA$MedianCrus, phenoMedianAA$MedianCerebrum,vcv(tree13))
lines(pGLS.mod4.ci$CI.plot$X,pGLS.mod4.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod4.ci$CI.plot$X,pGLS.mod4.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod4.pi <-gls.pi(phenoMedianAA$MedianCrus, phenoMedianAA$MedianCerebrum,vcv(tree13), 1)
lines(pGLS.mod4.pi$PI.plot$X,pGLS.mod4.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod4.pi$PI.plot$X,pGLS.mod4.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()


### -------------------------- ###
### Calculate R2s for PGLS models
### -------------------------- ###   

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

### -------------------------- ###
### PGLS Apes vs. Nonapes
### -------------------------- ###   
library(lsmeans)
library(psych)
library(data.table)

# Data
# 34 species
phenoMedian.dt <- phenoMedian
rownames(phenoMedian.dt) <- phenoMedian.dt$species
phenoMedian.dt <- phenoMedian.dt %>% dplyr::select(-MedianCrus, -species) %>% as.data.table()
# 13 species
phenoMedianAA.dt <- phenoMedianAA
rownames(phenoMedianAA.dt) <- phenoMedianAA.dt$species
phenoMedianAA.dt <- as.data.table(phenoMedianAA.dt)

## Models ##

sink("./5.phylogenetic_regressions/5.11b.BM-PGLS-cerebellumCerebrum-interaction.apes.txt")
# Mod1
mod1.interaction <- gls(MedianCerebellum ~ MedianCerebrum*Hominoidea, data=phenoMedian.dt, correlation=bm)
anova(mod1.interaction) # No significant interaction ape/nonape*Cerebral volume

# Obtain slopes
mod1.interaction$coefficients
mod1.lst <- lstrends(mod1.interaction, "Hominoidea", var="MedianCerebrum")
# Compare slopes
pairs(mod1.lst) # No significant interaction ape/nonape*Cerebral volume

# Difference in residual variances
# Calculate Pearson's R
mod1.correlation <- R2.lik(mod = mod1) # Calculate R2
print(mod1.correlation)
mod1.ape <- gls(MedianCerebellum ~ MedianCerebrum, data=subset(phenoMedian.dt, Hominoidea == "Ape"), correlation=bm)
mod1.ape.correlation <- R2.lik(mod = mod1.ape) # Calculate R2
print(mod1.ape.correlation)
mod1.nonape <- gls(MedianCerebellum ~ MedianCerebrum, data=subset(phenoMedian.dt, Hominoidea == "Nonape"), correlation=bm)
mod1.nonape.correlation <- R2.lik(mod = mod1.nonape) # Calculate R2
print(mod1.nonape.correlation)

# Compare R values with Fisher's R to Z
# Ape vs. Full
paired.r(mod1.correlation, mod1.ape.correlation, 
         n = phenoMedian.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # Splitting provides significantly worse correlation.
# Nonape vs. Full
paired.r(mod1.correlation, mod1.nonape.correlation, 
         n = phenoMedian.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # Splitting provides significantly worse correlation.
# Ape vs. Nonape
paired.r(mod1.nonape.correlation, mod1.ape.correlation, 
         n = phenoMedian.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # no difference
sink()

sink("./5.phylogenetic_regressions/5.11c.BM-PGLS-cerebellumCerebrum-interaction.apes13.txt")
# Mod2
mod2.interaction <- gls(MedianCerebellum ~ MedianCerebrum*Hominoidea, data=phenoMedianAA, correlation=bm)
anova(mod2.interaction) # No significant interaction ape/nonape*Cerebral volume

# Obtain slopes
mod2.interaction$coefficients
mod2.lst <- lstrends(mod2.interaction, "Hominoidea", var="MedianCerebrum")
# Compare slopes
pairs(mod2.lst) # No significant interaction ape/nonape*Cerebral volume

# Difference in residual variances
# Calculate Pearson's R
# mod1.correlation <- phenoMedian.dt[, cor(MedianCerebellum, MedianCerebrum), by = Hominoidea] #Only to create comparable structure.
mod2.correlation <- R2.lik(mod = mod2) # Calculate R2
print(mod2.correlation)
mod2.ape <- gls(MedianCerebellum ~ MedianCerebrum, data=subset(phenoMedianAA.dt, Hominoidea == "Ape"), correlation=bm)
mod2.ape.correlation <- R2.lik(mod = mod2.ape) # Calculate R2
print(mod2.ape.correlation)
mod2.nonape <- gls(MedianCerebellum ~ MedianCerebrum, data=subset(phenoMedianAA.dt, Hominoidea == "Nonape"), correlation=bm)
mod2.nonape.correlation <- R2.lik(mod = mod2.nonape) # Calculate R2
print(mod2.nonape.correlation)

# Compare R values with Fisher's R to Z
# Ape vs. Full
paired.r(mod2.correlation, mod2.ape.correlation, 
         n = phenoMedianAA.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # Splitting provides significantly worse correlation.
# Nonape vs. Full
paired.r(mod2.correlation, mod2.nonape.correlation, 
         n = phenoMedianAA.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # Splitting provides trend to worse correlation.
# Ape vs. Nonape
paired.r(mod2.nonape.correlation, mod2.ape.correlation, 
         n = phenoMedianAA.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # no difference
sink()


sink("./5.phylogenetic_regressions/5.12b.BM-PGLS-ansiformareaCerebellum-interaction.apes13.txt")
# Mod3
mod3.interaction <- gls(MedianCrus ~ MedianCerebellum*Hominoidea, data=phenoMedianAA.dt, correlation=bm)
anova(mod3.interaction) # No significant interaction ape/nonape*Cerebellar volume

# Obtain slopes
mod3.interaction$coefficients
mod3.lst <- lstrends(mod3.interaction, "Hominoidea", var="MedianCerebellum")
# Compare slopes
pairs(mod3.lst) # No significant interaction ape/nonape*Cerebellar volume

# Difference in residual variances
# Calculate Pearson's R
# mod1.correlation <- phenoMedian.dt[, cor(MedianCerebellum, MedianCerebrum), by = Hominoidea] #Only to create comparable structure.
mod3.correlation <- R2.lik(mod = mod3) # Calculate R2
print(mod3.correlation)
mod3.ape <- gls(MedianCrus ~ MedianCerebellum, data=subset(phenoMedianAA.dt, Hominoidea == "Ape"), correlation=bm)
mod3.ape.correlation <- R2.lik(mod = mod3.ape) # Calculate R2
print(mod3.ape.correlation)
mod3.nonape <- gls(MedianCrus~ MedianCerebellum, data=subset(phenoMedianAA.dt, Hominoidea == "Nonape"), correlation=bm)
mod3.nonape.correlation <- R2.lik(mod = mod3.nonape) # Calculate R2
print(mod3.nonape.correlation)

# Compare R values with Fisher's R to Z
# Ape vs. Full
paired.r(mod3.correlation, mod3.ape.correlation, 
         n = phenoMedianAA.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # No difference
# Nonape vs. Full
paired.r(mod3.correlation, mod3.nonape.correlation, 
         n = phenoMedianAA.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # No difference
# Ape vs. Nonape
paired.r(mod3.nonape.correlation, mod3.ape.correlation, 
         n = phenoMedianAA.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # Trend towards better correlation for splitting the data.
sink()

sink("./5.phylogenetic_regressions/5.12c.BM-PGLS-ansiformareaROC-interaction.apes13.txt")
# Mod3
mod3a.interaction <- gls(MedianCrus ~ MedianROC*Hominoidea, data=phenoMedianAA.dt, correlation=bm)
anova(mod3a.interaction) # No significant interaction ape/nonape*Cerebellar volume

# Obtain slopes
mod3a.interaction$coefficients
mod3a.lst <- lstrends(mod3a.interaction, "Hominoidea", var="MedianROC")
# Compare slopes
pairs(mod3a.lst) # No significant interaction ape/nonape*Cerebellar volume

# Difference in residual variances
# Calculate Pearson's R
# mod1.correlation <- phenoMedian.dt[, cor(MedianCerebellum, MedianCerebrum), by = Hominoidea] #Only to create comparable structure.
mod3a.correlation <- R2.lik(mod = mod3a) # Calculate R2
print(mod3a.correlation)
mod3a.ape <- gls(MedianCrus ~ MedianROC, data=subset(phenoMedianAA.dt, Hominoidea == "Ape"), correlation=bm)
mod3a.ape.correlation <- R2.lik(mod = mod3a.ape) # Calculate R2
print(mod3a.ape.correlation)
mod3a.nonape <- gls(MedianCrus~ MedianROC, data=subset(phenoMedianAA.dt, Hominoidea == "Nonape"), correlation=bm)
mod3a.nonape.correlation <- R2.lik(mod = mod3a.nonape) # Calculate R2
print(mod3a.nonape.correlation)

# Compare R values with Fisher's R to Z
# Ape vs. Full
paired.r(mod3a.correlation, mod3a.ape.correlation, 
         n = phenoMedianAA.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # No difference
# Nonape vs. Full
paired.r(mod3a.correlation, mod3a.nonape.correlation, 
         n = phenoMedianAA.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # No difference
# Ape vs. Nonape
paired.r(mod3a.nonape.correlation, mod3a.ape.correlation, 
         n = phenoMedianAA.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # Trend towards better correlation for splitting the data.
sink()


sink("./5.phylogenetic_regressions/5.13c.BM-PGLS-ansiformareaCerebrum-interaction.apes13.txt")
# Mod4
mod4.interaction <- gls(MedianCrus ~ MedianCerebrum*Hominoidea, data=phenoMedianAA.dt, correlation=bm)
anova(mod4.interaction) # No significant interaction ape/nonape*Cerebral volume

# Obtain slopes
mod4.interaction$coefficients
mod4.lst <- lstrends(mod4.interaction, "Hominoidea", var="MedianCerebrum")
# Compare slopes
pairs(mod4.lst) # No significant interaction ape/nonape*Cerebral volume

# Difference in residual variances
# Calculate Pearson's R
# mod1.correlation <- phenoMedian.dt[, cor(MedianCerebellum, MedianCerebrum), by = Hominoidea] #Only to create comparable structure.
mod4.correlation <- R2.lik(mod = mod4) # Calculate R2
print(mod4.correlation)
mod4.ape <- gls(MedianCrus ~ MedianCerebrum, data=subset(phenoMedianAA.dt, Hominoidea == "Ape"), correlation=bm)
mod4.ape.correlation <- R2.lik(mod = mod4.ape) # Calculate R2
print(mod4.ape.correlation)
mod4.nonape <- gls(MedianCrus~ MedianCerebrum, data=subset(phenoMedianAA.dt, Hominoidea == "Nonape"), correlation=bm)
mod4.nonape.correlation <- R2.lik(mod = mod4.nonape) # Calculate R2
print(mod4.nonape.correlation)

# Compare R values with Fisher's R to Z
# Ape vs. Full
paired.r(mod4.correlation, mod4.ape.correlation, 
         n = phenoMedianAA.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # No difference
# Nonape vs. Full
paired.r(mod4.correlation, mod4.nonape.correlation, 
         n = phenoMedianAA.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # No difference
# Ape vs. Nonape
paired.r(mod4.nonape.correlation, mod4.ape.correlation, 
         n = phenoMedianAA.dt[Hominoidea %in% c("Ape", "Nonape"), .N]) # No difference
sink()


### ----------------------------------------------------------- ###
### Regular linear regression.
### ----------------------------------------------------------- ### 

sink("5.phylogenetic_regressions/5.15.regular_regressions_allometries.txt")
print("Cerebellum regressed on cerebrum.")
reg.mod1 <- lm(MedianCerebellum ~ MedianCerebrum, data = phenoMedian)
summary(reg.mod1)

## Test what remains of this relationship in the smaller dataset, that only contains those species with full data (incl. Ansiform Area).
# This helps in interpretation of models 2 and 3, that depend on this more limited dataset.
print("Cerebellum regressed on cerebrum, 13 species, complete data (no NAs.")
reg.mod2<- lm(MedianCerebellum ~ MedianCerebrum, data = phenoMedianAA)
summary(reg.mod2)

print("Ansiform area regressed on cerebellum.")
reg.mod3 <- lm(MedianCrus ~ MedianCerebellum, data = phenoMedianAA)
summary(reg.mod3)

print("Ansiform area regressed on rest of cerebellum.")
reg.mod3a <- lm(MedianCrus ~ MedianROC, data = phenoMedianAA)
summary(reg.mod3a)

print("Ansiform area regressed on cerebrum.")
reg.mod4 <- lm(MedianCrus ~ MedianCerebrum, data = phenoMedianAA)
summary(reg.mod4)
sink()


### ---------------------------------------------------------------- ###
### GGPLOT Regular linear regressions plots, with allometric formulae
### ---------------------------------------------------------------- ### 

# Credit: Johnston, S: sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/.
# Define function
ggplotRegression <- function (fit, ylabel, xlabel) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "black") +
    labs(title = paste("y = ", 
                       signif(fit$coef[[1]],5 ), 
                       "+ ",
                       signif(fit$coef[[2]], 5),
                       "x",
                       "         R2 Adj. = ",signif(summary(fit)$adj.r.squared, 5),
                       " p-value =",signif(summary(fit)$coef[2,4], 5)), y = ylabel, x = xlabel)
}


tiff("./5.phylogenetic_regressions/5.16.regularRegression-cerebellumCerebrum.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplotRegression(reg.mod1, ylabel = "Median cerebellar volume (log10-transformed)", xlabel = "Median cerebral volume (log10-transformed)") + theme_classic()  +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
graphics.off()

tiff("./5.phylogenetic_regressions/5.16a.regularRegression-cerebellumCerebrum_13species.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplotRegression(reg.mod2, ylabel = "Median cerebellar volume (log10-transformed)", xlabel = "Median cerebral volume (log10-transformed)") + theme_classic()  +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
graphics.off()

tiff("./5.phylogenetic_regressions/5.17.regularRegression-ansiformCerebellum.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplotRegression(reg.mod3, ylabel = "Median ansiform area volume (log10-transformed)", xlabel = "Median cerebellar volume (log10-transformed)") + theme_classic()  +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
graphics.off()

tiff("./5.phylogenetic_regressions/5.17a.regularRegression-ansiformROC.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplotRegression(reg.mod3a, ylabel = "Median ansiform area volume (log10-transformed)", xlabel = "Median rest of cerebellar volume (log10-transformed)") + theme_classic()  +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
graphics.off()


tiff("./5.phylogenetic_regressions/5.18.regularRegression-ansiformCerebrum.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplotRegression(reg.mod4, ylabel = "Median ansiform area volume (log10-transformed)", xlabel = "Median cerebral volume (log10-transformed)") + theme_classic()  +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
graphics.off()



### -------------------------------------------------------- ###
### Analyse apes separately
### -------------------------------------------------------- ### 
phenoMedian.ape <- phenoMedian %>% dplyr::filter(Hominoidea == "Ape")
phenoMedian.nonape <- phenoMedian %>% dplyr::filter(Hominoidea == "Nonape")

phenoMedianAA.ape <- phenoMedianAA %>% dplyr::filter(Hominoidea == "Ape")
phenoMedianAA.nonape <- phenoMedianAA %>% dplyr::filter(Hominoidea == "Nonape")

reg.mod1.ape <- lm(MedianCerebellum ~ MedianCerebrum, data = phenoMedian.ape)
summary(reg.mod1.ape)
reg.mod1.nonape <- lm(MedianCerebellum ~ MedianCerebrum, data = phenoMedian.nonape)
summary(reg.mod1.nonape)

reg.mod2.ape <- lm(MedianCerebellum ~ MedianCerebrum, data = phenoMedianAA.ape)
summary(reg.mod2.ape)
reg.mod2.nonape <- lm(MedianCerebellum ~ MedianCerebrum, data = phenoMedianAA.nonape)
summary(reg.mod2.nonape)

reg.mod3.ape <- lm(MedianCrus ~ MedianCerebellum, data = phenoMedianAA.ape)
summary(reg.mod3.ape)
reg.mod3.nonape <- lm(MedianCrus ~ MedianCerebellum, data = phenoMedianAA.nonape)
summary(reg.mod3.nonape)

reg.mod3a.ape <- lm(MedianCrus ~ MedianROC, data = phenoMedianAA.ape)
summary(reg.mod3a.ape)
reg.mod3a.nonape <- lm(MedianCrus ~ MedianROC, data = phenoMedianAA.nonape)
summary(reg.mod3a.nonape)

reg.mod4.ape <- lm(MedianCrus ~ MedianCerebrum, data = phenoMedianAA.ape)
summary(reg.mod4.ape)
reg.mod4.nonape <- lm(MedianCrus ~ MedianCerebrum, data = phenoMedianAA.nonape)
summary(reg.mod4.nonape)

### ---------------------------------------------------------------- ###
### GGPLOT
### ---------------------------------------------------------------- ### 
tiff("./5.phylogenetic_regressions/5.19.regularRegression-cerebellumCerebrum-apes.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplot(phenoMedian, aes(y=MedianCerebellum, x=MedianCerebrum, colour = Hominoidea)) + 
  geom_point(na.rm=T, aes(color=phenoMedian$Hominoiodea)) +
  geom_smooth(method="lm", na.rm=T, se=F, formula=y~x, aes(color=phenoMedian$Hominoidea)) +
  theme_classic() +
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25)) +
  labs(y = "", x = "") +
  theme(legend.position="top") +
  theme(legend.text = element_text(size=30)) +
  theme(legend.title = element_text(size=30)) +
  scale_colour_discrete("Hominoidea") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), size = 9)
graphics.off()


tiff("./5.phylogenetic_regressions/5.19a.regularRegression-cerebellumCerebrum_13species-apes.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplot(phenoMedianAA, aes(y=MedianCerebellum, x=MedianCerebrum, colour = Hominoidea)) + 
  geom_point(na.rm=T, aes(color=phenoMedianAA$Hominoiodea)) +
  geom_smooth(method="lm", na.rm=T, se=F, formula=y~x, aes(color=phenoMedianAA$Hominoidea)) +
  theme_classic() +
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25)) +
  labs(y = "", x = "") +
  theme(legend.position="top") +
  theme(legend.text = element_text(size=30)) +
  theme(legend.title = element_text(size=30)) +
  scale_colour_discrete("Hominoidea") +
  guides(color = guide_legend(override.aes = list(size = 5)))  +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), size = 9)
graphics.off()

tiff("./5.phylogenetic_regressions/5.20.regularRegression-ansiformCerebellum-apes.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplot(phenoMedianAA, aes(y=MedianCrus, x=MedianCerebellum, colour = Hominoidea)) + 
  geom_point(na.rm=T, aes(color=phenoMedianAA$Hominoiodea)) +
  geom_smooth(method="lm", na.rm=T, se=F, formula=y~x, aes(color=phenoMedianAA$Hominoidea)) +
  theme_classic() +
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25)) +
  labs(y = "", x = "") +
  theme(legend.position="top") +
  theme(legend.text = element_text(size=30)) +
  theme(legend.title = element_text(size=30)) +
  scale_colour_discrete("Hominoidea") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), size = 9)
graphics.off()

tiff("./5.phylogenetic_regressions/5.20a.regularRegression-ansiformROC-apes.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplot(phenoMedianAA, aes(y=MedianCrus, x=MedianROC, colour = Hominoidea)) + 
  geom_point(na.rm=T, aes(color=phenoMedianAA$Hominoiodea)) +
  geom_smooth(method="lm", na.rm=T, se=F, formula=y~x, aes(color=phenoMedianAA$Hominoidea)) +
  theme_classic() +
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25)) +
  labs(y = "", x = "") +
  theme(legend.position="top") +
  theme(legend.text = element_text(size=30)) +
  theme(legend.title = element_text(size=30)) +
  scale_colour_discrete("Hominoidea") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), size = 9)
graphics.off()

tiff("./5.phylogenetic_regressions/5.21.regularRegression-ansiformCerebrum-apes.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplot(phenoMedianAA, aes(y=MedianCrus, x=MedianCerebrum, colour= Hominoidea)) + 
  geom_point(na.rm=T, aes(color=phenoMedianAA$Hominoiodea)) +
  geom_smooth(method="lm", na.rm=T, se=F, formula=y~x, aes(color=phenoMedianAA$Hominoidea)) +
  theme_classic() +
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25)) +
  labs(y = "", x = "") +
  theme(legend.position="top") +
  theme(legend.text = element_text(size=30)) +
  theme(legend.title = element_text(size=30)) +
  scale_colour_discrete("Hominoidea") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), size = 9)
graphics.off()


### -------------------------------------------------------- ###
### Although we do not correct all analyses for body mass for
### several reasons, (volatile body weight data that was spar-
### sely sampled; potential loss of signal, focus on relative
### organization of the cerebello-cerebral system) it seems
### that primate ansiform area scaling becomes larger as a di-
### rect result of larger cerebella, which are closely-related
### to larger brains and body mass. To show that the relative
### size of the ansiform area increases with body mass, we cor-
### relate these variables.
### -------------------------------------------------------- ### 

## Data input
CCmass <- read_csv2("./0.cleanInput/PhenotypesMedian.csv") %>% 
  dplyr::rename(mass = `Body mass species mean`) %>%
  filter(!is.na(mass))
CCmass$CCrat <- log10(CCmass$MedianCerebroCerebellar) 

AAmass <- CCmass %>% filter(!is.na(MedianCerebellarCrus))
AAmass <- AAmass[-c(9,10),]
AAmass$AArat <- log10(AAmass$MedianCerebellarCrus) 

AAmass <- read_csv2("./0.cleanInput/PhenotypesMedian.csv") %>% 
  dplyr::rename(mass = `Body mass species mean`) %>%
  dplyr::filter(!is.na(mass)) %>% filter(!is.na(MedianCerebellarCrus))

## Correlation - CC Ratio - BodyMass)
cor.test(CCmass$MedianCerebroCerebellar, CCmass$mass, method=c("pearson"))

## Correlation - AA Ratio - BodyMass)
cor.test(AAmass$MedianCerebellarCrus, AAmass$mass, method=c("pearson"))

### -------------------------- ###
### Plotting regressions
### -------------------------- ### 
library("ggpubr")
ggscatter(CCmass, x = "mass", y = "MedianCerebroCerebellar",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Body Mass (in kilograms; log10-transformed)", ylab = "Ratio cerebellum/cerebrum")

ggscatter(AAmass, x = "mass", y = "MedianCerebellarCrus", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Body Mass (in kilograms; log10-transformed)", ylab = "Ratio ansiform area/cerebellum")


## log10-transformed ratios ## 
ggscatter(CCmass, x = "mass", y = "CCrat", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Body Mass (in kilograms; log10-transformed)", ylab = "Ratio cerebellum/cerebrum (loge-transformed")

ggscatter(AAmass, x = "mass", y = "AArat", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Body Mass (in kilograms; log10-transformed)", ylab = "Ratio ansiform area/cerebellum (log10-transformed)")


## Normality ## 
sink("5.phylogenetic_regressions/5.22.normality_BWtraits.txt")
# Shapiro-Wilk normality test CC
shapiro.test(CCmass$MedianCerebroCerebellar) # => p = 0.07 = NORMAL
# Shapiro-Wilk normality test AA
shapiro.test(AAmass$MedianCerebellarCrus) # => p = 0.14 = NORMAL
# Shapiro-Wilk normality test for mass
shapiro.test(CCmass$mass) # => p = 0.59 = NORMAL
shapiro.test(AAmass$mass) # => p = 0.46 = NORMAL
sink()

### -------------------------- ###
### With lm. Separating apes/
### non-apes.
### -------------------------- ###

# Model construction
reg.mod1.MASS <- lm(MedianCerebroCerebellar ~ mass, data = CCmass)
reg.mod2.MASS <- lm(MedianCerebellarCrus ~ mass, data = AAmass)


sink("5.phylogenetic_regressions/5.23.regular_regressions_BW.txt")
print("Ratios regressed on body mass.")
summary(reg.mod1.MASS)
summary(reg.mod2.MASS)
sink()

tiff("./5.phylogenetic_regressions/5.24.regularRegression-CCmass.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplotRegression(reg.mod1.MASS, ylabel = "Cerebellar/cerebral volume ratio", xlabel = "Body mass (log10-transformed)") + theme_classic()  +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
graphics.off()

tiff("./5.phylogenetic_regressions/5.25.regularRegression-AAmass.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplotRegression(reg.mod2.MASS, ylabel = "Ansiform Area/cerebellar volume ratio", xlabel = "Body mass (log10-transformed)") + theme_classic()  +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
graphics.off()


### -------------------------- ###
### END of script
### -------------------------- ###   