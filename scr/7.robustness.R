### ------------------------------------- ###
### Robustness of results, 
# replication in Stephan et al. (1981) data.
### Neville Magielse, 2023
### ------------------------------------- ###        

### SET UP ###
getwd()
setwd("~/Documents/PhD/EvolutionPrimates/Final_code/")

outputdir <- "7.robustness"
dir.create(outputdir)

library(tidyverse)
library(dplyr)
library(magrittr)
library(ape)
library(data.table)

lut <- read_csv2("./input/10kTrees_34PrimateSpecies_adapted.csv")  #look-up table
Stephan <- read_csv2("./input/Stephan1981-34species.csv") %>% dplyr::select(-ID, -n) %>% as.data.frame()
#%>% dplyr::rename(species = WilsonReederName) #Cerebellar and cerebral volumes integrated with Stephan et al. (1981).
Stephan.log <- Stephan
Stephan.log$CerebellarVol <- log10(Stephan$CerebellarVol)
Stephan.log$CerebralVol <- log10(Stephan$CerebralVol)


# Now, coincidentally, there are also 34 primate species that match between the Stephan collection and 10kTrees
# .. as of January 13th, 2023. We will call this tree 'tree34S'.
tree34S <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_Stephan34.nex"), method="extend")
plotTree(tree34S)


### -------------------------- ###
# 1. Evolutionary model selection
# Roberto Toro, Katja Heuer 2018
# Script adapted by N. Magielse for inclusion of cerebellar data & robustness analysis 2022
### -------------------------- ###     

## SET-UP ## 
library(Rphylopars)
library(phytools)
library(magrittr)
library(tidyverse)

source("./input/catn.R")

## Fit Brownian Model model.
p_BM <- phylopars(trait_data=Stephan.log, tree=tree34S)

# Testing BM vs. Star-model.
catn("Fit and compare Pagel's lambda model for lambda=1 (BM) and lambda=0 (star)")
p_lambda <- phylopars(trait_data=Stephan.log, tree=tree34S, usezscores = F)
p_star <- phylopars(trait_data=Stephan.log, tree=tree34S, model="star", usezscores = F)
sink("7.robustness/7.1.test-pagel-lamdba=1-vs-lambda=0_full.txt")
p_lambda
p_star
catn('Log-likelihood full-BM', logLik(p_lambda), sep='\t')
catn('Log-likelihood star',logLik(p_star), sep='\t')
chi_square <- as.double(2*(logLik(p_lambda) - logLik(p_star))) # 2*(logLik_alt - logLik_null)
catn('chi_square', chi_square, sep='\t')
degrees_freedom <- p_lambda$npars - p_star$npars # df = difference in number of model parameters
catn('df', degrees_freedom, sep='\t')
p <- pchisq(q = chi_square,df = degrees_freedom,lower.tail = FALSE) # p-value
catn ('p-value', p, sep='\t')
sink()

## Fit the Ornstein-Uhlenbeck model ##
#------------------------------------
# Univariate
# Save phylogenetic and phenotypic trait variance-covariance matrices, variance explained by phylogeny
catn("Fit Ornstein-Uhlenbeck model (OU, univariate: a single alpha value shared by all traits)")
p_OU <- phylopars(trait_data = Stephan.log, tree=tree34S, model = "OU", usezscores = F) # Something is off here, veryh high AIC values (seems unrealistic).
sink("7.robustness/7.2.OU-univariate_full.txt")
p_OU
sink()

# Multivariate #DOES NOT WORK
catn("Fit multivariate OU (a diagonal matrix of alpha values)")
p_mvOU <- phylopars(trait_data = Stephan.log,tree = tree34S, model = "mvOU", full_alpha = F, usezscores = F)
sink("7.robustness/7.3.OU-multivariate_full.txt")
p_mvOU
sink()

#Diagonal Multivariate #DOES NOT WORK
catn("Fit multivariate OU (a diagonal matrix of alpha values)")
p_mvOU_diag <- phylopars(trait_data = Stephan.log, tree = tree34S, model = "mvOU", usezscores = F)
sink("7.robustness/7.4.OU-multivariate-diagonal_full.txt")
p_mvOU_diag
sink()

## Fit the Early-Burst model ##
#------------------------------------
catn("Fit early-burst model (EB)")
p_EB <- phylopars(trait_data = Stephan.log, tree = tree34S, model = "EB", usezscores = F)
sink("7.robustness/7.5.EB_full.txt")
p_EB # Estimated trait covariance and EB rate parameter
sink()


# Model selection  
#------------------------------------
catn("Model selection using AIC")
sink("7.robustness/7.6.model-selection_full.txt")
catn("Brownian motion", AIC(p_BM), sep='\t')
catn("Ornstein-Uhlenbeck, single alpha", AIC(p_OU), sep='\t')
catn("Ornstein-Uhlenbeck, diagonal alpha matrix", AIC(p_mvOU), sep='\t')
catn("Ornstein-Uhlenbeck, full alpha matrix", AIC(p_mvOU_diag), sep='\t')
catn("Early burst", AIC(p_EB), sep='\t')
catn("Star model, lambda = 0", AIC(p_star), sep='\t')
sink()



### ---------------------------------- ###
# 2. Ancestral state reconstructions
# Roberto Toro, Katja Heuer 2018
# Script adapted by Neville Magielse 
# for inclusion of cerebellar and 
# ansiform area data, and inclusion
# of stephan et al. (1981) values - 2023
### --------------------------------- ###     

## We reuse the trees and data from part 1. ## 
#tree = tree34S
#Data = Stephan

# Create Brownian motion model
p.BM <- phylopars(trait_data=Stephan.log, tree=tree34S)
# Create Brownian motion model with complete star-phylogeny.
p.BM.star <- phylopars(trait_data=Stephan.log, tree=tree34S, model = 'star')

#------------------------------------------
# SAVE ANCESTRAL STATE ESTIMATIONS TO CSV
#------------------------------------------
# ntips Number of species we have phenotypes for
ntips <- length(tree34S$tip.label)
# Number of nodes inside the tree (ntips-1) + including the tips
nnodes <- ntips + tree34S$Nnode

## Put into csv the info for all the estimated phenotypes at each internal node (=all ancestors)
catn("save data in csv")
csv <- cbind(p.BM$anc_recon[(ntips+1):nnodes,])
# lower and upper confidence interval
lci <- csv - sqrt(p.BM$anc_var[(ntips+1):nnodes,])*1.96
uci <- csv + sqrt(p.BM$anc_var[(ntips+1):nnodes,])*1.96
colnames(lci) <- paste(colnames(lci), "-95%CI", sep = "")
colnames(uci) <- paste(colnames(uci), "+95%CI", sep = "")
csv <- cbind(csv,lci)
csv <- cbind(csv,uci)
write.csv(csv, file="./7.robustness/7.7.ACE_fullmodel.csv")

# Create Brownian motion model wit raw values
p.BM.raw <- phylopars(trait_data=Stephan, tree=tree34S)
# Create Brownian motion model with complete star-phylogeny.
p.BM.star.raw <- phylopars(trait_data=Stephan, tree=tree34S, model = 'star')

# Make a table with raw values as well.
catn("save data in csv")
csv <- cbind(p.BM.raw$anc_recon[(ntips+1):nnodes,])
# lower and upper confidence interval
lci <- csv - sqrt(p.BM.raw$anc_var[(ntips+1):nnodes,])*1.96
uci <- csv + sqrt(p.BM.raw$anc_var[(ntips+1):nnodes,])*1.96
colnames(lci) <- paste(colnames(lci), "-95%CI", sep = "")
colnames(uci) <- paste(colnames(uci), "+95%CI", sep = "")
csv <- cbind(csv,lci)
csv <- cbind(csv,uci)
write.csv(csv, file="./7.robustness/7.7a.ACE_fullmodel.raw.csv")

#------------------------------------------
# PLOT ANCESTRAL STATE ESTIMATIONS
#------------------------------------------

# Create plotting functions again.
myplotanc.ci.Taxo <- function(tree, tips, tipnames, states, legend.title) {
  names(tips) <- tipnames
  obj <- contMap(tree,tips,method="user",anc.states=states, plot="F")
  obj <- setMap(obj, colors=c("blue", "white", "red"))
  plot.contMap(obj, legend=FALSE,ylim=c(1-0.09*(Ntip(obj$tree)-1),Ntip(obj$tree)),
               mar=c(5.1,0.4,0.4,0.4))
  add.color.bar(36.502,obj$cols,title= legend.title, y=3,
                lims=obj$lims,digits=2,prompt=FALSE,x=0,
                y=1-0.08*(Ntip(obj$tree)-1),lwd=4,fsize=1,subtitle="")
  axisPhylo()
  title(xlab="Time to present (in million years)", adj = 0.17)
  errorbar.contMap(obj, lwd=2)
}

myplotanc.ci.noTaxo <- function(tree, tips, tipnames, states, legend.title) {
  names(tips) <- tipnames
  obj <- contMap(tree,tips,method="user",anc.states=states, plot="F", ftype = "off")
  obj <- setMap(obj, colors=c("blue", "white", "red"))
  plot.contMap(obj, legend=FALSE,ylim=c(1-0.09*(Ntip(obj$tree)-1),Ntip(obj$tree)),
               mar=c(5.1,0.4,0.4,0.4), ftype = "off")
  add.color.bar(36.502,obj$cols,title= legend.title, y=3,
                lims=obj$lims,digits=2,prompt=FALSE,x=0,
                y=1-0.08*(Ntip(obj$tree)-1),lwd=4,fsize=0.8,subtitle="")
  errorbar.contMap(obj, lwd=2)
}

catn("plot ancestral values")
tipnames <- rownames(p.BM$anc_recon)[1:ntips]

# Cerebellum
tiff("./7.robustness/7.8.BM-anc-cerebellarVolume-ci.tiff", width = 4, height = 8, units = 'in', res = 300)
myplotanc.ci.noTaxo(tree34S, p.BM$anc_recon[1:ntips,"CerebellarVol"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"CerebellarVol"], legend.title = "Mean cerebellar volume \n (log10 mm3)")
graphics.off()

# Cerebrum
tiff("./7.robustness/7.9.BM-anc-cerebralVolume-ci.tiff", width = 4, height = 8, units = 'in', res = 300)
myplotanc.ci.noTaxo(tree34S, p.BM$anc_recon[1:ntips,"CerebralVol"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"CerebralVol"], legend.title = "Mean cerebral volume \n (log10 mm3)")
graphics.off()

# Cerebello-cerebral ratio.
tiff("./7.robustness/7.10.BM-anc-cerebrocerebellar-ci.tiff", width = 8, height = 8, units = 'in', res = 300)
myplotanc.ci.Taxo(tree34S, p.BM$anc_recon[1:ntips,"CC"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"CC"], legend.title = "Cerebellar/ cerebral \n volume (%)")
graphics.off()


### -------------------------- ###
# 3. Phylogenetic regressions 
# (PGLS/ PIC Regressions)
# N. Magielse, 2022
### -------------------------- ### 

## Load packages
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(corrplot)
library(tidyverse)
library(magrittr)
library(evomap)
# library(car) #maybe not necessary (trying some things)
library(rr2) # For PGLS R2, also see R2s for Correlated Data: Phylogenetic Models, LMMs, and GLMMs from Anthony Ives.
#library(phylolm)
#library(relaimpo)

# Phenotypes are grouped by species name, as the tips of the tree (this step can be skipped here as the table contains one species mean as is.)
phenolist <- split(Stephan.log,Stephan.log$species)
phenolist <- phenolist[tree34S$tip.label]

# Get phenotypes
C <- lapply(phenolist,"[[","CerebellarVol")
V <- lapply(phenolist,"[[","CerebralVol")

#------------------------------------
# Phylogenetic independent contrasts
#------------------------------------

# Compute contrasts
# Full
pic.C <- pic(C,tree34S)
pic.V <- pic(V, tree34S)

#------------------------------------
# PIC-regression
#------------------------------------
# Compute a linear regression of C on V
result <- lm(pic.C ~ pic.V)
sink("7.robustness/7.10.lm-C-V.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on V
pdf("7.robustness/7.11.lm-C-V-mean-nointercept-plot.pdf")
plot(pic.C ~ pic.V)
abline(a = 0, b = coef(result))
graphics.off()

## PGLS under Brownian Motion (best supported model) should be the same as PIC-regression forced through 0.

## Stephan 34 species, Brownian Motion.
bm<-corBrownian(1,tree34S)
pglsModel<-gls(CerebellarVol~CerebralVol,data=Stephan.log,correlation=bm)
summary(pglsModel)
PICmodel <- lm(pic.C ~pic.V + 0)
sink("7.robustness/7.12.PGLS-PIC-C-V.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

## We also allow Pagel's lambda to vary, using different configurations of corPagel.
## We test three models, starting at lambda = 0.8, and comparing it with ANOVAs to models starting at lamda = 0.0 and 1.0, respectively.
sink("7.robustness/7.12a.PGLS-PIC-C-V-variable-lambda.txt")
PGLS.CV <-gls(CerebellarVol~CerebralVol,correlation=corPagel(value=0.8,phy=tree34S),data=Stephan.log) 
intervals(PGLS.CV,which="var-cov")
PGLS.CV.l0 <- gls(CerebellarVol~CerebralVol,correlation=corPagel(value=0.0,phy=tree34S, fixed = T),data=Stephan.log) 
PGLS.CV.l1 <- gls(CerebellarVol~CerebralVol,correlation=corPagel(value=1.0,phy=tree34S, fixed = T),data=Stephan.log) 
anova(PGLS.CV, PGLS.CV.l0)
anova(PGLS.CV, PGLS.CV.l1)
sink()

pdf("7.robustness/7.12b.PGLS-PIC-C-V-lambda-distribution.pdf")
lambda<-seq(0,1,length.out=500)
lik<-sapply(lambda,function(lambda)logLik(gls(CerebellarVol~CerebralVol, correlation=corPagel(value=lambda,phy=tree34S,fixed=TRUE), data=Stephan.log)))
plot(lik~lambda,type="l",main=expression(paste("LikelihoodPlotfor", lambda)),ylab="LogLikelihood",xlab=expression(lambda))
abline(v=PGLS.CV$modelStruct,col="black")
graphics.off()


### -------------------------- ###
### PGLS Plotting with Caper
### -------------------------- ###   

# Model 1: Cerebellar volume regressed on cerebral volume.
tiff("7.robustness/7.13.BM-PGLS-cerebellumCerebrum-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod1 <- gls(CerebellarVol ~ CerebralVol, data=Stephan.log, correlation=bm)
plot(mod1)
plot(Stephan.log$CerebellarVol ~ Stephan.log$CerebralVol, xlab = "Cerebral volume (log-transformed)", ylab = "Cerebellar volume (log-tranformed)")
abline(mod1)
abline(a=-0.4137426, b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod1.ci <-gls.ci(Stephan.log$CerebellarVol, Stephan.log$CerebralVol,vcv(tree34S))
lines(pGLS.mod1.ci$CI.plot$X,pGLS.mod1.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod1.ci$CI.plot$X,pGLS.mod1.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod1.pi <-gls.pi(Stephan.log$CerebellarVol, Stephan.log$CerebralVol,vcv(tree34S), 1)
lines(pGLS.mod1.pi$PI.plot$X,pGLS.mod1.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod1.pi$PI.plot$X,pGLS.mod1.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()


### -------------------------- ###
### END of script
### -------------------------- ### 

