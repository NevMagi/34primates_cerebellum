### -------------------------------------- ###
# Evolutionary model selection
# Roberto Toro, Katja Heuer 2018
# Script adapted by N. Magielse
# for inclusion of cerebellar data - 2022
### -------------------------------------- ###     

## SET-UP ## 
library(ape)
library(Rphylopars)
library(phytools)
library(magrittr)
library(tidyverse)

source("./input/catn.R")

# Set working directory
getwd()
setwd("~/Documents/PhD/EvolutionPrimates/Final_code/")

# Set output directory
outputdir <- "3.evo_model_selection"
dir.create(outputdir)


## LOAD DATA ##
## Load 10k tree, force ultrametric by extending branches
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
#lut <- read.csv("input/10kTrees_34PrimateSpecies_adapted.csv", sep = ";")

## Load the phenotypic data ##
# We try both raw and log-transformed)
## All measures from Heuer et al. 2019, and cerebellar and ansiform area volumes
pheno <- read.csv2("./0.cleanInput/Phenotypes_raw_outliersRM.csv", sep=';')  %>% dplyr::rename(species = WilsonReederName) %>% dplyr::select(species, CerebellarVol, SurfaceArea, CerebralVol, AbsGI, FoldingLength, FoldingNumber, Lambda, Delta, CrusVol_1_3)


## Try with only cerebellar, cerebral, and ansiform area volumes.
pheno_raw_MOI <- read.csv("./0.cleanInput/Phenotypes_raw_outliersRM.csv", sep = ";", dec = ",") %>% dplyr::rename(species = WilsonReederName) %>% dplyr::select(species, CerebellarVol, CerebralVol, CrusVol_1_3)
# pheno_MOI <- read.csv("./0.cleanInput/Phenotypes_outliersRM.csv", sep = ";", dec = ",") %>% rename(species = WilsonReederName) %>% dplyr::select(species, CerebellarVol, CerebralVol, CrusVol_1_3)

## Explore the effect of only using species median values
phenoMedian_raw_MOI <- read.csv("./0.cleanInput/PhenotypesMedian_raw_outliersRM.csv", sep = ";", dec = ",") %>% dplyr::rename(species = WilsonReederName) %>% dplyr::select(species, MedianCerebellum, MedianCerebrum, MedianCrus)
# phenoMedian_MOI <- read.csv("./0.cleanInput/PhenotypesMedian_outliersRM.csv", sep = ";", dec = ",") %>% rename(species = WilsonReederName) %>% dplyr::select(species, MedianCerebellum, MedianCerebrum, MedianCrus)

## Only measurements from Heuer et al. 2019
pheno2 <- pheno %>% select(species, SurfaceArea, CerebralVol, AbsGI, FoldingLength, FoldingNumber, Lambda, Delta)

#-------------------------------
options(width=as.integer(200))


## Fit a Brownian Motion Model for all neocortical measurements from Heuer et al. 2019, and Cerebellar + Ansiform Area volume.
p_BM <- phylopars(trait_data=pheno, tree=tree)


## TESTING BEST-SUPPORTED EVOLUTIONARY MODEL BASED ON THE TRAIT DATA ##
## Using Pagel's correction, we can find evidence on whether a model fully explained by phylogenetic signal (lambda =1) or one devoid of phylogenetic signal (star model; lambda =0 better explains the data.
## Test Pagel's lambda 1 (no transformation of branch lengths) versus 0 (complete star phylogeny, with all tip branches equal in length and all internal branches of length 0): 
## We construct phylogenetic and phenotypic trait variance-covariance matrices for both cases.
## We chose to compare models based on full phenotypic data.
catn("Fit and compare Pagel's lambda model for lambda=1 (BM) and lambda=0 (star)")
p_lambda <- phylopars(trait_data=pheno, tree=tree, usezscores = F)
p_star <- phylopars(trait_data=pheno, tree=tree, model="star", usezscores = F)
sink("3.evo_model_selection/3.1.test-pagel-lamdba=1-vs-lambda=0_full.txt")
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
p_OU <- phylopars(trait_data = pheno, tree=tree, model = "OU", usezscores = F)
sink("3.evo_model_selection/3.2.OU-univariate_full.txt")
p_OU
sink()

# Multivariate 
catn("Fit multivariate OU (a diagonal matrix of alpha values)")
p_mvOU <- phylopars(trait_data = pheno,tree = tree, model = "mvOU", full_alpha = F, usezscores = F)
sink("3.evo_model_selection/3.3.OU-multivariate_full.txt")
p_mvOU
sink()

# Diagonal Multivariate
catn("Fit multivariate OU (a diagonal matrix of alpha values)")
p_mvOU_diag <- phylopars(trait_data = pheno,tree = tree, model = "mvOU", usezscores = F)
sink("3.evo_model_selection/3.4.OU-multivariate-diagonal_full.txt")
p_mvOU_diag
sink()

## Fit the Early-Burst model ##
#------------------------------------
catn("Fit early-burst model (EB)")
p_EB <- phylopars(trait_data = pheno, tree = tree, model = "EB", usezscores = F)
sink("3.evo_model_selection/3.5.EB_full.txt")
p_EB # Estimated trait covariance and EB rate parameter
sink()

# Model selection  
#------------------------------------
catn("Model selection using AIC")
sink("3.evo_model_selection/3.6.model-selection_full.txt")
catn("Brownian motion", AIC(p_BM), sep='\t')
catn("Ornstein-Uhlenbeck, single alpha", AIC(p_OU), sep='\t')
catn("Ornstein-Uhlenbeck, diagonal alpha matrix", AIC(p_mvOU), sep='\t')
catn("Ornstein-Uhlenbeck, full alpha matrix", AIC(p_mvOU_diag), sep='\t')
catn("Early burst", AIC(p_EB), sep='\t')
catn("Star model, lambda = 0", AIC(p_star), sep='\t')
sink()


## TESTING BEST-SUPPORTED EVOLUTIONARY MODEL BASED ON ONLY CEREBELLAR, CEREBRAL, AND ANSIFORM AREA VOLUME ##
# These were ultimately used to draw the conclusions in the paper, because they were the volumes used in the primary analyses.


p_BM_a <- phylopars(trait_data=pheno_raw_MOI, tree=tree)

## Using Pagel's correction, we can find evidence on whether a model fully explained by phylogenetic signal (lambda =1) or one devoid of phylogenetic signal (star model; lambda =0 better explains the data.
## Test Pagel's lambda 1 (no transformation of branch lengths) versus 0 (complete star phylogeny, with all tip branches equal in length and all internal branches of length 0): 
## We construct phylogenetic and phenotypic trait variance-covariance matrices for both cases.
## We chose to construct both models based on phenotypic data from the current study: cerebellar, cerebral, and ansiform area volumes,# as well as ratios between the volumes #.
# ?? scales the tree between a constant-rates model (??=1) to one where every species is statistically independent of every other species in the tree (??=0).
catn("Fit and compare Pagel's lambda model for lambda=1 (BM) and lambda=0 (star)")
p_lambda_a <- phylopars(trait_data=pheno_raw_MOI, tree=tree, usezscores = F)
p_star_a <- phylopars(trait_data=pheno_raw_MOI, tree=tree, model="star", usezscores = F)
sink("3.evo_model_selection/3.1a.test-pagel-lamdba=1-vs-lambda=0_MOI.txt")
p_lambda_a
p_star_a
catn('Log-likelihood full-BM', logLik(p_lambda_a), sep='\t')
catn('Log-likelihood star',logLik(p_star_a), sep='\t')
chi_square_a <- as.double(2*(logLik(p_lambda_a) - logLik(p_star_a))) # 2*(logLik_alt - logLik_null)
catn('chi_square', chi_square_a, sep='\t')
degrees_freedom_a <- p_lambda_a$npars - p_star_a$npars # df = difference in number of model parameters
catn('df', degrees_freedom_a, sep='\t')
p_a <- pchisq(q = chi_square_a,df = degrees_freedom_a,lower.tail = FALSE) # p-value
catn ('p-value', p, sep='\t')
sink()

## Fit the Ornstein-Uhlenbeck model ##
#------------------------------------
# Univariate
# Save phylogenetic and phenotypic trait variance-covariance matrices, variance explained by phylogeny
catn("Fit Ornstein-Uhlenbeck model (OU, univariate: a single alpha value shared by all traits)")
p_OU_a <- phylopars(trait_data = pheno_raw_MOI, tree=tree, model = "OU", usezscores = F)
sink("3.evo_model_selection/3.2a.OU-univariate_MOI.txt")
p_OU_a
sink()

# Multivariate
catn("Fit multivariate OU (a diagonal matrix of alpha values)")
p_mvOU_a <- phylopars(trait_data = pheno_raw_MOI,tree = tree, model = "mvOU", full_alpha = F, usezscores = F)
sink("3.evo_model_selection/3.3a.OU-multivariate_MOI.txt")
p_mvOU_a
sink()

#Diagonal Multivariate
catn("Fit multivariate OU (a diagonal matrix of alpha values)")
p_mvOU_diag_a <- phylopars(trait_data = pheno_raw_MOI,tree = tree, model = "mvOU", usezscores = F)
sink("3.evo_model_selection/3.4a.OU-multivariate-diagonal_MOI.txt")
p_mvOU_diag_a
sink()

## Fit the Early-Burst model ##
#------------------------------------
catn("Fit early-burst model (EB)")
p_EB_a <- phylopars(trait_data = pheno_raw_MOI, tree = tree, model = "EB", usezscores = F)
sink("3.evo_model_selection/3.5a.EB_MOI.txt")
p_EB # Estimated trait covariance and EB rate parameter
sink()

# Model selection  
#------------------------------------
catn("Model selection using AIC")
sink("3.evo_model_selection/3.6a.model-selection_MOI.txt")
catn("Brownian motion", AIC(p_BM_a), sep='\t')
catn("Ornstein-Uhlenbeck, single alpha", AIC(p_OU_a), sep='\t')
catn("Ornstein-Uhlenbeck, diagonal alpha matrix", AIC(p_mvOU), sep='\t')
catn("Ornstein-Uhlenbeck, full alpha matrix", AIC(p_mvOU_diag), sep='\t')
catn("Early burst", AIC(p_EB_a), sep='\t')
catn("Star model, lambda = 0", AIC(p_star_a), sep='\t')
sink()


## Only Heuer et al. 2019 measurements ##

p_BM_b <- phylopars(trait_data=pheno2, tree=tree)

## Using Pagel's correction, we can find evidence on whether a model fully explained by phylogenetic signal (lambda =1) or one devoid of phylogenetic signal (star model; lambda =0 better explains the data.
## Test Pagel's lambda 1 (no transformation of branch lengths) versus 0 (complete star phylogeny, with all tip branches equal in length and all internal branches of length 0): 
## We construct phylogenetic and phenotypic trait variance-covariance matrices for both cases.
# ?? scales the tree between a constant-rates model (??=1) to one where every species is statistically independent of every other species in the tree (??=0).
catn("Fit and compare Pagel's lambda model for lambda=1 (BM) and lambda=0 (star)")
p_lambda_b <- phylopars(trait_data=pheno2, tree=tree, usezscores = F)
p_star_b <- phylopars(trait_data=pheno2, tree=tree, model="star", usezscores = F)
sink("3.evo_model_selection/3.1b.test-pagel-lamdba=1-vs-lambda=0_old.txt")
p_lambda_b
p_star_b
catn('Log-likelihood full-BM', logLik(p_lambda_b), sep='\t')
catn('Log-likelihood star',logLik(p_star_b), sep='\t')
chi_square_b <- as.double(2*(logLik(p_lambda_b) - logLik(p_star_b))) # 2*(logLik_alt - logLik_null)
catn('chi_square', chi_square_b, sep='\t')
degrees_freedom_b <- p_lambda_b$npars - p_star_b$npars # df = difference in number of model parameters
catn('df', degrees_freedom_b, sep='\t')
p_b <- pchisq(q = chi_square_b,df = degrees_freedom_b,lower.tail = FALSE) # p-value
catn ('p-value', p, sep='\t')
sink()

## Fit the Ornstein-Uhlenbeck model ##
#------------------------------------
# Univariate
# Save phylogenetic and phenotypic trait variance-covariance matrices, variance explained by phylogeny
catn("Fit Ornstein-Uhlenbeck model (OU, univariate: a single alpha value shared by all traits)")
p_OU_b <- phylopars(trait_data = pheno2, tree=tree, model = "OU", usezscores = F)
sink("3.evo_model_selection/3.2b.OU-univariate_old.txt")
p_OU_b
sink()

# Multivariate
catn("Fit multivariate OU (a diagonal matrix of alpha values)")
p_mvOU_b <- phylopars(trait_data = pheno,tree = tree, model = "mvOU", full_alpha = F, usezscores = F)
sink("3.evo_model_selection/3.3b.OU-multivariate_old.txt")
p_mvOU_b
sink()

#Diagonal Multivariate
catn("Fit multivariate OU (a diagonal matrix of alpha values)")
p_mvOU_diag_b <- phylopars(trait_data = pheno,tree = tree, model = "mvOU", usezscores = F)
sink("3.evo_model_selection/3.4b.OU-multivariate-diagonal_old.txt")
p_mvOU_diag_b
sink()

## Fit the Early-Burst model ##
#------------------------------------
catn("Fit early-burst model (EB)")
p_EB_b <- phylopars(trait_data = pheno, tree = tree, model = "EB", usezscores = F)
sink("3.evo_model_selection/3.5b.EB_old.txt")
p_EB_b # Estimated trait covariance and EB rate parameter
sink()

# Model selection  
#------------------------------------
catn("Model selection using AIC")
sink("3.evo_model_selection/3.6b.model-selection_old.txt")
catn("Brownian motion", AIC(p_BM_b), sep='\t')
catn("Ornstein-Uhlenbeck, single alpha", AIC(p_OU_b), sep='\t')
catn("Ornstein-Uhlenbeck, diagonal alpha matrix", AIC(p_mvOU_b), sep='\t')
catn("Ornstein-Uhlenbeck, full alpha matrix", AIC(p_mvOU_diag_b), sep='\t')
catn("Early burst", AIC(p_EB_b), sep='\t')
catn("Star model, lambda = 0", AIC(p_star_b), sep='\t')
sink()


### -------------------------- ###
### END of script
### -------------------------- ###   
