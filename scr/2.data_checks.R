### ---------------------------------- ###
# Two important additional data checks
# N. Magielse, 2022
### ---------------------------------- ###     

# We are using R version 4.1.0 (2021-05-18)

### ---------------------------------------------------------------------------- ###
# PART I:
# Testing parametric assumptions about the data by checking if the residual error, 
# corrected for main effects as well as the phylogenetic tree structure, is normal.
### ---------------------------------------------------------------------------- ###     

### LOAD PACKAGES ### 
library(phytools)
library(nlme)
library(nortest)
library(tidyverse)
library(ggplot2)

## Set directories ##
setwd("~/Documents/PhD/EvolutionPrimates/Final_code/")
outputdir <- "2.data_checks"
dir.create(outputdir)

## Load data ##
tree34 <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
tree13 <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_13.nex"), method="extend")
phenoMedianOutliers <- read.csv2("./0.cleanInput/PhenotypesMedian_outliersRM.csv")
source('./input/catn.R')

## Make vectors with dependent and independent variables, as well as residual error based on the tree structure ##

## First specifify whic specimens have ansiform area observations.
ansiform <- phenoMedianOutliers$MedianCrus
names(ansiform) <- phenoMedianOutliers$WilsonReederName
ansiform <- ansiform[!is.na(ansiform)]

# Cerebellum, cerebrum, and specimens with complete observations (13 species; include ansiform area.)
cerebellum <- phenoMedianOutliers$MedianCerebellum
names(cerebellum) <- phenoMedianOutliers$WilsonReederName
cerebellum13 <- cerebellum %>%
  subset(names(cerebellum) %in% names(ansiform))
cerebrum <-  phenoMedianOutliers$MedianCerebrum
names(cerebrum) <- phenoMedianOutliers$WilsonReederName
cerebrum13 <- cerebrum %>%
  subset(names(cerebrum) %in% names(ansiform))
pheno13 <- data.frame(cerebellum13, cerebrum13, ansiform) 

# Calculate residual errors based on the best-known tree structure.
ResidualError34 = fastBM(tree34)
ResidualError13 = fastBM(tree13)

# We used the Lilliefors test, since variance and mean of the distribution are unknown.
# Lilliefors, Hubert W. (1967-06-01). "On the Kolmogorov-Smirnov Test for Normality with Mean and Variance Unknown". Journal of the American Statistical Association. 62 (318): 399–402. doi:10.1080/01621459.1967.10482916. ISSN 0162-1459
# The raw data need not be distributed normally. We care instead about multivariate normality of the residual error in y given X (http://blog.phytools.org/2013/02/a-comment-on-distribution-of-residuals.html).
lillie_cerebellum <- lillie.test(cerebellum) #NORMAL
lillie_cerebrum <- lillie.test(cerebrum) #NOT NORMAL
lillie_ansiform <- lillie.test(ansiform) #NORMAL

## We thus test whether the residual errors in the different comparisons are normally distributed. 
## We do this both for raw (read: not phylogenetically-corrected) phenotypes and phylogenetic independent contrasts. Phylogenetic independent contrasts are made both from species median data, and all specimens (excl. outliers).
## Since (3.evo_model_selection) points out that a Brownian Motion is best supported (significantly more so than a star-model) we draw out conclusions based on the PIC-models.
## This is important, as phylogenetic and non-phylogenetic models make opposite assumptions about the data's correlation structure.

## Non-phylogenetic data
# Cerebro-cerebellar
nonPhylo_cerebrocerebellar <- lm(MedianCerebellum ~ MedianCerebrum, data = phenoMedianOutliers)
nonPhylo_cerebrocerebellar
lillie_nonPhyloCC <- lillie.test(residuals(nonPhylo_cerebrocerebellar)) #NORMAL

# Ansiform Area as a factor of cerebellar volume
nonPhylo_aac <- lm(MedianCrus ~ MedianCerebellum, data = phenoMedianOutliers)
nonPhylo_aac
lillie_nonPhyloAAC <- lillie.test(residuals(nonPhylo_aac)) #NORMAL

# Ansiform Area as a factor of cerebral volume
nonPhylo_aacerebrum <- lm(MedianCrus ~ MedianCerebrum, data = phenoMedianOutliers)
nonPhylo_aacerebrum
lillie_nonPhyloAACerebrum <- lillie.test(residuals(nonPhylo_aacerebrum)) #NORMAL


## Phylogenetically-corrected data 
# Cerebro-cerebellar
fitCC<-gls(MedianCerebellum~MedianCerebrum, data = phenoMedianOutliers, correlation=corBrownian(1,tree34)) ## A Brownian Motion model of evolution was best supported (see 3.evo_model_selection)
fitCC
lillie_fitCC <- lillie.test(residuals(fitCC)) #NORMAL

# Ansiform Area as a factor of cerebellar volume
fitAAC<-gls(ansiform~cerebellum13, data = pheno13, correlation=corBrownian(1,tree13)) ## A Brownian Motion model of evolution was best supported (see 3.evo_model_selection)
fitAAC
lillie_fitAAC <- lillie.test(residuals(fitAAC)) #NORMAL

# Ansiform Area as a factor of cerebral volume
fitAACerebrum<-gls(ansiform~cerebrum13, data = pheno13, correlation=corBrownian(1,tree13)) ## A Brownian Motion model of evolution was best supported (see 3.evo_model_selection)
fitAACerebrum
lillie_fitAACerebrum<- lillie.test(residuals(fitAACerebrum)) #NORMAL

## Save the parametric tests as csv ##
catn("Test of parametric assumptions")
sink("2.data_checks/2.1.parametric_raw_phenotypes.txt")
lillie_cerebellum
lillie_cerebrum
lillie_ansiform
sink()

sink("2.data_checks/2.2.parametric_nonPhylo_ols.txt")
lillie_nonPhyloCC
lillie_nonPhyloAAC
lillie_nonPhyloAACerebrum
sink()

sink("2.data_checks/2.3.parametric_phylo_gls.txt")
lillie_fitCC
lillie_fitAAC
lillie_fitAACerebrum
sink()

### --------------------------------------------------------------- ###
# PART II:
# Testing inter observer reliability of Ansiform Area segmentations. 
# Segmentations were made by N. Magielse (obs1)  & M. Abbaspour (obs2)
### --------------------------------------------------------------- ### 

## SET-UP ##
library(nlme)
library(phytools)
#library(geiger)
library(ape)
library(tidyverse)

### DATA ###
interobserver <- read_csv2("./0.cleanInput/InterObserver_reliability.csv") %>% gather(Observer,Ansiform_Area_Volume, CrusVol_obs1:CrusVol_obs2) %>% dplyr::select(c(WilsonReederName, Observer, Ansiform_Area_Volume))
interobserver[order(interobserver$WilsonReederName), ] # Measurements for the 2 different observers for 6 species

## GENERAL TWO-WAY ANOVA ##
sink("2.data_checks/2.4.interobserver_tw_anova.txt")
interobs.aov <- aov(Ansiform_Area_Volume ~ WilsonReederName + Observer, data = interobserver)
summary(interobs.aov)
sink()

## Highly significant effect of species (as expected). No effect of observer, good.

## Also calculate intraclass correlations (ICCs)
library(irr)
interobserverICC <- read_csv2("./0.cleanInput/InterObserver_reliability.csv") %>% dplyr::select(CrusVol_obs1, CrusVol_obs2)
interobserverICC$CrusVol_obs1 <- 10^interobserverICC$CrusVol_obs1
interobserverICC$CrusVol_obs2 <- 10^interobserverICC$CrusVol_obs2
sink("2.data_checks/2.5.interobserver_ICC.txt")
icc(interobserverICC,model='twoway',type='agreement')
sink()

### -------------------------------------------------------- ###
# PART III:
# Testing manual cerebellar and ansiform area segmentations 
# versus CERES automated segmentations (for humans only). 
# Segmentations were made by N. Magielse  & V. Steigauf (obs3).
### -------------------------------------------------------- ### 

manual.v.auto <- read_csv2("./input/manual_automated_segmentations.csv")
whole <- manual.v.auto %>% dplyr::select(c(SpecimenID, Manual_full, Automated_full)) %>% gather(Method,Cerebellar_volume, Manual_full:Automated_full) 
ansiform <- manual.v.auto %>% dplyr::select(c(SpecimenID, Manual_crus, Automated_crus)) %>% gather(Method,Ansiform_volume, Manual_crus:Automated_crus)

# Whole cerebellar volume
sink("2.data_checks/2.6.human_manual-v-automated_cerebellum_tw_anova.txt")
whole.aov <- aov(Cerebellar_volume ~  Method, data = whole)
summary(whole.aov)
sink()

# Ansiform area volume
sink("2.data_checks/2.7.human_manual-v-automated_ansiform_tw_anova.txt")
ansiform.aov <- aov(Ansiform_volume ~  Method, data = ansiform)
summary(ansiform.aov)
sink()

# It seems that mostly the means differ.
# Just to check if segmentations otherwise resemble each other, 
# we equalize the means and perform anova, too.
colMeans(manual.v.auto[,2:5])
(146.2790-125.8500)/146.2790*100
# = 13.96578% difference
(43.9377-41.2130)/43.9377*100
# = 6.20128% difference

manual.v.auto2 <- manual.v.auto 
manual.v.auto2$Automated_full<-1.1396578*(manual.v.auto2$Automated_full)
manual.v.auto2$Automated_crus<-1.0620128*(manual.v.auto2$Automated_crus)

sink("2.data_checks/2.8.meansEqualized_manual-v-automated_tw_anova.txt")
print("Do not use, nor trust, this data!")
whole2 <- manual.v.auto2 %>% dplyr::select(c(SpecimenID, Manual_full, Automated_full)) %>% gather(Method,Cerebellar_volume, Manual_full:Automated_full) 
ansiform2 <- manual.v.auto2 %>% dplyr::select(c(SpecimenID, Manual_crus, Automated_crus)) %>% gather(Method,Ansiform_volume, Manual_crus:Automated_crus)
whole.aov2 <- aov(Cerebellar_volume ~  Method + SpecimenID, data = whole2)
summary(whole.aov2)
ansiform.aov2 <- aov(Ansiform_volume ~  Method, data = ansiform2)
summary(ansiform.aov2)
sink()

## Calculate manual vs. automated ICCs for whole cerebellar volume and ansiform area volume.
cerebellumICC.MA <- manual.v.auto %>% dplyr::select(Manual_full, Automated_full)
ansiformICC.MA <- manual.v.auto %>% dplyr::select(Manual_crus, Automated_crus)

library(irr)
sink("2.data_checks/2.9.Manual-v-Automated_ICC.txt")
icc(cerebellumICC.MA,model='oneway',type='agreement')
icc(ansiformICC.MA,model='oneway',type='agreement')
sink()


### -------------------------- ###
### END of script
### -------------------------- ###   
