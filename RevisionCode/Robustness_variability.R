############################################
### Additional analysis taking 100 random
### multiplications of the volumetric data,
### establishing a measure of robustness of
### our main conclusions in light of the 
### great variability that might be possible.
### This volumetric variability might come from
### unaccounted for sex, ages, (tissue-specific
### shrinkage, and most importantly, general
### intraspecific variation.
###
### We assume that this variation may maximally
### cause a 2.0x (or 0.5x) miss-estimation of the 
### actual volumes, so transformations of the matrix
### will be randomly sampled between these values.
############################################

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
setwd("~/Documents/PhD/EvolutionPrimates/RevisionCode") ## Provide path to folder with scripts and input data.

## Set output directory
outputdir <- "3.Robustness_variability"
dir.create(outputdir)

## Load data
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
tree13 <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_13.nex"), method="extend") # To be used for those species with full values (no NAs for Ansiform Area Volume)
pheno <- read.csv2("./input/Phenotypes_outliersRM.csv", sep=';')  %>% dplyr::rename(species = WilsonReederName) %>% dplyr::select(species, Hominoidea, CerebellarVol, SurfaceArea, CerebralVol, AbsGI, FoldingLength, FoldingNumber, Lambda, Delta, CrusVol_1_3, ROC)
phenoAA <- pheno %>% drop_na()
phenoMedian <- read.csv2("./input/PhenotypesMedian_outliersRM.csv", sep=';')  %>% dplyr::rename(species = WilsonReederName) %>% dplyr::select(species, Hominoidea, MedianCerebellum, MedianCerebrum, MedianCrus, MedianROC)
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
V <- lapply(phenolist,"[[","CerebralVol")
AA <- lapply(phenolist,"[[","CrusVol_1_3")
ROC <- lapply(phenolist,"[[","ROC")

# Get phenotypes for species without NAs anywhere.
C13 <- lapply(phenolist13,"[[","CerebellarVol")
V13 <- lapply(phenolist13,"[[","CerebralVol")
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


#------------------------------------
# Phylogenetic independent contrasts
#------------------------------------

# Compute contrasts
pic.C <- pic.ortho(C,tree,intra=TRUE)
pic.V <- pic.ortho(V,tree,intra=TRUE)
pic.AA <- pic.ortho(AA, tree, intra=TRUE) ## This one does not really make sense. We are only able to calculate PIC for nodes 51 and 52. PICs are more sensibly calculated from the 13-species tree, that has complete data.
pic.ROC <- pic.ortho(ROC, tree, intra=TRUE) ## Ibidem.

# Compute contrasts for 13-species data
pic.C13 <- pic.ortho(C13,tree13,intra=TRUE)
pic.V13 <- pic.ortho(V13,tree13,intra=TRUE)
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

### Saving PGLS and PIC coefficients as separate columns will facilitate downstream analyses.

################################################
### We now want to make 1000 random transformations
### on the original phenotype dataframe (phenoMedian)
### and specifically its MedianCerebellum and Median-
### Cerebrum columns. These random transformations 
### are constraint to range between 0.5 and 2.0 based
### on our assumption on maximum error in volumes
### due to unaccounted for (unaccountable-for) natural
### intraspecific variability especially (also differ-
### ences in data quality, in shrinkage). 
###
################################################

### First, we examine the cerebello-cerebral PGLS regression.

########################################  
## MAIN model. This simulation introduces
## only ratio differences between cerebellar
## and cerebral volumes.
########################################  

set.seed(123) # for reproducibility

# Create a function that performs the transformations and regressions
pgls_regression <- function(df, tree) {
  
  # Generate a random ratio between 0.91 and 1.1
  ratio <- runif(1, min = 0.91, max = 1.1)
  
  # Apply the random ratio transformation to MedianCerebellum and MedianCerebrum
  df$MedianCerebellum <- df$MedianCerebellum * ratio
  df$MedianCerebrum <- df$MedianCerebrum / ratio
  
  # Fit PGLS model under Brownian Motion
  bm <- corBrownian(1, tree)
  pgls_model <- gls(MedianCerebellum ~ MedianCerebrum, data = df, correlation = bm)
  
  # Extract coefficients from the PGLS model
  pgls_intercept <- pgls_model$coefficients[[1]]
  pgls_med_cerebrum <- pgls_model$coefficients[["MedianCerebrum"]]
  
  # Calculate R-squared for PGLS model
  pgls_residuals <- residuals(pgls_model)
  pgls_tss <- sum((df$MedianCerebellum - mean(df$MedianCerebellum))^2)
  pgls_rss <- sum(pgls_residuals^2)
  pgls_R_squared <- 1 - (pgls_rss / pgls_tss)
  
  # Return the results as a list
  return(list(pgls_intercept = pgls_intercept,
              pgls_med_cerebrum = pgls_med_cerebrum,
              pgls_R_squared = pgls_R_squared))
}


# Perform the transformations and regressions for 10000 randomly transformed dataframes
results <- lapply(1:10000, function(i) pgls_regression(phenoMedian, tree))

# Combine the results into a data frame
results_df <- do.call(rbind, results)

# save the results to a CSV file
write.csv(results_df, "3.Robustness_variability/3.1.PGLS-C-V-10000.csv", row.names = FALSE)


########################################  
## This model also introduces up to 2x 
## variability in between rows of the dataframe,
##.not just ratio changes between the columns).
########################################  
set.seed(123) # for reproducibility

# Create a function that performs the transformations and regressions
pgls_regression_intraspecific_conserv <- function(df, tree) {
  
  # Apply random transformations to MedianCerebellum and MedianCerebrum
  ratio <- runif(nrow(df), min = 0.91, max = 1.1)
  transform_factor <- runif(nrow(df), min = 0.5, max = 2)
  df$MedianCerebellum <- df$MedianCerebellum * transform_factor * ratio
  df$MedianCerebrum <- df$MedianCerebrum * transform_factor / ratio
  
  # Fit PGLS model under Brownian Motion
  bm <- corBrownian(1, tree)
  pgls_model <- gls(MedianCerebellum ~ MedianCerebrum, data = df, correlation = bm)
  
  # Extract coefficients from the PGLS model
  pgls_intercept <- pgls_model$coefficients[[1]]
  pgls_med_cerebrum <- pgls_model$coefficients[["MedianCerebrum"]]
  
  # Calculate R-squared for PGLS model
  pgls_residuals <- residuals(pgls_model)
  pgls_tss <- sum((df$MedianCerebellum - mean(df$MedianCerebellum))^2)
  pgls_rss <- sum(pgls_residuals^2)
  pgls_R_squared <- 1 - (pgls_rss / pgls_tss)
  
  # Return the results as a list
  return(list(pgls_intercept = pgls_intercept,
              pgls_med_cerebrum = pgls_med_cerebrum,
              pgls_R_squared = pgls_R_squared))
}

# Perform the transformations and regressions for 10000 randomly transformed dataframes
results_intraspecific_conservative <- lapply(1:10000, function(i) pgls_regression_intraspecific_conserv(phenoMedian, tree))

# Combine the results into a data frame
results_df_intraspecific_conservative <- do.call(rbind, results_intraspecific_conservative)

# save the results to a CSV file
write.csv(results_df_intraspecific_conservative, "3.Robustness_variability/3.2.PGLS-C-V-10000-intraspecific-conservative.csv", row.names = FALSE)


# This version assumes a smaller error; a factor between 2/3 and 1.5 for unaccounted for intraspecefic variation.
pgls_regression_intraspecific <- function(df, tree) {
  
  # Apply random transformations to MedianCerebellum and MedianCerebrum
  ratio <- runif(nrow(df), min = 0.91, max = 1.1)
  transform_factor <- runif(nrow(df), min = 2/3, max = 1.5)
  df$MedianCerebellum <- df$MedianCerebellum * transform_factor * ratio
  df$MedianCerebrum <- df$MedianCerebrum * transform_factor / ratio
  
  # Fit PGLS model under Brownian Motion
  bm <- corBrownian(1, tree)
  pgls_model <- gls(MedianCerebellum ~ MedianCerebrum, data = df, correlation = bm)
  
  # Extract coefficients from the PGLS model
  pgls_intercept <- pgls_model$coefficients[[1]]
  pgls_med_cerebrum <- pgls_model$coefficients[["MedianCerebrum"]]
  
  # Calculate R-squared for PGLS model
  pgls_residuals <- residuals(pgls_model)
  pgls_tss <- sum((df$MedianCerebellum - mean(df$MedianCerebellum))^2)
  pgls_rss <- sum(pgls_residuals^2)
  pgls_R_squared <- 1 - (pgls_rss / pgls_tss)
  
  # Return the results as a list
  return(list(pgls_intercept = pgls_intercept,
              pgls_med_cerebrum = pgls_med_cerebrum,
              pgls_R_squared = pgls_R_squared))
}

# Perform the transformations and regressions for 10000 randomly transformed dataframes
results_intraspecific <- lapply(1:10000, function(i) pgls_regression_intraspecific(phenoMedian, tree))

# Combine the results into a data frame
results_df_intraspecific <- do.call(rbind, results_intraspecific)

# save the results to a CSV file
write.csv(results_intraspecific, "3.Robustness_variability/3.3.PGLS-C-V-10000-intraspecific.csv", row.names = FALSE)


###########################################################
#### Ansiform area regression on rest of cerebellar volume
###########################################################

set.seed(123) # for reproducibility

# Create a function that performs the transformations and regressions
pgls_regression_conserv_AAROC <- function(df, tree) {
  
  # Generate a random ratio between 0.91 and 1.1
  ratio <- runif(1, min = 0.91, max = 1.1)
  
  # Apply random transformations to MedianCrus and MedianROC
  df$MedianCrus <- df$MedianCrus * ratio
  df$MedianROC <- df$MedianROC  / ratio
  
  # Fit PGLS model under Brownian Motion
  bm <- corBrownian(1, tree13)
  pgls_model <- gls(MedianCrus ~ MedianROC, data = df, correlation = bm)
  
  # Extract coefficients from the PGLS model
  pgls_intercept <- pgls_model$coefficients[[1]]
  pgls_med_ROC <- pgls_model$coefficients[["MedianROC"]]
  
  # Calculate R-squared for PGLS model
  pgls_residuals <- residuals(pgls_model)
  pgls_tss <- sum((df$MedianCrus - mean(df$MedianCrus))^2)
  pgls_rss <- sum(pgls_residuals^2)
  pgls_R_squared <- 1 - (pgls_rss / pgls_tss)
  
  # Return the results as a list
  return(list(pgls_intercept = pgls_intercept,
              pgls_med_ROC = pgls_med_ROC,
              pgls_R_squared = pgls_R_squared))
}

# Perform the transformations and regressions for 10000 randomly transformed dataframes
results_AAROC <- lapply(1:10000, function(i) pgls_regression_conserv_AAROC(phenoMedianAA, tree13))

# Combine the results into a data frame
results_df_AAROC <- do.call(rbind, results_AAROC)

# save the results to a CSV file
write.csv(results_df_AAROC, "3.Robustness_variability/3.4.PGLS-AA-ROC-10000.csv", row.names = FALSE)

###########################################################
#### Ansiform area regression on cerebral volume
###########################################################

set.seed(123) # for reproducibility

# Create a function that performs the transformations and regressions
pgls_regression_conserv_AAV <- function(df, tree13) {
  
  # Apply random transformations to MedianCrus and MedianROC
  ratio <- runif(1, min = 0.91, max = 1.1)
  df$MedianCrus <- df$MedianCrus * ratio
  df$MedianCerebrum <- df$MedianCerebrum / ratio
  
  # Fit PGLS model under Brownian Motion
  bm <- corBrownian(1, tree13)
  pgls_model <- gls(MedianCrus ~ MedianCerebrum, data = df, correlation = bm)
  
  # Extract coefficients from the PGLS model
  pgls_intercept <- pgls_model$coefficients[[1]]
  pgls_med_Cer <- pgls_model$coefficients[["MedianCerebrum"]]
  
  # Calculate R-squared for PGLS model
  pgls_residuals <- residuals(pgls_model)
  pgls_tss <- sum((df$MedianCrus - mean(df$MedianCrus))^2)
  pgls_rss <- sum(pgls_residuals^2)
  pgls_R_squared <- 1 - (pgls_rss / pgls_tss)
  
  # Return the results as a list
  return(list(pgls_intercept = pgls_intercept,
              pgls_med_Cer = pgls_med_Cer,
              pgls_R_squared = pgls_R_squared))
}

# Perform the transformations and regressions for 10000 randomly transformed dataframes
results_AAV <- lapply(1:10000, function(i) pgls_regression_conserv_AAV(phenoMedianAA, tree13))

# Combine the results into a data frame
results_df_AAV <- do.call(rbind, results_AAV)

# save the results to a CSV file
write.csv(results_df_AAV, "3.Robustness_variability/3.5.PGLS-AA-V-10000.csv", row.names = FALSE)


################################################
### We now visualize how the intercepts and the 
### slopes of the PGLS analyses are divided over
### 10.000 simulations with random transformations.
### 
### We also illustrate how our original PGLS analysis
### maps on these distributions.
################################################
library(ggplot2)

## Cerebellum vs cerebrum, only cross-tissue differences.
results_df <- results_df %>% as.data.frame() 
results_df$pgls_intercept <- as.numeric(results_df$pgls_intercept)
results_df$pgls_med_cerebrum <- as.numeric(results_df$pgls_med_cerebrum)
results_df$pgls_R_squared <- as.numeric(results_df$pgls_R_squared)

## Intercept
pdf("3.Robustness_variability/3.6.PGLS.intercept.distribution10000.CC.cross-tissue.pdf")
ggplot(results_df, aes(pgls_intercept)) +
  geom_histogram(color = "#000000", fill = "#0099F8") +
  theme_classic() +
  geom_vline(xintercept = -0.627)
graphics.off()

## Slope
pdf("3.Robustness_variability/3.6a.PGLS.slope.distribution10000.CC.cross-tissue.pdf")
ggplot(results_df, aes(pgls_med_cerebrum)) +
  geom_histogram(color = "#000000", fill = "#0099F8") +
  theme_classic() +
  geom_vline(xintercept = 0.955)
graphics.off()

## Cerebellum vs cerebrum, conservative; introducing 2x intraspecific variation.
results_df_intraspecific_conservative <- results_df_intraspecific_conservative %>% as.data.frame() 
results_df_intraspecific_conservative$pgls_intercept <- as.numeric(results_df_intraspecific_conservative$pgls_intercept)
results_df_intraspecific_conservative$pgls_med_cerebrum <- as.numeric(results_df_intraspecific_conservative$pgls_med_cerebrum)
results_df_intraspecific_conservative$pgls_R_squared <- as.numeric(results_df_intraspecific_conservative$pgls_R_squared)

## Intercept
pdf("3.Robustness_variability/3.7.PGLS.intercept.distribution10000.CC.conservative.2xintraspecific.pdf")
ggplot(results_df_intraspecific_conservative, aes(pgls_intercept)) +
  geom_histogram(color = "#000000", fill = "#0099F8") +
  theme_classic() +
  geom_vline(xintercept = -0.627)
graphics.off()

## Slope
pdf("3.Robustness_variability/3.7a.PGLS.slope.distribution10000.CC.conservative.2xintraspecific.pdf")
ggplot(results_df_intraspecific_conservative, aes(pgls_med_cerebrum)) +
  geom_histogram(color = "#000000", fill = "#0099F8") +
  theme_classic() +
  geom_vline(xintercept = 0.955)
graphics.off()

## Cerebellum vs cerebrum; introducing 1.5x intraspecific variation.
results_df_intraspecific <- results_df_intraspecific %>% as.data.frame() 
results_df_intraspecific$pgls_intercept <- as.numeric(results_df_intraspecific$pgls_intercept)
results_df_intraspecific$pgls_med_cerebrum <- as.numeric(results_df_intraspecific$pgls_med_cerebrum)
results_df_intraspecific$pgls_R_squared <- as.numeric(results_df_intraspecific$pgls_R_squared)

## Intercept
pdf("3.Robustness_variability/3.8.PGLS.intercept.distribution10000.CC.1.5xintraspecific.pdf")
ggplot(results_df_intraspecific, aes(pgls_intercept)) +
  geom_histogram(color = "#000000", fill = "#0099F8") +
  theme_classic() +
  geom_vline(xintercept = -0.627)
graphics.off()

## Slope
pdf("3.Robustness_variability/3.8a.PGLS.slope.distribution10000.CC.1.5xintraspecific.pdf")
ggplot(results_df_intraspecific, aes(pgls_med_cerebrum)) +
  geom_histogram(color = "#000000", fill = "#0099F8") +
  theme_classic() +
  geom_vline(xintercept = 0.955)
graphics.off()


## Ansiform area vs rest of cerebellum, conservative
results_df_AAROC <- results_df_AAROC %>% as.data.frame() 
results_df_AAROC$pgls_intercept <- as.numeric(results_df_AAROC$pgls_intercept)
results_df_AAROC$pgls_med_ROC <- as.numeric(results_df_AAROC$pgls_med_ROC)
results_df_AAROC$pgls_R_squared <- as.numeric(results_df_AAROC$pgls_R_squared)


## Intercept
pdf("3.Robustness_variability/3.9.PGLS.intercept.distribution10000.AAROC.cross-tissue.pdf")
ggplot(results_df_AAROC, aes(pgls_intercept)) +
  geom_histogram(color = "#000000", fill = "#0099F8") +
  theme_classic() +
  geom_vline(xintercept = -1.833)
graphics.off()

## Slope
pdf("3.Robustness_variability/3.9a.PGLS.slope.distribution10000.AAROC.cross-tissue.pdf")
ggplot(results_df_AAROC, aes(pgls_med_ROC)) +
  geom_histogram(color = "#000000", fill = "#0099F8") +
  theme_classic() +
  geom_vline(xintercept = 1.297)
graphics.off()

## Ansiform area vs cerebrum, conservative
results_df_AAV <- results_df_AAV %>% as.data.frame() 
results_df_AAV$pgls_intercept <- as.numeric(results_df_AAV$pgls_intercept)
results_df_AAV$pgls_med_Cer <- as.numeric(results_df_AAV$pgls_med_Cer)
results_df_AAV$pgls_R_squared <- as.numeric(results_df_AAV$pgls_R_squared)

## Intercept
pdf("3.Robustness_variability/3.10.PGLS.intercept.distribution10000.AAV.cross-tissue.pdf")
ggplot(results_df_AAV, aes(pgls_intercept)) +
  geom_histogram(color = "#000000", fill = "#0099F8") +
  theme_classic() +
  geom_vline(xintercept = -2.784)

## Slope
pdf("3.Robustness_variability/3.10a.PGLS.slope.distribution10000.AAV.cross-tissue.pdf")
ggplot(results_df_AAV, aes(pgls_med_Cer)) +
  geom_histogram(color = "#000000", fill = "#0099F8") +
  theme_classic() +
  geom_vline(xintercept = 1.245)
graphics.off()

#################################
#### END OF SCRIPT
#################################
