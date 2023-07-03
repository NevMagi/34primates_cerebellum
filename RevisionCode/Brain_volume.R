############################################
### Analysis to examine scaling relative
### to whole brain volume. This might provide
### even better understanding of the data.
### 
### There are 7 main analyses:
### 1. Q: Check outliers in BrainVol. Also see how constant ratios (vs. BrainVol) are intraspecifically. A: No outliers beyond the originally found ones.
### 2. Q: Cerebellar ~ Cerebral partial regression. A: no notable imrpovement in fit of the PGLS model.
### 3. Cerebellar &  Cerebral ~ whole brain scaling differences for strepsirrhines/haplorhines (reviewer 3).
### 4. Perform ACE for cerebellum normalized against BrainVol.
### 5. See if we can replicate encephalization in our data - check if different in strepsirrhines/haplorrhines.
############################################

library(ape)
#library(picante)
#library(caper)
library(geiger)
library(Rphylopars)
library(phytools)
library(tidyverse)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(plyr)
library(dispRity)
library(nlme)
library(rr2)
library(evomap)
library(ggrepel)
library(ggpubr)
library(viridis)
library(hrbrthemes)
hrbrthemes::import_roboto_condensed()

## Set working directory
getwd()
setwd("~/Documents/PhD/EvolutionPrimates/RevisionCode") ## Provide path to folder with scripts and input data.

## Set output directory
outputdir <- "6.Brain_volume"
dir.create(outputdir)

## Read data
specimens <- read.csv2("./input/Phenotypes_split_BrainVolume.csv", sep=';')
specimens <- specimens[!(specimens$SpecimenID=="MacaqueCrabierMa_8f9e" | specimens$SpecimenID=="MacaqueRhesusMac_bdf8"),]
species <- read.csv2("./input/PhenotypesMedian_split_BrainVolume.csv", sep=';')
#lut <- read_csv2("./input/10kTrees_34PrimateSpecies_adapted.csv")  #look-up table
#english_names <- lut$`English Name`[match(species$WilsonReederName, lut$WilsonReederName)]
#species$EnglishName <- english_names
#species <- merge(species, lut[c("WilsonReederName", "English Name")], by = "WilsonReederName", all.x = T)
species2 <- read.csv2("./input/PhenotypesMedian_split_BrainVolume_EQplot.csv", sep=';')
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
tree13 <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_13.nex"), method="extend") # To be used for those species with full values (no NAs for Ansiform Area Volume)

##################################
### SOME PREPARATION
##################################

# Some back-and-forth, normalize cerebellar and cerebral volumes
species$MedianCerebellum <- 10^(species$MedianCerebellum) # Unlog
species$MedianCerebrum <- 10^(species$MedianCerebrum) # Unlog
species$MedianCrus <- 10^(species$MedianCrus) # Unlog
species$MedianROC <- 10^(species$MedianROC) # Unlog
species$BrainVol <- 10^(species$BrainVol) # Unlog
species$MedianCerebellumNormalized <- species$MedianCerebellum/species$BrainVol # New
species$MedianCerebrumNormalized <- species$MedianCerebrum/species$BrainVol # New
species$MedianCrusNormalized <- species$MedianCrus/species$BrainVol # New
species$MedianROCNormalized <- species$MedianROC/species$BrainVol # New
species$MedianCerebellum <- log10(species$MedianCerebellum) # Reassign
species$MedianCerebrum <- log10(species$MedianCerebrum) # Reassign
species$BrainVol <- log10(species$BrainVol) # Reassign
species$MedianCrus <- log10(species$MedianCrus) # Reassign
species$MedianROC <- log10(species$MedianROC) # Reassign

# Also for the per-specimen data
specimens$CerebellarVol <- 10^(specimens$CerebellarVol) # Unlog
specimens$CerebralVol <- 10^(specimens$CerebralVol) # Unlog
specimens$CrusVol_1_3 <- 10^(specimens$CrusVol_1_3) # Unlog
specimens$ROC <- 10^(specimens$ROC) # Unlog
specimens$BrainVol <- 10^(specimens$BrainVol) # Unlog
specimens$CerebellumNormalized <- specimens$CerebellarVol/specimens$BrainVol # New
specimens$CerebrumNormalized <- specimens$CerebralVol/specimens$BrainVol # New
specimens$CrusNormalized <- specimens$CrusVol_1_3/specimens$BrainVol # New
specimens$ROCNormalized <- specimens$ROC/specimens$BrainVol # New
specimens$CerebellarVol <- log10(specimens$CerebellarVol) # Reassign
specimens$CerebralVol <- log10(specimens$CerebralVol) # Reassign
specimens$BrainVol <- log10(specimens$BrainVol) # Reassign
specimens$CrusVol_1_3 <- log10(specimens$CrusVol_1_3) # Reassign
specimens$ROC <- log10(specimens$ROC) # Reassign


## Do the same for the per-specimen data?

## Shorter both dataframes to only include rows that have data for ansiform area volume.
speciesAA <- species %>% drop_na(MedianCrus)
speciesAA <- speciesAA[!(speciesAA$WilsonReederName=="Macaca_fascicularis"),]
specimensAA <- specimens %>% drop_na(CrusVol_1_3)

## Make phenolists, observations grouped per species.
phenolist.species <- split(species,species$WilsonReederName)
phenolist.specimens <- split(specimens, specimens$WilsonReederName)
phenolist.speciesAA <- split(speciesAA, speciesAA$WilsonReederName)
phenolist.specimensAA <- split(specimensAA, specimensAA$WilsonReederName)

# Reorder to match tip label order
phenolist.species <- phenolist.species[tree$tip.label];
phenolist.specimens <- phenolist.specimens[tree$tip.label]
phenolist.speciesAA <- phenolist.speciesAA[tree13$tip.label]
phenolist.specimensAA <- phenolist.specimensAA[tree13$tip.label]

# Get phenotypes
C <- lapply(phenolist.species,"[[","MedianCerebellum")
V <- lapply(phenolist.species,"[[","MedianCerebrum")
AA <- lapply(phenolist.species,"[[","MedianCrus")
ROC <- lapply(phenolist.species,"[[","MedianROC")
BV <- lapply(phenolist.species,"[[","BrainVol" )
C.norm <- lapply(phenolist.species,"[[","MedianCerebellumNormalized")
V.norm <- lapply(phenolist.species,"[[","MedianCerebrumNormalized")
AA.norm <- lapply(phenolist.species,"[[","MedianCrusNormalized")
ROC.norm <- lapply(phenolist.species,"[[","MedianROCNormalized")

# per specimen
C.specimens <- lapply(phenolist.specimens,"[[","CerebellarVol")
V.specimens <- lapply(phenolist.specimens,"[[","CerebralVol")
AA.specimens <- lapply(phenolist.specimens,"[[","CrusVol_1_3")
ROC.specimens <- lapply(phenolist.specimens,"[[","ROC")
BV.specimens <- lapply(phenolist.specimens,"[[","BrainVol")
C.specimens.norm <- lapply(phenolist.specimens,"[[","CerebellumNormalized")
V.specimens.norm <- lapply(phenolist.specimens,"[[","CerebrumNormalized")
AA.specimens.norm <- lapply(phenolist.specimens,"[[","CrusNormalized")
ROC.specimens.norm <- lapply(phenolist.specimens,"[[","ROCNormalized")


# Get phenotypes that are complete (N=13)
C.AA <- lapply(phenolist.speciesAA,"[[","MedianCerebellum")
V.AA <- lapply(phenolist.speciesAA,"[[","MedianCerebrum")
AA.AA <- lapply(phenolist.speciesAA,"[[","MedianCrus")
ROC.AA <- lapply(phenolist.speciesAA,"[[","MedianROC")
BV.AA <- lapply(phenolist.speciesAA,"[[","BrainVol" )
C.normAA <- lapply(phenolist.speciesAA,"[[","MedianCerebellumNormalized")
V.normAA <- lapply(phenolist.speciesAA,"[[","MedianCerebrumNormalized")
AA.normAA <- lapply(phenolist.speciesAA,"[[","MedianCrusNormalized")
ROC.normAA <- lapply(phenolist.speciesAA,"[[","MedianROCNormalized")

C.specimens.AA <- lapply(phenolist.specimensAA,"[[","CerebellarVol")
V.specimens.AA <- lapply(phenolist.specimensAA,"[[","CerebralVol")
AA.specimens.AA <- lapply(phenolist.specimensAA,"[[","CrusVol_1_3")
ROC.specimens.AA <- lapply(phenolist.specimensAA,"[[","ROC")
BV.specimens.AA <- lapply(phenolist.specimensAA,"[[","BrainVol" )
C.specimens.normAA <- lapply(phenolist.specimensAA,"[[","CerebellumNormalized")
V.specimens.normAA <- lapply(phenolist.specimensAA,"[[","CerebrumNormalized")
AA.specimens.normAA <- lapply(phenolist.specimensAA,"[[","CrusNormalized")
ROC.specimens.normAA <- lapply(phenolist.specimensAA,"[[","ROCNormalized")

# Compute contrasts
pic.C <- pic.ortho(C,tree,intra=TRUE)
pic.V <- pic.ortho(V,tree,intra=TRUE)
#pic.AA <- pic.ortho(AA, tree, intra=TRUE) ## This one does not really make sense. We are only able to calculate PIC for nodes 51 and 52. PICs are more sensibly calculated from the 13-species tree, that has complete data.
#pic.ROC <- pic.ortho(ROC, tree, intra=TRUE) ## Ibidem.
pic.BV <- pic.ortho(BV, tree, intra =TRUE)
pic.C.norm <- pic.ortho(C.norm, tree, intra =TRUE)
pic.V.norm  <- pic.ortho(V.norm, tree, intra =TRUE)
pic.AA.norm <- pic.ortho(AA.norm, tree, intra =TRUE)
pic.ROC.norm <- pic.ortho(ROC.norm, tree, intra =TRUE)

# Contrast for complete data (N=13)
pic.C.AA <- pic.ortho(C.AA,tree13,intra=TRUE)
pic.V.AA <- pic.ortho(V.AA,tree13,intra=TRUE)
pic.AA.AA <- pic.ortho(AA.AA, tree13, intra=TRUE)
pic.ROC.AA <- pic.ortho(ROC.AA, tree13, intra=TRUE) 
pic.BV.AA <- pic.ortho(BV.AA, tree13, intra =TRUE)
pic.C.normAA <- pic.ortho(C.normAA, tree13, intra =TRUE)
pic.V.normAA  <- pic.ortho(V.normAA, tree13, intra =TRUE)
pic.AA.normAA <- pic.ortho(AA.normAA, tree13, intra =TRUE)
pic.ROC.normAA <- pic.ortho(ROC.normAA, tree13, intra =TRUE)



## To be able to label the outliers in individual trait plots, we create a function to identify them ##
find_outlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

###########################################################
### PART 1: Data distribution
### Part 1a: Intraspecific distribution of Brain Volume
### Part 1b: ... Does this make sense? Maybe not.
### 
###########################################################

### ---------------------------- ###
# Box plots, highlighting outliers
### ---------------------------- ###    

# Create a viridis 6-color scale vector
sixcol <- viridis(6, alpha = 0.6, begin = 0, end = 1, direction = 1, option = "D")
names(sixcol) <- c('Black spider monkey', 'Central chimpanzee', 'Common squirrel monkey', 'Crab-eating macaque', 'Human', 'Rhesus monkey')

# Define a custom theme
custom_theme <- theme_pubr() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5, size = 15, margin = margin(r = 40)),
    axis.text.x = element_text(vjust = 0.4, size = 15, angle = 45),
    legend.position = "top",
    legend.text = element_text(),
    legend.key.size = unit(3, 'cm'),
    plot.margin = margin(20, 40, 40, 40),  # Adjust plot margins
    plot.background = element_rect(fill = "white"),  # Set white plot background
    panel.border = element_blank(),  # Remove panel borders
    panel.grid.major = element_line(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

# Create the plot
pdf('6.Brain_volume/6.1.boxplots_intraspecific.pdf', width = 12, height = 8)

box_brainvolume <- specimens[specimens$English.Name %in% c('Human', 'Central chimpanzee', 'Crab-eating macaque', 'Rhesus monkey', 'Black spider monkey', 'Common squirrel monkey'), ] %>%
  group_by(English.Name) %>%
  mutate(outlier = ifelse(find_outlier(BrainVol), BrainVol , NA)) %>%
  ggplot(aes(y = BrainVol, x = English.Name, fill = English.Name)) +
  geom_raster(fill = "white", interpolate = TRUE) +  # Add raster background with white color
  geom_boxplot(outlier.colour = "red", outlier.shape = "triangle", outlier.size = 3) +
  scale_fill_manual(values = sixcol) +
  geom_point(color = "black", size = 1.4, alpha = 0.9) +
  custom_theme +  # Apply custom theme
  ggtitle("Whole brain volume distribution") +
  labs(y = "Brain volumes (in mm3, log-transformed)", x = "") +
  geom_text(aes(label = outlier), na.rm = TRUE, vjust = 2)

print(box_brainvolume)
dev.off()



###########################################################
### PART 2: MAIN PGLS analyses, with Brain volume as covariate.
###########################################################

## PGLS under Brownian Motion should be the same as PIC-regression forced through 0.
# We use species median data here, so that neither the PGLS and PIC-model incorporate intraspecific variance.

# Model 1: Cerebellar volume regressed on cerebral volume. Accounting for brain volume.
bm<-corBrownian(1,tree)
pglsModel1 <-gls(MedianCerebellum~MedianCerebrum+BrainVol,data=species,correlation=bm)
summary(pglsModel1)
PICmodel1 <- lm(pic.C ~pic.V + pic.BV + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("6.Brain_volume/6.2.PGLS-PIC-C-V-accounting_BV.txt")
print(summary(pglsModel1))
print(coef(pglsModel1))
print("---------------------------")
print(summary(PICmodel1))
print(coef(PICmodel1))
sink()

# Model 1a: Cerebellar volume regressed on cerebral volume. Original.
bm<-corBrownian(1,tree)
pglsModel1a <-gls(MedianCerebellum~MedianCerebrum,data=species,correlation=bm)
summary(pglsModel1a)
PICmodel1a <- lm(pic.C ~pic.V + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("6.Brain_volume/6.2a.PGLS-PIC-C-V.txt")
print(summary(pglsModel1a))
print(coef(pglsModel1a))
print("---------------------------")
print(summary(PICmodel1a))
print(coef(PICmodel1a))
sink()

# Model 2: Ansiform area regressed on rest of cerebellar volume. Accounting for brain volume.
bm13<-corBrownian(1,tree13)
pglsModel2 <-gls(MedianCrus~MedianROC+BrainVol,data=speciesAA,correlation=bm13)
summary(pglsModel2)
PICmodel2 <- lm(pic.AA.AA ~pic.ROC.AA + pic.BV.AA + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("6.Brain_volume/6.3.PGLS-PIC-AA-ROC-accounting_BV.txt")
print(summary(pglsModel2))
print(coef(pglsModel2))
print("---------------------------")
print(summary(PICmodel2))
print(coef(PICmodel2))
sink()

# Model 2a: Ansiform area regressed on rest of cerebellar volume. Original.
bm13<-corBrownian(1,tree13)
pglsModel2a <-gls(MedianCrus~MedianROC,data=speciesAA,correlation=bm13)
summary(pglsModel2a)
PICmodel2a <- lm(pic.AA.AA ~pic.ROC.AA + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("6.Brain_volume/6.3a.PGLS-PIC-AA-ROC.txt")
print(summary(pglsModel2a))
print(coef(pglsModel2a))
print("---------------------------")
print(summary(PICmodel2a))
print(coef(PICmodel2a))
sink()

# Model 3: Ansiform area regressed on cerebral volume. Accounting for brain volume.
bm13<-corBrownian(1,tree13)
pglsModel3 <- gls(MedianCrus~MedianCerebrum+BrainVol,data=speciesAA,correlation=bm13)
summary(pglsModel3)
PICmodel3 <- lm(pic.AA.AA ~pic.V.AA + pic.BV.AA + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("6.Brain_volume/6.4.PGLS-PIC-AA-V-accounting_BV.txt")
print(summary(pglsModel3))
print(coef(pglsModel3))
print("---------------------------")
print(summary(PICmodel3))
print(coef(PICmodel3))
sink()

# Model 3a: Ansiform area regressed on cerebral volume. Original.
bm13<-corBrownian(1,tree13)
pglsModel3a <- gls(MedianCrus~MedianCerebrum,data=speciesAA,correlation=bm13)
summary(pglsModel3a)
PICmodel3a <- lm(pic.AA.AA ~pic.V.AA + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("6.Brain_volume/6.4a.PGLS-PIC-AA-V.txt")
print(summary(pglsModel3a))
print(coef(pglsModel3a))
print("---------------------------")
print(summary(PICmodel3a))
print(coef(PICmodel3a))
sink()

### -------------------------- ###
### Calculate R2s for PGLS models
### Now especially important because
### the R2lik function provides one
### of the few appropriate ways to
### test improvement of fit by intro-
### ducing a covariate in PGLS.
### -------------------------- ###   

sink("./6.Brain_volume/6.5.r2_phylogenetic_regressions.txt")
print("Cerebellar volume regressed on cerebral volume. With Brain Volume.")
mod1.r2 <- R2(mod = pglsModel1, phy = tree, sigma2_d = "s2w") #check out the meanings of the different R2s
mod1.r2
print("Cerebellar volume regressed on cerebral volume. Without Brain Volume.")
mod1a.r2 <- R2(mod = pglsModel1a, phy = tree, sigma2_d = "s2w") #check out the meanings of the different R2s
mod1a.r2 

print("Compare the models - difference cerebellar-cerebral regressions")
mod1.r2.compare <- R2(mod = pglsModel1, pglsModel1a, sigma2_d = "s2w") #check out the meanings of the different R2s
mod1.r2.compare

print("Ansiform area volume regressed on rest of cerebellar volume. With Brain Volume.")
mod2.r2 <- R2(mod = pglsModel2, phy = tree13, sigma2_d = "s2w")
mod2.r2 
print("Ansiform area volume regressed on rest of cerebellar volume. Without Brain Volume.")
mod2a.r2 <- R2(mod =  pglsModel2a, phy = tree, sigma2_d = "s2w") #check out the meanings of the different R2s
mod2a.r2

print("Compare the models - difference ansiform-ROC regressions")
mod2.r2.compare <- R2(mod = pglsModel2, pglsModel2a, sigma2_d = "s2w") #check out the meanings of the different R2s
mod2.r2.compare


print("Ansiform area volume regressed on cerebral volume. With Brain Volume.")
mod3.r2 <- R2(mod = pglsModel3, phy = tree13, sigma2_d = "s2w")
mod3.r2
print("Ansiform area volume regressed on cerebral volume. Without Brain Volume.")
mod3a.r2 <- R2(mod = pglsModel3a, phy = tree13, sigma2_d = "s2w")
mod3a.r2

print("Compare the models - difference - ansiform - cerebral regressions")
mod2.r2.compare <- R2(mod = pglsModel2, pglsModel2a, sigma2_d = "s2w") #check out the meanings of the different R2s
mod2.r2.compare
sink()

### VISUALIZE THE PGLS ANALYSES ###
## How to visualize partial PGLS regression? Probably no appropriate way to do it.

###########################################################
### PART 3: PGLS Regressions of phenotypes of interest
### versus Brain volume.
### This also adds to supplemental analyses 1, where we checked
### if cerebellar - to - cerebral scaling, and cerebellar and cerebral
### scaling versus body mass (Isler et al. 2008) were significantly
### different for strepsirrhines/haplorhines. However, the most
### appropriate comparison is seeing if cerebellar and cerebral
### volumes scale differently between both groups when compared
### to whole brain volumes recorded in the same specimens.
###########################################################

## Data preparation
#One could manually assign groupings and filter based on that.
rownames(species) <- species$WilsonReederName
STREP <- species %>% filter(Strepsirrhine== "Yes")
tree_strep <- treedata(tree,STREP, sort=T,warnings=T)$phy # Match tree to the strepsirrhine data

# The haplorhine data contains a species that had no body mass data, that we will remove from the tree as well.
HAPL <- species %>% filter(Strepsirrhine== "No")
tree_hapl<- treedata(tree,HAPL, sort=T,warnings=T)$phy # Match tree to the data
# This is an alternative way to subset the trees :). Very convenient
Strepsirrhines <-getTips(tree,findMRCA(tree,c("Galago_demidoff","Microcebus_murinus")))

#####################################
### FIRST: strepsirrhines
#####################################

Cstrep <- STREP$MedianCerebellum
Vstrep <- STREP$MedianCerebrum
BVstrep <- STREP$BrainVol 

# Compute contrasts
# Full
pic.C.strep <- pic(Cstrep,tree_strep)
pic.V.strep <- pic(Vstrep, tree_strep)
pic.BV.strep <- pic(BVstrep, tree_strep)

#------------------------------------
# PIC-regressions
#------------------------------------
## Cerebellum ##
# Compute a linear regression of C on BV
result <- lm(pic.C.strep ~ pic.BV.strep)
sink("6.Brain_volume/6.6.lm-C-BV.strepsirrhine.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on BV
pdf("6.Brain_volume/6.7.lm-C-BV-mean-nointercept-plot.strepsirrhine.pdf")
plot(pic.C.strep ~ pic.BV.strep)
abline(a = 0, b = coef(result))
graphics.off()

## Cerebellar PGLS
bm.strep <-corBrownian(1,tree_strep)
pglsModel<-gls(MedianCerebellum~BrainVol,data=STREP,correlation=bm.strep)
summary(pglsModel)
PICmodel <- lm(pic.C.strep ~pic.BV.strep + 0)
sink("6.Brain_volume/6.8.PGLS-PIC-C-BV.strepsirrhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()


## Cerebrum ##
# Compute a linear regression of V on BV
result <- lm(pic.V.strep ~ pic.BV.strep)
sink("6.Brain_volume/6.9.lm-V-BV.strepsirrhine.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of V on BV
pdf("6.Brain_volume/6.10.lm-V-BV-mean-nointercept-plot.strepsirrhine.pdf")
plot(pic.V.strep ~ pic.BV.strep)
abline(a = 0, b = coef(result))
graphics.off()

## Cerebellar PGLS
bm.strep <-corBrownian(1,tree_strep)
pglsModel<-gls(MedianCerebrum~BrainVol,data=STREP,correlation=bm.strep)
summary(pglsModel)
PICmodel <- lm(pic.V.strep ~pic.BV.strep + 0)
sink("6.Brain_volume/6.11.PGLS-PIC-V-BV.strepsirrhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

#####################################
### NEXT: Haplorhines
#####################################

Chapl <- HAPL$MedianCerebellum
Vhapl <- HAPL$MedianCerebrum
BVhapl <- HAPL$BrainVol 

# Compute contrasts
# Full
pic.C.hapl <- pic(Chapl,tree_hapl)
pic.V.hapl <- pic(Vhapl, tree_hapl)
pic.BV.hapl <- pic(BVhapl, tree_hapl)

#------------------------------------
# PIC-regressions
#------------------------------------
## Cerebellum ##
# Compute a linear regression of C on BV
result <- lm(pic.C.hapl ~ pic.BV.hapl)
sink("6.Brain_volume/6.12.lm-C-BV.haplorhine.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of C on BV
pdf("6.Brain_volume/6.13.lm-C-BV-mean-nointercept-plot.haplorhine.pdf")
plot(pic.C.hapl ~ pic.BV.hapl)
abline(a = 0, b = coef(result))
graphics.off()

## Cerebellar PGLS
bm.hapl <-corBrownian(1,tree_hapl)
pglsModel<-gls(MedianCerebellum~BrainVol,data=HAPL,correlation=bm.hapl)
summary(pglsModel)
PICmodel <- lm(pic.C.hapl ~pic.BV.hapl + 0)
sink("6.Brain_volume/6.14.PGLS-PIC-C-BV.haplorhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

## Cerebrum ##
# Compute a linear regression of V on BV
result <- lm(pic.V.hapl ~ pic.BV.hapl)
sink("6.Brain_volume/6.15.lm-V-BV.haplorhine.txt")
summary(result)
sink()

# Plot PIC regression
# Regression of V on BV
pdf("6.Brain_volume/6.16.lm-V-BV-mean-nointercept-plot.haplorhine.pdf")
plot(pic.V.hapl ~ pic.BV.hapl)
abline(a = 0, b = coef(result))
graphics.off()

## Cerebellar PGLS
bm.hapl<-corBrownian(1,tree_hapl)
pglsModel<-gls(MedianCerebrum~BrainVol,data=HAPL,correlation=bm.hapl)
summary(pglsModel)
PICmodel <- lm(pic.V.hapl ~pic.BV.hapl + 0)
sink("6.Brain_volume/6.17.PGLS-PIC-V-BV.haplorhine.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

###########################################################
### Use pANCOVA to test deviations in intercept or slope
### and with that deviation from the general allometric 
### relationship for cerebellar and cerebellar
### to brain volume scaling,
###########################################################
library(evomap)
library(geiger)

# This is an alternative way to subset the trees :). Very convenient
Strepsirrhines <-getTips(tree,findMRCA(tree,c("Galago_demidoff","Microcebus_murinus")))
Haplorhines <-getTips(tree,findMRCA(tree,c("Saimiri_sciureus","Cercopithecus_cephus_cephus")))

#For differences in slope: 
grpS<-rep("A",length(rownames(species))) 
grpS[Strepsirrhines]<-"B" 
grpS<-as.factor(grpS) 
names(grpS)<-rownames(species)

#For differences in intercept: 
grpI<-rep("A",length(rownames(species))) 
grpI[Strepsirrhines]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(species)

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
data <- data[1:34,]
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
tree <-treedata(tree,species, sort=T,warnings=T)$ph
rownames(data) <- tree[["tip.label"]]

## Again, extremely hacky, but works.
data$Dependent1 <- species$MedianCerebellum
data$Dependent2 <- species$MedianCerebrum
data$Independent <- species$BrainVol

## Plot separate PGLS regressions for Strepsirrhines and Haplorhines
## Cerebellum vs. BrainVolume
pdf("6.Brain_volume/6.17.PGLS.Strepsirrhine.cerebellum-brainvol.pdf")
plot(data$Dependent1~data$Independent,col="black",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent1","Independent",data,tree,model="BM",group=Strepsirrhines,col="pink",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent1","Independent",data,tree,model="BM",group=Haplorhines,col="azure4",lwd=5,cex=1,pch=19)
graphics.off()

## Model construction
Model.CBV <-model.matrix(as.formula(Dependent1~Independent),data)
Model.CBV.s <-model.matrix(as.formula(Dependent1~grpS:Independent),data) #slopes differ
Model.CBV.i <-model.matrix(as.formula(Dependent1~grpI + Independent),data) #intercepts differ
Model.CBV.si <-model.matrix(as.formula(Dependent1~grpI + grpS:Independent),data) #slopes & intercepts differ

sink("6.Brain_volume/6.18.ModelTesting.SlopesIntercepts.cerebellumbrainvol.txt")
## Testing versus baseline model (common intercept and slope)
print("(1) Differences in slopes, holding intercept constant:")
gls.ancova(Dependent1~Independent,vcv(tree),Model.CBV,Model.CBV.s)

print("(2) Differences in intercept, holding slopes constant:")
gls.ancova(Dependent1~Independent,vcv(tree),Model.CBV,Model.CBV.i)

print("(3) Differences in slopes and differences in intercept:")
gls.ancova(Dependent1~Independent,vcv(tree),Model.CBV,Model.CBV.si)
sink()


## Cerebrum vs. BrainVolume
pdf("6.Brain_volume/6.19.PGLS.Strepsirrhine.cerebrum-brainvol.pdf")
plot(data$Dependent2~data$Independent,col="black",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent2","Independent",data,tree,model="BM",group=Strepsirrhines,col="pink",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent2","Independent",data,tree,model="BM",group=Haplorhines,col="azure4",lwd=5,cex=1,pch=19)
graphics.off()

## Model construction
Model.VBV <-model.matrix(as.formula(Dependent2~Independent),data)
Model.VBV.s <-model.matrix(as.formula(Dependent2~grpS:Independent),data) #slopes differ
Model.VBV.i <-model.matrix(as.formula(Dependent2~grpI + Independent),data) #intercepts differ
Model.VBV.si <-model.matrix(as.formula(Dependent2~grpI + grpS:Independent),data) #slopes & intercepts differ

sink("6.Brain_volume/6.20.ModelTesting.SlopesIntercepts.cerebrum-brainvol.txt")
## Testing versus baseline model (common intercept and slope)
print("(1) Differences in slopes, holding intercept constant:")
gls.ancova(Dependent2~Independent,vcv(tree),Model.VBV,Model.VBV.s)

print("(2) Differences in intercept, holding slopes constant:")
gls.ancova(Dependent2~Independent,vcv(tree),Model.VBV,Model.VBV.i)

print("(3) Differences in slopes and differences in intercept:")
gls.ancova(Dependent2~Independent,vcv(tree),Model.VBV,Model.VBV.si)
sink()

## As a reminder: cerebellar vs. cerebrum scaling between both groups. 
## Here we also added the one species that we had missing body mass data for (excluded in supplemental 5.strepsirrhines)
## We also add all three allometries together in one plot.

## Cerebrum vs. BrainVolume
pdf("6.Brain_volume/6.21.PGLS.Strepsirrhine.cerebellum-cerebrum.pdf")
plot(data$Dependent1~data$Dependent2,col="black",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
pGLS.plotGrade("Dependent1","Dependent2",data,tree,model="BM",group=Strepsirrhines,col="pink",lwd=5,cex=1,pch=19)
pGLS.plotGrade("Dependent1","Dependent2",data,tree,model="BM",group=Haplorhines,col="azure4",lwd=5,cex=1,pch=19)
graphics.off()

## Model construction
Model.CC<-model.matrix(as.formula(Dependent1~Dependent2),data)
Model.CC.s <-model.matrix(as.formula(Dependent1~grpS:Dependent2),data) #slopes differ
Model.CC.i <-model.matrix(as.formula(Dependent1~grpI + Dependent2),data) #intercepts differ
Model.CC.si <-model.matrix(as.formula(Dependent1~grpI + grpS:Dependent2),data) #slopes & intercepts differ

sink("6.Brain_volume/6.22.ModelTesting.SlopesIntercepts.cerebellum-cerebrum.txt")
## Testing versus baseline model (common intercept and slope)
print("(1) Differences in slopes, holding intercept constant:")
gls.ancova(Dependent1~Dependent2,vcv(tree),Model.CC,Model.CC.s)

print("(2) Differences in intercept, holding slopes constant:")
gls.ancova(Dependent1~Dependent2,vcv(tree),Model.CC,Model.CC.i)

print("(3) Differences in slopes and differences in intercept:")
gls.ancova(Dependent1~Dependent2,vcv(tree),Model.CC,Model.CC.si)
sink()


## Plot them together
pdf("6.Brain_volume/6.23.PGLS.Strepsirrhine.composite.pdf", width = 16, height = 8)
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


###########################################################
### PART 4: Recalculate ancestral character estimations for 
### the cerebellum, normalized against the whole brain volume.
### We can then compare this with the cerebellar/cerebral 
### ancestral character estimation from the original analyses (part 4).
###########################################################

#------------------------------------
# FIT BROWNIAN MOTION MODELS
#------------------------------------
source("./input/catn.R")

## Construct the same BM model as in the main analysis; only added Brain Volume.
#pheno <- read.csv2("../Final_code/0.cleanInput/Phenotypes_outliersRM.csv", sep=';')  %>% dplyr::rename(species = WilsonReederName) ## This is the relative path to the original code folder. Change accordingly.
pheno <- specimens %>% dplyr::rename(species=WilsonReederName)
pheno <- pheno[c("species","CerebellarVol", "SurfaceArea" , "CerebralVol","AbsGI","FoldingLength","FoldingNumber","Lambda","Delta", "CrusVol_1_3", "CerebroCerebellar", "CerebellarCrus", "BrainVol", "CerebellumNormalized", "CerebrumNormalized", "CrusNormalized")]

phenoAA <- specimensAA %>% dplyr::rename(species=WilsonReederName)
phenoAA <- phenoAA[c("species","CerebellarVol", "SurfaceArea" , "CerebralVol","AbsGI","FoldingLength","FoldingNumber","Lambda","Delta", "CrusVol_1_3", "CerebroCerebellar", "CerebellarCrus", "BrainVol", "CerebellumNormalized", "CerebrumNormalized", "CrusNormalized")]

library(Rphylopars)
# The full table of measurements from the current study and Heuer et al., 2019.
catn("fit a Brownian Motion model")
p.BM <- Rphylopars::phylopars(trait_data=pheno, tree=tree)
sink("6.Brain_volume/6.24.BM_full.txt")
p.BM
sink()

## For the 13 complete-data species only
catn("fit a Brownian Motion model")
p.BM13 <- Rphylopars::phylopars(trait_data=phenoAA, tree=tree13)
sink("6.Brain_volume/6.24a.BM_13species.txt")
p.BM13
sink()

#------------------------------------------
# PLOT ANCESTRAL STATE ESTIMATIONS
#------------------------------------------

## Define a two-color scaled function for plotting ancestral state estimations along the phylogenetic tree. ##
myplotanc <- function(tree, tips, tipnames, states, mytitle) {
  names(tips) <- tipnames
  obj <- contMap(tree,tips,method="user",anc.states=states, plot="F")
  obj <- setMap(obj, colors=c("white", "darkblue"))
  plot.contMap(obj)
  axisPhylo()
  title(mytitle)
}

myplotanc.ci.Taxo <- function(tree, tips, tipnames, states, legend.title, rangeCIs) {
  names(tips) <- tipnames
  obj <- contMap(tree,tips,method="user",anc.states=states, plot="F", , lims = rangeCIs)
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

myplotanc.ci.noTaxo <- function(tree, tips, tipnames, states, legend.title, rangeCIs) {
  names(tips) <- tipnames
  obj <- contMap(tree,tips,method="user",anc.states=states, plot="F", ftype = "off", lims = rangeCIs)
  obj <- setMap(obj, colors=c("blue", "white", "red"))
  plot.contMap(obj, legend=FALSE,ylim=c(1-0.09*(Ntip(obj$tree)-1),Ntip(obj$tree)),
               mar=c(5.1,0.4,0.4,0.4), ftype = "off")
  add.color.bar(36.502,obj$cols,title= legend.title, y=3, fsize=1.0,
                lims=obj$lims,digits=2,prompt=FALSE,x=0,
                y=1-0.08*(Ntip(obj$tree)-1), fsize=0.8,subtitle="")
  axisPhylo()
  title(xlab="Time to present (in million years)", adj = 0.17)
  errorbar.contMap(obj, lwd=2)
}

myplotanc.ci.TaxoAA <- function(tree, tips, tipnames, states, legend.title) {
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
}


catn("plot ancestral values")

# ntips Number of species we have phenotypes for
ntips <- length(tree$tip.label)
# Number of nodes inside the tree (ntips-1) + including the tips
nnodes <- ntips + tree$Nnode

# ntips Number of species we have phenotypes for
ntips13 <- length(tree13$tip.label)
# Number of nodes inside the tree (ntips-1) + including the tips
nnodes13 <- ntips13 + tree13$Nnode

tipnames <- rownames(p.BM$anc_recon)[1:ntips]
tipnames13 <- rownames(p.BM13$anc_recon)[1:ntips13]

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
write.csv(csv, file="./6.Brain_volume/6.25.ACE_fullmodel-inclBrain.csv")

# Ancestral characters brain volume (10logged mm3)
pdf("./6.Brain_volume/6.26.BM-anc-BrainVolume.pdf", width = 10, height = 10)
myplotanc.ci.Taxo(tree, p.BM$anc_recon[1:ntips,"BrainVol"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"BrainVol"], "Reconstruction brain volume")
graphics.off()

# Ancestral characters cerebellum normalized against brain volume (ratio)
pdf("./6.Brain_volume/6.27.BM-anc-cerebellum-brainNormalized.pdf", width = 8, height = 10)
myplotanc.ci.noTaxo(tree, p.BM$anc_recon[1:ntips,"CerebellumNormalized"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"CerebellumNormalized"], "Reconstruction cerebellar volume \n (normalised against brain volume)", rangeCIs = c(0.072,0.165))
graphics.off()


## Additionally; let's see how cerebral and crura ACEs look like.
# Ancestral characters cerebral volume (10logged mm3)
pdf("./6.Brain_volume/6.28.BM-anc-cerebrum-brainNormalized.pdf", width = 8, height = 10)
myplotanc.ci.noTaxo(tree, p.BM$anc_recon[1:ntips,"CerebrumNormalized"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"CerebrumNormalized"], "Reconstruction cerebral volume \n (normalised against brain volume)", rangeCIs = c(0.636,0.858))
graphics.off()

# Ancestral characters cerebellum normalized against brain volume (ratio) ## JUST CURIOSITY. MAKES LITTLE SENSE.
pdf("./6.Brain_volume/6.29.BM-anc-ansiform-brainNormalized.pdf", width = 8, height = 10)
myplotanc.ci.TaxoAA(tree, p.BM$anc_recon[1:ntips,"CrusNormalized"], tipnames, p.BM$anc_recon[(ntips+1):nnodes,"CrusNormalized"], "Reconstruction ansiform area volume \n (normalised against brain volume)")
graphics.off()


###########################################################
### PART 5: Replicate encephalisation quotient in ordinary
### and phylogenetic regressions.
###########################################################

specimens2 <- specimens %>% drop_na(BodyMass) 
species2 <- species2 %>% drop_na(Body.mass.species.mean)

## Ordinary regression
sink("6.Brain_volume/6.30.regular_regressions_encephalisation.txt")
print("Brain volume regressed on body mass.")
reg.mod1 <- lm(BrainVol ~ BodyMass, data = specimens2)
summary(reg.mod1)
sink()

## Ordinary regression; species median
sink("6.Brain_volume/6.30a.regular_regressions_encephalisation-median.txt")
print("Brain volume regressed on body mass.")
reg.mod2 <- lm(BrainVol ~ Body.mass.species.mean, data = species2)
summary(reg.mod2)
sink()


### ---------------------------------------------------------------- ###
### GGPLOT Regular linear regressions plots, with allometric formulae
### ---------------------------------------------------------------- ### 

# Credit: Johnston, S: sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/.
# Function was adapted to specific purpose of the paper.

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

ggplotRegression2 <- function(fit, ylabel, xlabel, specimens) {
  require(ggplot2)
  
  # Step 1: Fit the linear regression model
  x <- fit$model[, 2]
  y <- fit$model[, 1]
  
  # Retrieve species information
  species_data <- specimens2$WilsonReederName
  
  # Step 2: Calculate the residuals
  residuals <- fit$residuals
  
  # Step 3: Calculate mean and standard deviation of residuals
  mean_residual <- mean(residuals)
  std_residual <- sd(residuals)
  
  # Step 4: Calculate number of standard deviations away from the mean for each point
  std_deviations <- residuals / std_residual
  
  # Filter out missing or inconsistent data points
  valid_indices <- !is.na(x) & !is.na(std_deviations) & !is.na(species_data)
  x <- x[valid_indices]
  std_deviations <- std_deviations[valid_indices]
  species_data <- species_data[valid_indices]
  
  # Create a new data frame with the required variables
  data <- data.frame(
    x = x,
    y = std_deviations,
    species = species_data
  )
  
  # Define the color palette for species
  color_palette <- rainbow(length(unique(species_data)))
  
  # Step 5: Modify the ggplot code to use the merged data frame
  ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = species, size = 2)) +
    scale_color_manual(values = color_palette) +
    geom_smooth(method = "lm", col = "black") +
    labs(
      title = paste(
        "y = ",
        signif(fit$coef[[1]], 5),
        " + ",
        signif(fit$coef[[2]], 5),
        "x",
        "         R2 Adj. = ",
        signif(summary(fit)$adj.r.squared, 5),
        " p-value =",
        signif(summary(fit)$coef[2, 4], 5)
      ),
      y = "Number of Standard Deviations",
      x = xlabel,
      color = "Species"  # Set the color legend label
    )
}


ggplotRegression3 <- function(fit, ylabel, xlabel, species) {
  require(ggplot2)
  
  # Step 1: Fit the linear regression model
  x <- fit$model[, 2]
  y <- fit$model[, 1]
  
  # Retrieve species information
  species_data <- species2$Eqplot
  
  # Step 2: Calculate the residuals
  residuals <- fit$residuals
  
  # Step 3: Calculate mean and standard deviation of residuals
  mean_residual <- mean(residuals)
  std_residual <- sd(residuals)
  
  # Step 4: Calculate number of standard deviations away from the mean for each point
  std_deviations <- residuals / std_residual
  
  # Filter out missing or inconsistent data points
  valid_indices <- !is.na(x) & !is.na(std_deviations) & !is.na(species_data)
  x <- x[valid_indices]
  std_deviations <- std_deviations[valid_indices]
  species_data <- species_data[valid_indices]
  
  # Create a new data frame with the required variables
  data <- data.frame(
    x = x,
    y = std_deviations,
    species = species_data
  )
  
  # Define the color palette for species
  color_palette <- rainbow(length(unique(species_data)))
  
  # Step 5: Modify the ggplot code to use the merged data frame
  ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = species), size = 7) +
    scale_color_manual(values = color_palette) +
    geom_text(aes(label = species_data), hjust = 1, vjust = -1, size=2.0, position = position_jitter(width = 0.07, height = 0.07)) +
    geom_smooth(method = "lm", col = "black") +
    labs(
      title = paste(
        "y = ",
        signif(fit$coef[[1]], 5),
        " + ",
        signif(fit$coef[[2]], 5),
        "x",
        "         R2 Adj. = ",
        signif(summary(fit)$adj.r.squared, 5),
        " p-value =",
        signif(summary(fit)$coef[2, 4], 5)
      ),
      y = "Number of Standard Deviations",
      x = xlabel,
      color = "Species"  # Set the color legend label
    )
}

ggplotRegression4 <- function(fit, ylabel, xlabel, species) {
  require(ggplot2)
  
  # Step 1: Fit the linear regression model
  x <- fit$model[, 2]
  y <- fit$model[, 1]
  
  # Retrieve species information
  species_data <- species2$EnglishName
  
  # Step 2: Calculate the residuals
  residuals <- fit$residuals
  
  # Step 3: Calculate mean and standard deviation of residuals
  mean_residual <- mean(residuals)
  std_residual <- sd(residuals)
  
  # Step 4: Calculate number of standard deviations away from the mean for each point
  std_deviations <- residuals / std_residual
  
  # Filter out missing or inconsistent data points
  valid_indices <- !is.na(x) & !is.na(std_deviations) & !is.na(species_data)
  x <- x[valid_indices]
  std_deviations <- std_deviations[valid_indices]
  species_data <- species_data[valid_indices]
  
  # Create a new data frame with the required variables
  data <- data.frame(
    x = x,
    y = std_deviations,
    species = species_data
  )
  
  # Define the color palette for species
  color_palette <- rainbow(length(unique(species_data)))
  
  # Step 5: Modify the ggplot code to use the merged data frame
  ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = species_data), size = 7) +
    scale_color_manual(values = color_palette) +
    geom_smooth(method = "lm", col = "black") +
    labs(
      title = paste(
        "y = ",
        signif(fit$coef[[1]], 5),
        " + ",
        signif(fit$coef[[2]], 5),
        "x",
        "\nR2 Adj. = ",
        signif(summary(fit)$adj.r.squared, 5),
        " p-value =",
        signif(summary(fit)$coef[2, 4], 5)
      ),
      y = ylabel,
      x = xlabel,
      color = "Species"  # Set the color legend label
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14, margin = margin(t = 5, r = 5, b = 5, l = 5)),
      axis.title = element_text(size = 16, margin = margin(t = 20, r = 20, b = 20, l = 20)),
      legend.title = element_text(size = 14),  # Adjust the legend title size
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, hjust = 0.5, vjust = -0.5, margin = margin(b = 20)),
      plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines")
    )
  # Return the dataframe with number of standard deviations
  return(data)
}


ggplotRegression5 <- function(fit, ylabel, xlabel, species) {
  require(ggplot2)
  
  # Step 1: Fit the linear regression model
  x <- fit$model[, 2]
  y <- fit$model[, 1]
  
  # Retrieve species information
  species_data <- species2$Eqplot
  
  # Step 2: Calculate the residuals
  residuals <- fit$residuals
  
  # Step 3: Calculate mean and standard deviation of residuals
  mean_residual <- mean(residuals)
  std_residual <- sd(residuals)
  
  # Step 4: Calculate number of standard deviations away from the mean for each point
  std_deviations <- residuals / std_residual
  
  # Filter out missing or inconsistent data points
  valid_indices <- !is.na(x) & !is.na(std_deviations) & !is.na(species_data)
  x <- x[valid_indices]
  std_deviations <- std_deviations[valid_indices]
  species_data <- species_data[valid_indices]
  
  # Create a new data frame with the required variables
  data <- data.frame(
    x = x,
    y = std_deviations,
    species = species_data
  )
  
  # Define the color palette for species
  #color_palette <- rainbow(length(unique(species_data)))
  color_palette <- c('lightskyblue', 'azure4','plum3')
  
  # Step 5: Modify the ggplot code to use the merged data frame
  ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = species_data), size = 7) +
    scale_color_manual(values = color_palette) +
    geom_smooth(method = "lm", col = "black") +
    labs(
      title = paste(
        "y = ",
        signif(fit$coef[[1]], 5),
        " + ",
        signif(fit$coef[[2]], 5),
        "x",
        "\nR2 Adj. = ",
        signif(summary(fit)$adj.r.squared, 5),
        " p-value =",
        signif(summary(fit)$coef[2, 4], 5)
      ),
      y = ylabel,
      x = xlabel,
      color = "Species"  # Set the color legend label
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14, margin = margin(t = 5, r = 5, b = 5, l = 5)),
      axis.title = element_text(size = 16, margin = margin(t = 20, r = 20, b = 20, l = 20)),
      legend.title = element_text(size = 14),  # Adjust the legend title size
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, hjust = 0.5, vjust = -0.5, margin = margin(b = 20)),
      plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines")
    )
}

## Mostly to resemble the original "Encephalisation Quotient" figures.
tiff("./6.Brain_volume/6.31.regularRegression-brainvol-bodymass.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplotRegression(reg.mod1, ylabel = "Brain volume (log10-transformed)", xlabel = "Body mass (species mean; log10-transformed)") + theme_classic()  +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
graphics.off()

tiff("./6.Brain_volume/6.31a.regularRegression-brainvol-bodymass-sd.tiff", width = 24, height = 8, units = "in", res = 300)
ggplotRegression2(fit = reg.mod1, ylabel = "Brain volume (log10-transformed)", xlabel = "Body mass (species mean; log10-transformed)", specimens = specimens2) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),  # Adjust the legend title size
        legend.text = element_text(size = 12))
graphics.off()


# Only species median values
tiff("./6.Brain_volume/6.32.regularRegression-brainvol-bodymass-median.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplotRegression(reg.mod2, ylabel = "Brain volume (species median, log10-transformed)", xlabel = "Body mass (species mean; log10-transformed)") + theme_classic()  +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
graphics.off()

tiff("./6.Brain_volume/6.32a.regularRegression-brainvol-bodymass-median-sd.tiff", width = 20, height = 8, units = "in", res = 300)
ggplotRegression3(reg.mod2, ylabel = "Brain volume (species, median log10-transformed)", xlabel = "Body mass (species mean; log10-transformed)", species = species2) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),  # Adjust the legend title size
        legend.text = element_text(size = 12))
graphics.off()

sink("./6.Brain_volume/6.32b.EQ-Species.txt")
ggplotRegression4(reg.mod2, ylabel = "Number of standard deviations\nin brain volume (species median; log10-transformed)", xlabel = "Body mass (species mean; log10-transformed)", species = species2)
sink()

tiff("./6.Brain_volume/6.32c.regularRegression-brainvol-bodymass-median-sd-notext.tiff", width = 20, height = 8, units = "in", res = 300)
ggplotRegression5(reg.mod2, ylabel = "Number of standard deviations\nin brain volume (species median; log10-transformed)", xlabel = "Body mass (species mean; log10-transformed)", species = species2)
graphics.off()

### ---------------------------------------------------------------- ###
### PGLS REGRESSION and VISUALISATION
### ---------------------------------------------------------------- ### 

#####################################
### CALCULATE PICs
#####################################
tree33 <- drop.tip(tree, tip = "Trachypithecus_(Trachypithecus)_germaini")
BV <- species2$BrainVol
BM <- species2$Body.mass.species.mean

# Compute contrasts
pic.BV <- pic(BV,tree33)
pic.BM <- pic(BM, tree33)

## PGLS under Brownian Motion should be the same as PIC-regression forced through 0.
# We use species median data here, so that neither the PGLS and PIC-model incorporate intraspecific variance.
bm<-corBrownian(1,tree33)
pglsModel<-gls(BrainVol~Body.mass.species.mean,data=species2,correlation=bm)
summary(pglsModel)
PICmodel <- lm(pic.BV ~pic.BM + 0) ## NOT YET THE SAME AS PGLS, check what is up.
sink("6.Brain_volume/6.34.PGLS-PIC-C-V.txt")
print(summary(pglsModel))
print(coef(pglsModel))
print("---------------------------")
print(summary(PICmodel))
print(coef(PICmodel))
sink()

### -------------------------- ###
### PGLS Plotting with Caper
### -------------------------- ###   
library(caper)
library(evomap)

# Model 1: Brain Volume Regressed on Body Mass
tiff("./6.Brain_volume/6.35.BM-PGLS-BrainVol-BodyMass-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod1 <- gls(BrainVol~Body.mass.species.mean,data=species2, correlation=bm)
plot(mod1)
plot(species2$BrainVol~species2$Body.mass.species.mean, ylab = "Brain volume (in mm3; log10-transformed)", xlab = "Body mass (in kg; log10-transformed)")
abline(mod1)
abline(a=4.139046 , b=1.0, col= "blue") #MANUAL INPUT
# Confidence intervals
pGLS.mod1.ci <-gls.ci(species2$BrainVol, species2$Body.mass.species.mean,vcv(tree33))
lines(pGLS.mod1.ci$CI.plot$X,pGLS.mod1.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod1.ci$CI.plot$X,pGLS.mod1.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod1.pi <-gls.pi(species2$BrainVol, species2$Body.mass.species.mean,vcv(tree33), 1)
lines(pGLS.mod1.pi$PI.plot$X,pGLS.mod1.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod1.pi$PI.plot$X,pGLS.mod1.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()


## Encephalisation quotient formula
## MAMMAL EQ formula: EQ = Brain weight (g) / 0.12*bodyweight (g)^(2/3) (Jerison 1979)
# Brain density: 1,100 gram/cm3, https://doi.org/10.1111/j.1600-0404.1970.tb05606.
EQ <- (((10^species$BrainVol)/1000)*1.1) / (0.12*((10^species$Body.mass.species.mean)*1000)^(2/3))
species3 <- species
species3$EQ <- EQ
names(EQ) <- species$WilsonReederName

###########################################################
### END OF SCRIPT
###########################################################
