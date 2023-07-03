
### ----------------------------------- ###
### Data Cleaning - Full phenotype tables
### Neville Magielse, 2022
### ----------------------------------- ###        

# We are using R version 4.1.0 (2021-05-18)
# Nota bene: we specified dplyr functions specifically in our code i.e., dplyr::select. Its functions otherwise clashed with other packages in our environment.

### SET UP ###
getwd()
setwd("~/Documents/PhD/EvolutionPrimates/Final_code/")

outputdir <- "0.cleanInput"
dir.create(outputdir)

library(tidyverse)
library(magrittr)

### INPUT ####

# Load in all the required data. All data can be found in the Github Repository accompanying this project.
# Load in the species lookup table  cerebellar (and cerebral) volumes (https://colab.research.google.com/drive/1DIZgL63Am4uyv5nu41cwxcTBfEfgsHcG#scrollTo=6zDtFHuiP7lW)
# And the Crus I/II volumes (https://colab.research.google.com/drive/1CkPf6SmP3MsSA7ovafPN6wewh3Nmtm6Q#scrollTo=84YQivOtOqF7)
# Stats2 is an adapted version of (https://github.com/neuroanatomy/34primates/blob/master/data/derived/stats/stats.csv) and includes all neocortical measurements from Heuer et al., 2019
lut <- read_csv2("./input/10kTrees_34PrimateSpecies_adapted.csv")  #look-up table
Cerebellum <- read_csv2("./input/dfCerebellum.csv") #Cerebellar volumes
Neville <- read_csv2("./input/dfCrusNeville.csv") # Combined Crus I/II segmentations made by Neville Magielse
Mahta <- read_csv2("./input/dfCrusMahta.csv") # Combined Crus I/II segmentations made by Mahta Abbaspour
Vanessa <- read_csv2("./input/dfCrusVanessa.csv") #Combined Crus I/II segmentations made by Vanessa Steigauf
Neocortical <- read_csv2("./input/stats2.csv") # Neocortical measurements from Heuer et al. 2019

## Body weight data
## Was ultimately not used. See discussion in the manuscript.
BodyWeight <- read_csv2("./input/bodyweight-isleretal2008.csv") %>% dplyr::select(c("Species_Full","Body mass species mean")) #Placeholder average body weight data collected by Isler et al. 2008 
names(BodyWeight)[1] <- 'WilsonReederName' #Change column name for merging purposes
BodyWeight$WilsonReederName <- sub(" ", "_", BodyWeight$WilsonReederName) # And sub all spaces in the species names to underscores, also for matching purposes.
BodyWeight$`Body mass species mean` <- (BodyWeight$`Body mass species mean`)/1000 # Converts to kilos
# Change occurrences of names that are slightly different in the Body Weight data from the Wilson & Reeder names
BodyWeight$WilsonReederName[BodyWeight$WilsonReederName == 'Cercopithecus_cephus'] <- 'Cercopithecus_cephus_cephus'
BodyWeight$WilsonReederName[BodyWeight$WilsonReederName == 'Gorilla_beringei'] <- 'Gorilla_beringei_graueri'
BodyWeight$WilsonReederName[BodyWeight$WilsonReederName == 'Gorilla_gorilla'] <- 'Gorilla_gorilla_gorilla'
BodyWeight$WilsonReederName[BodyWeight$WilsonReederName == 'Varecia_variegata'] <- 'Varecia_variegata_variegata'
# Trachypithecus_(Trachypithecus)_germaini has no Body Weight data available

## Write a new Body Weight table
write_csv2(BodyWeight, "./0.cleanInput/BodyWeight.csv")

### MERGING and DATA CLEAN-UP ###

## Merge all the data together, after which we will do some cleaning
phenotypes <- left_join(lut, Cerebellum,by="SpecimenID") %>% left_join(Neocortical, by="SpecimenID") %>% left_join(Neville, by="SpecimenID") %>% left_join(Mahta, by="SpecimenID") %>% left_join(Vanessa, by="SpecimenID")
  
# We change the name of some columns, and delete the old and redundant cerebellar volume column. 
# We also add a column where we will put Crus observations from obs. 1&3 together.
# Additionally, we add body weight data (placeholder data, obtained from the collection of Isler et al. 2008 @https://doi.org/10.1016/j.jhevol.2008.08.004), that can now be merged because we changed the WilsonReederName column name.
phenotypes <- phenotypes %>% dplyr::rename(SurfaceArea =`Surface Area`) %>% dplyr::rename(FoldingLength =`Folding Length`) %>% dplyr::rename(FoldingNumber =`Folding number`) %>% dplyr::rename(CerebralVol=Volume) %>% dplyr::rename(CrusVol_obs1=CrusVol.x) %>% dplyr::rename(CrusVol_obs2=CrusVol.y) %>%  dplyr::rename(CrusVol_obs3=CrusVol) %>% add_column(CrusVol_1_3 = 0, .after = "CrusVol_obs3") %>% dplyr::select(-c(Cerebellum)) %>% left_join(BodyWeight, by="WilsonReederName") %>% dplyr::rename(BodyMass =`Body mass species mean`)  

# Then, let's combine the Crus I/II measurements from observer 1 and observer 3 in one column, as well. Values from observer 2 are only used to validate the manual segmentation method, and are therefore omitted here.
# This first necessitates changing the 0.0s into NAs ###Or does it? Remove from final script if not necessary.###
phenotypes[phenotypes==0]<-NA
phenotypes <- phenotypes %>% mutate(CrusVol_1_3 = coalesce(CrusVol_obs1,CrusVol_obs3))

## Next, we add columns for the ratios of Cerebellar/Cerebral volume, and Crus/Cerebellar Volume
# And some more switching between 0s and NaNs #phenotypes[is.na(phenotypes)] = 0 ###Potentially remove as well###
phenotypes <- phenotypes %>%
  mutate(CerebroCerebellar = (CerebellarVol/CerebralVol)*100) %>% mutate(CerebellarCrus = (CrusVol_1_3/CerebellarVol)*100)

## We also add columns for cerebellar, cerebral, and Crus Volumes relative to body mass.
phenotypes <- phenotypes %>%
  mutate(CerebellarMass = CerebellarVol/BodyMass) %>% mutate(CerebralMass = CerebralVol/BodyMass) %>% mutate(CrusMass = CrusVol_1_3/BodyMass) ## %>% mutate(CerebroCerebellarMass = CerebroCerebellar/BodyMass) %>% mutate(CerebellarCrusMass = CerebellarCrus/BodyMass)

# Lastly, we add a columns where we subtract Crura I-II volumes from whole cerebellar volume. This then allows up to regress Crura I-II versus 'rest of cerebellum'.
phenotypes <- phenotypes %>%
  mutate(ROC = CerebellarVol - CrusVol_1_3) %>%
  mutate(ROCCrus = ((CrusVol_1_3)/(ROC))*100)


## Check if the data is well-structured and tidy.
class(phenotypes)
str(phenotypes)

### RAW PHENOTYPES ###
phenotypes_raw <- phenotypes

## Write the raw phenotypes into a file
write_csv2(phenotypes_raw, "./0.cleanInput/Phenotypes_raw.csv")

### LOG-TRANSFORMATION ###

## Perform logarithmic transformation for all measurements that differ across orders of magnitude
phenotypes$CerebellarVol <- log10(phenotypes$CerebellarVol)
phenotypes$SurfaceArea <- log10(phenotypes$SurfaceArea)
phenotypes$CerebralVol <- log10(phenotypes$CerebralVol)
#phenotypes_tidy$AbsGI <- log10(phenotypes_tidy$AbsGI) # not necessary to transform (only 1 order of magnitude)
phenotypes$FoldingLength <- log10(phenotypes$FoldingLength)
phenotypes$FoldingNumber <- log10(phenotypes$FoldingNumber)
phenotypes$Delta[phenotypes$Delta<0] <- 0 #Change delta values smaller than 0 to 0
phenotypes$CrusVol_obs1 <- log10(phenotypes$CrusVol_obs1)
phenotypes$CrusVol_obs2 <- log10(phenotypes$CrusVol_obs2)
phenotypes$CrusVol_obs3 <- log10(phenotypes$CrusVol_obs3)
phenotypes$CrusVol_1_3 <- log10(phenotypes$CrusVol_1_3)
phenotypes$BodyMass <- log10(phenotypes$BodyMass)
#phenotypes$CerebroCerebellar <- log10(phenotypes$CerebroCerebellar) #not necessary to transform (percentages (1 order of magnitude))
#phenotypes$CerebellarCrus <- log10(phenotypes$CerebellarCrus) #not necessary to transform (percentages (1 order of magnitude))
phenotypes$CerebellarMass <- log10(phenotypes$CerebellarMass)
phenotypes$CerebralMass <- log10(phenotypes$CerebralMass)
phenotypes$CrusMass <- log10(phenotypes$CrusMass)
phenotypes$ROC <- log10(phenotypes$ROC)

### SPLIT THE HOMINOIDEA GROUP UP INTO HUMAN, CHIMP, AND REST ###
phenotypes_split <- phenotypes %>% mutate(Clade = ifelse(WilsonReederName == "Homo_sapiens" & Clade == "Hominoidea", "Hominoidea_Human", Clade)) %>% mutate(Clade = ifelse(WilsonReederName == "Pan_troglodytes_troglodytes" & Clade == "Hominoidea", "Hominoidea_Chimp", Clade))


### INTEROBSERVER RELIABILITY ###
# Make a separate table for evaluating the inter-observer reliability of the manual segmentations 
# As segmented by Neville Magielse (Cerebellum_lab) and Mahta Abbaspour (lobule_mahta) (see: https://colab.research.google.com/drive/1CkPf6SmP3MsSA7ovafPN6wewh3Nmtm6Q#scrollTo=84YQivOtOqF7)
# These data will be used to assess the influence of the observer via (phylogenetic) MANOVA, to assess if observer identity has a significant impact on volumes or on the interaction between volume and species.
Observeridentity <- Mahta %>%  left_join(Neville, by="SpecimenID") %>% left_join(lut, by="SpecimenID") %>% dplyr::rename(CrusVol_obs1=CrusVol.x) %>% dplyr::rename(CrusVol_obs2=CrusVol.y) %>% filter(CrusVol_obs1>"0.0000" & CrusVol_obs2>"0.0000")

#Log-transform these values, too
Observeridentity$CrusVol_obs1 <- log10(Observeridentity$CrusVol_obs1)
Observeridentity$CrusVol_obs2 <- log10(Observeridentity$CrusVol_obs2)
                                         
#Write the newly minted files into the output folder, ready for analysis.                                                                                            
write_csv2(phenotypes, "./0.cleanInput/Phenotypes.csv")
write_csv2(phenotypes_split, "./0.cleanInput/Phenotypes_split.csv")
write_csv2(Observeridentity, "./0.cleanInput/InterObserver_reliability.csv")


### SPECIES MEDIANS ###

## We also calculate species median values for cerebellar and cerebral volume, in order to calculate and plot species' raw measurement regressions or correlations.
# Cerebellar Volume, Cerebral Volume, Ansiform Area Volume
phenotypesMedian <- phenotypes_raw %>%
  group_by(WilsonReederName) %>%
  dplyr::summarize(MedianCerebellum = median(CerebellarVol, na.rm=TRUE)) %>% left_join(phenotypes_raw %>%
                                                                                      group_by(WilsonReederName) %>%
                                                                                      dplyr::summarize(MedianCerebrum = median(CerebralVol, na.rm=TRUE)), by = "WilsonReederName") %>% left_join(phenotypes_raw %>%
                                                                                                                                                                                                   group_by(WilsonReederName) %>%
                                                                                                                                                                                                   dplyr::summarize(MedianCrus = median(CrusVol_1_3, na.rm=TRUE)), by = "WilsonReederName") %>% left_join(phenotypes_raw %>%
                                                                                                                                                                                                                                                                                                            group_by(WilsonReederName) %>%
                                                                                                                                                                                                                                                                                                            dplyr::summarize(MedianROC = median(ROC, na.rm=TRUE)), by = "WilsonReederName") %>% left_join(BodyWeight, by = "WilsonReederName")



## Append species information
phenotypesMedian <- phenotypesMedian %>%
  left_join(lut, by="WilsonReederName") %>% dplyr::select(c("WilsonReederName","Hominoidea", "Clade", "MedianCerebellum", "MedianCerebrum", "MedianCrus", "MedianROC", "Body mass species mean")) %>% distinct()

## We also add columns for the ratios of Cerebellar/Cerebral volume, and Crus/Cerebellar Volume (as well as Crus/Rest of Cerebellum), using only species medians.
phenotypesMedian <- phenotypesMedian %>%
  mutate(MedianCerebroCerebellar = ((MedianCerebellum)/(MedianCerebrum))*100) %>% mutate(MedianCerebellarCrus = ((MedianCrus)/(MedianCerebellum))*100) %>% mutate(MedianROCCrus = ((MedianCrus/MedianROC))*100)

## We also add columns for volumes corrected for body mass.
phenotypesMedian <- phenotypesMedian  %>%
  mutate(MedianCerebellarMass = MedianCerebellum/`Body mass species mean`) %>% mutate(MedianCerebralMass = MedianCerebrum/`Body mass species mean`) %>% mutate(MedianCrusMass = MedianCrus/`Body mass species mean`)

phenotypesMedian_raw <- phenotypesMedian

## Write the median file (raw)
write_csv2(phenotypesMedian_raw, "./0.cleanInput/PhenotypesMedian_raw.csv")


## Log-transform the Median Volumes
phenotypesMedian$MedianCerebellum <- log10(phenotypesMedian$MedianCerebellum)
phenotypesMedian$MedianCerebrum <- log10(phenotypesMedian$MedianCerebrum)
phenotypesMedian$MedianCrus <- log10(phenotypesMedian$MedianCrus)
phenotypesMedian$MedianROC <- log10(phenotypesMedian$MedianROC)
phenotypesMedian$`Body mass species mean` <- log10(phenotypesMedian$`Body mass species mean`)
# Cerebellum/Cerebral and Ansiform Area/Cerebellum ratio do not need to be log-transformed, as we turned them into percentages earlier (in same order of magnitude).
phenotypesMedian$MedianCerebellarMass <- log10(phenotypesMedian$MedianCerebellarMass)
phenotypesMedian$MedianCerebralMass <- log10(phenotypesMedian$MedianCerebralMass)
phenotypesMedian$MedianCrusMass <- log10(phenotypesMedian$MedianCrusMass)

### SPLIT THE HOMINOIDEA GROUP UP INTO HUMAN, CHIMP, AND REST ###
phenotypesMedian_split <- phenotypesMedian %>% mutate(Clade = ifelse(WilsonReederName == "Homo_sapiens" & Clade == "Hominoidea", "Hominoidea_Human", Clade)) %>% mutate(Clade = ifelse(WilsonReederName == "Pan_troglodytes_troglodytes" & Clade == "Hominoidea", "Hominoidea_Chimp", Clade))

## Write the median file (log-transformed)
write_csv2(phenotypesMedian_split, "./0.cleanInput/PhenotypesMedian_split.csv")
write_csv2(phenotypesMedian, "./0.cleanInput/PhenotypesMedian.csv")



## Also for mean. This is only to compare to previous datasets that reported means only.
phenotypesMean <- phenotypes_raw %>%
  group_by(WilsonReederName) %>%
  dplyr::summarize(MeanCerebellum = mean(CerebellarVol, na.rm=TRUE)) %>% left_join(phenotypes_raw %>%
                                                                                     group_by(WilsonReederName) %>%
                                                                                     dplyr::summarize(MeanCerebrum = mean(CerebralVol, na.rm=TRUE)), by = "WilsonReederName") %>% left_join(phenotypes_raw %>%
                                                                                                                                                                                              group_by(WilsonReederName) %>%
                                                                                                                                                                                              dplyr::summarize(MeanCrus = mean(CrusVol_1_3, na.rm=TRUE)), by = "WilsonReederName") %>% left_join(phenotypes_raw %>%
                                                                                                                                                                                                                                                                                               group_by(WilsonReederName) %>%
                                                                                                                                                                                                                                                                                               dplyr::summarize(MeanROC = mean(ROC, na.rm=TRUE)), by = "WilsonReederName")
write_csv2(phenotypesMean, "./0.cleanInput/PhenotypesMean_raw.csv")


# Log-transform
phenotypesMean$MeanCerebellum <- log10(phenotypesMean$MeanCerebellum)
phenotypesMean$MeanCerebrum <- log10(phenotypesMean$MeanCerebrum)
phenotypesMean$MeanCrus <- log10(phenotypesMean$MeanCrus)
phenotypesMean$MeanROC <- log10(phenotypesMean$MeanROC)

write_csv2(phenotypesMean, "./0.cleanInput/PhenotypesMean.csv")

### ----------------------------- ###
# Note that 1 file 
# (PhenotypesMedian_raw_outliersRM_periods) 
# was added to the input directory
# manually. Make sure that it is
# also in your 0.cleanInput
# folder in order to proceed.
# i.e. clone the GitHub Repository
# folder structure and content.
### ----------------------------- ###

### ---------------------------- ###
### END of script
### -------------------------- ###                                                                                                                                
