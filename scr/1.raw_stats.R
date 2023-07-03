
### ------------------------------------------- ###
# Statistics and visualization of the raw 
# (i.e., non-phylogenetically corrected) data
# Roberto Toro, Katja Heuer 2018
# Script adapted by N. Magielse for inclusion of 
# cerebellar data and additional visualizations/ 
# normality checks - 2022
### ------------------------------------------- ###     

# Notes
# Enable corrplot in your 'Packages' tab if you use RStudio
# We are using R version 4.1.0 (2021-05-18)

### -------------------------- ###
# Set-up
### -------------------------- ###    

library(corrplot) # necessary?
library(nortest) # necessary?
library(tidyverse)
library(magrittr)
library(ggplot2)
library(nortest)

# Set working directory
## Set working directory
dir="/Documents/PhD/EvolutionPrimates/Final_code/" 
setwd(dir)

## Create outputdir
outputdir <- "1.raw_descriptives"
dir.create(outputdir)

### LOAD IN DATA ###
# Load neuroanatomical phenotypes
# Which should all have Wilson&Reeder names appended from running the '0.cleanInput' script.
pheno <- read_csv2("./0.cleanInput/Phenotypes.csv")
pheno_split <- read_csv2("./0.cleanInput/Phenotypes_split.csv")
pheno_raw <- read_csv2("./0.cleanInput/Phenotypes_raw.csv")
phenoMedian <- read_csv2("./0.cleanInput/PhenotypesMedian.csv")
phenoMedian_split <- read_csv2("./0.cleanInput/PhenotypesMedian_split.csv")
phenoMedian_raw <- read_csv2("./0.cleanInput/PhenotypesMedian_raw.csv")
lut <- read_csv2("./input/10kTrees_34PrimateSpecies_adapted.csv")  #look-up table


#phenoMedian_raw <- phenoMedian
#phenoMedian_raw$MedianCerebellum <- 10^phenoMedian_raw$MedianCerebellum
#phenoMedian_raw$MedianCerebrum <- 10^phenoMedian_raw$MedianCerebrum
#phenoMedian_raw$MedianCrus <- 10^phenoMedian_raw$MedianCrus


### -------------------------- ###
# Plotting set-up
### -------------------------- ###    

### Create plotting functions ###
# Plot matrix of scatter plots of the raw data.
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE, breaks=10)
  breaks <- h$breaks;
  nB <- length(breaks)
  y <- h$counts;
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

panel.cor <- function(x, y) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  test <- cor.test(x,y)
  rlow <- round(test$conf.int[1], digits=2)
  rup <- round(test$conf.int[2], digits=2)
  #p <- format.pval(test$p.value, eps=1e-6)
  txt <- paste0(r, "\n[", rlow, ", ",rup, "]")
  text(0.25, 0.25, txt, cex = 0.5)
}

loessCI <- function(x, y, myTitle, xlab, ylab, expx=FALSE, expy=FALSE) {
  plx <- predict(loess(y~x, span=0.7), se=TRUE)
  a <- order(x)
  if(expx) {
    x <- 10^x
  }
  plot(x,y,pch=19,col=col,cex=3)
  lines(x[a],plx$fit[a])
  lines(x[a],plx$fit[a]-qt(0.975,plx$df)*plx$se[a], lty=2)
  lines(x[a],plx$fit[a]+qt(0.975,plx$df)*plx$se[a], lty=2)
  title(myTitle, xlab = xlab, ylab = ylab)
}

loessCI_text <- function(x, y, myTitle, xlab, ylab, expx=FALSE, expy=FALSE) {
  plx <- predict(loess(y~x, span=0.7), se=TRUE)
  a <- order(x)
  if(expx) {
    x <- 10^x
  }
  plot(x,y,pch=19,col=col,cex=3)
  text(x,y,cex=0.2)
  lines(x[a],plx$fit[a])
  lines(x[a],plx$fit[a]-qt(0.975,plx$df)*plx$se[a], lty=2)
  lines(x[a],plx$fit[a]+qt(0.975,plx$df)*plx$se[a], lty=2)
  title(myTitle, xlab = xlab, ylab = ylab)
}

### Asssign colors to specimens ###

# Color the dots in the plot according to family
colour <- list()
colour["Galagonidae"] <- "#a5b7b77f"
colour["Loridae"] <- "#efe0ce7f"
colour["Lemuriformes"] <- "#b69b9d7f"
colour["Cebidae"] <- "#f4e9be7f"
colour["Atelidae"] <- "#caddef7f"
colour["Hominoidea"] <- "#e9c3927f"
colour["Hominoidea_Human"] <- "#d9c3927f"
colour["Hominoidea_Chimp"] <- "#ffc3927f"
colour["Colobinae"] <- "#c1cfb67f"
colour["Papionini"] <- "#939db77f"
colour["Cercopithecini"] <- "#ab91817f"
colour["red"] <- "#ff00007f"

# Add family info to each specimen
group <- list()
group["Human_51058"] <- "Hominoidea_Human"
group["Human_51060"] <- "Hominoidea_Human"
group["Human_51061"] <- "Hominoidea_Human"
group["Human_51062"] <- "Hominoidea_Human"
group["Human_51063"] <- "Hominoidea_Human"
group["Human_51114"] <- "Hominoidea_Human"
group["Human_51116"] <- "Hominoidea_Human"
group["Human_51118"] <- "Hominoidea_Human"
group["Human_51152"] <- "Hominoidea_Human"
group["Human_51153"] <- "Hominoidea_Human"
group["Chimp_Abby_Yerkes"] <- "Hominoidea_Chimp"
group["Chimp_Amanda_Yerkes"] <- "Hominoidea_Chimp"
group["ApelleCebusApell_f4b9"] <- "Cebidae"
group["Chimp_Artemus_Yerkes"] <- "Hominoidea_Chimp"
group["Chimp_Arthur_Yerkes"] <- "Hominoidea_Chimp"
group["Atele_6717"] <- "Atelidae"
group["AteleAtelesAter_7245"] <- "Atelidae"
group["AyeAyeDaubentoni_cb56"] <- "Lemuriformes"
group["Chimp_Azalea_Yerkes"] <- "Hominoidea_Chimp"
group["BabouinCercopith_a723"] <- "Papionini"
group["Chimp_Barbara_Yerkes"] <- "Hominoidea_Chimp"
group["Chimp_Barney_Yerkes"] <- "Hominoidea_Chimp"
group["Bonobo_brian"] <- "Hominoidea"
group["Chimp_Carl_Yerkes"] <- "Hominoidea_Chimp"
group["ColobeColobusPol_6118"] <- "Colobinae"
group["Chimp_David_Yerkes"] <- "Hominoidea_Chimp"
group["DouroucouliAotus_091c"] <- "Cebidae"
group["GalagoDemidoviiG_6cf3"] <- "Galagonidae"
group["Gibbon_buddy"] <- "Hominoidea"
group["Gorilla_kinyani"] <- "Hominoidea"
group["GorillaBeringeiG_0854"] <- "Hominoidea"
group["LangurIndienPith_8390"] <- "Colobinae"
group["LemurMongosLemur_b941"] <- "Lemuriformes"
group["LepilemurAQueueR_309f"] <- "Lemuriformes"
group["MacacaFascicular_3a21"] <- "Papionini"
group["MacacaFascicularis_32200"] <- "Papionini"
group["MacacaFascicularis_32202"] <- "Papionini"
group["MacacaFascicularis_32204"] <- "Papionini"
group["MacacaFascicularis_32205"] <- "Papionini"
group["MacacaFascicularis_mf01"] <- "Papionini"
group["MacacaFascicularis_mf02"] <- "Papionini"
group["MacacaMulatta_32198"] <- "Papionini"
group["MacacaMulatta_32199"] <- "Papionini"
group["MacacaMulatta_32201"] <- "Papionini"
group["MacacaMulatta_32203"] <- "Papionini"
group["MacacaMulatta_mm01"] <- "Papionini"
group["MacaqueCrabierMa_8f9e"] <- "Papionini"
group["MacaqueRhesusMac_bdf8"] <- "Papionini"
group["MakiCattaLemurCa_dec8"] <- "Lemuriformes"
group["MangabeyCercoceb_4458"] <- "Papionini"
group["MangabeyCouronne_6b65"] <- "Papionini"
group["MarsousetOedipom_e834"] <- "Cebidae"
group["MicrocebeDeCoque_e3a9"] <- "Lemuriformes"
group["MicrocebeMignonM_ca6c"] <- "Lemuriformes"
group["MoustacLasiopygi_23cc"] <- "Cercopithecini"
group["NycticebusTardig_08cb"] <- "Loridae"
group["Orangutan_80bb"] <- "Hominoidea"
group["OuistitiAPinceau_b4aa"] <- "Cebidae"
group["PithecusGermaini_6a2c"] <- "Colobinae"
group["SaimiriCassiquia_7176"] <- "Cebidae"
group["SaimiriCassiquia_eb44"] <- "Cebidae"
group["SapajouCapucinCe_d971"] <- "Cebidae"
group["SingeLaineuxLago_0f84"] <- "Atelidae"
group["VariNoirEtBlancL_4d0e"] <- "Lemuriformes"
group["VervetCercopithe_72ae"] <- "Cercopithecini"

col=unlist(colour[unlist(group[pheno_split$SpecimenID])])
col_ape=col[grepl("Hominoidea",names(col))]
col_nonape <- col[-pmatch(col_ape,col)]


### ---------------------------- ###
# Raw descriptives - all specimens
### ---------------------------- ###    

## Save figure of scatter plots
pdf('1.raw_descriptives/1.1.scatterplots.pdf')
pairs(~CerebellarVol+CerebralVol+CrusVol_1_3+SurfaceArea+AbsGI+FoldingLength+FoldingNumber+Lambda+Delta,
      data=pheno_raw,
      labels=c("Cerebellar \n Volume","Cerebral \n Volume","Ansiform \n Area Volume","Cerebral \n Surface Area","Absolute \n Gyrification","Folding \n Length", "Folding \n Number", "Fold \n Wavelength","Fold \nDepth"),
      panel=panel.smooth,
      cex.labels = 0.6,
      upper.panel = NULL)
graphics.off()

## Plot & save correlation matrix
M <- cbind(pheno_raw$CerebellarVol, pheno_raw$CerebralVol, pheno_raw$CrusVol_1_3, pheno_raw$SurfaceArea, pheno_raw$AbsGI,pheno_raw$FoldingLength,pheno_raw$FoldingNumber,pheno_raw$Lambda,pheno_raw$Delta)
R <- cor(na.omit(M))
rownames(R) <- colnames(R) <- c("Cerebellar \n Volume","Cerebral \n Volume","Ansiform \ Area","Cerebral \n Surface Area","Absolute \n Gyrification","Folding \n Length", "Folding \n Number", "Fold \n Wavelength","Fold \nDepth")
pdf('1.raw_descriptives/1.2.correlationsEllipse.pdf')
corrplot(R, method="ellipse", order="AOE")
dev.off()

# With original data order
M <- cbind(pheno_raw$CerebellarVol, pheno_raw$CerebralVol, pheno_raw$CrusVol_1_3, pheno_raw$SurfaceArea, pheno_raw$AbsGI,pheno_raw$FoldingLength,pheno_raw$FoldingNumber,pheno_raw$Lambda,pheno_raw$Delta)
R <- cor(na.omit(M))
rownames(R) <- colnames(R) <- c("Cerebellar \n Volume","Cerebral \n Volume","Ansiform Area \n Volume","Cerebral \n Surface Area","Absolute \n Gyrification","Folding \n Length", "Folding \n Number", "Fold \n Wavelength","Fold \nDepth")
pdf('1.raw_descriptives/1.2a.correlationsEllipse.pdf')
corrplot(R, method="ellipse", order="original")
dev.off()

# With numerical correlation representation
pdf('1.raw_descriptives/1.3.correlationsNumerical.pdf')
corrplot(R, method="number", order="AOE")
dev.off()

# And with original data order
pdf('1.raw_descriptives/1.3a.correlationsNumerical.pdf')
corrplot(R, method="number", order="original",
         tl.cex = 0.6,
         type = "upper",
         tl.col = "black",
         diag = FALSE,
         col = colorRampPalette(c("darkblue", "white", "darkred"))(100),
         tl.srt = 45,
         tl.pos = "n",
         tl.offset = 1.0)
dev.off()

## Create vectors to test normality of the data. 
C <- pheno_raw$CerebellarVol 
names(C) <- pheno_raw$SpecimenID
V <- pheno_raw$CerebralVol
names(V) <- pheno_raw$SpecimenID
lillie.test(C) # not normal (as expected under assumption of phylogenetic signal in traits)
lillie.test(V) # not normal

## Select only those specimens that have complete data (based on AA segmentation availability)
pheno_raw_omit <- pheno_raw %>% dplyr::select("SpecimenID", "CerebellarVol", "CerebralVol", "CrusVol_1_3", "ROC") %>% na.omit()
C_omit <- pheno_raw_omit$CerebellarVol
names(C_omit) <- pheno_raw_omit$SpecimenID
V_omit <- pheno_raw_omit$CerebralVol
names(V_omit) <- pheno_raw_omit$SpecimenID
AA_omit <- pheno_raw_omit$CrusVol_1_3
names(AA_omit) <- pheno_raw_omit$SpecimenID

## Perform LOESS (non-parametric, does not assume data-normality) between main variables of interest in the current study.

## Save figures
# Cerebellum vs Cerebrum
pdf('1.raw_descriptives/1.4.cerebellar-cerebral.pdf')
loessCI(V, C, "Cerebellar Volume versus Cerebral Volume", xlab = "Cerebral Volume", ylab = "Cerebellar Volume")
graphics.off()

pdf('1.raw_descriptives/1.4a.cerebellar-cerebral.pdf')
loessCI_text(V, C, "Cerebellar Volume versus Cerebral Volume", xlab = "Cerebral Volume", ylab = "Cerebellar Volume")
graphics.off()

# Ansiform Area vs Cerebellum
pdf('1.raw_descriptives/1.5.crus-cerebellar.pdf')
loessCI(C_omit, AA_omit, "Ansiform Area Volume versus Cerebellar Volume", xlab = "Cerebellar Volume", ylab = "Ansiform Area Volume")
graphics.off()

pdf('1.raw_descriptives/1.5a.crus-cerebellar.pdf')
loessCI_text(C_omit, AA_omit, "Ansiform Area Volume versus Cerebellar Volume", xlab = "Cerebellar Volume", ylab = "Ansiform Area Volume")
graphics.off()

#Ansiform Area vs Cerebrum
pdf('1.raw_descriptives/1.6.crus-cerebral.pdf')
loessCI(V_omit, AA_omit, "Ansiform Area Volume versus Cerebral Volume", xlab = "Cerebral Volume", ylab = "Ansiform Area Volume")
graphics.off()

pdf('1.raw_descriptives/1.6a.crus-cerebral.pdf')
loessCI_text(V_omit, AA_omit, "Ansiform Area Volume versus Cerebral Volume", xlab = "Cerebral Volume", ylab = "Ansiform Area Volume")
graphics.off()


### -------------------------- ###
# Log-transformed data
### -------------------------- ###    

## Save figure of scatter plots
pdf('1.raw_descriptives/1.7.logscatterplots.pdf')
pairs(~CerebellarVol+CerebralVol+CrusVol_1_3+SurfaceArea+AbsGI+FoldingLength+FoldingNumber+Lambda+Delta,
      data=pheno,
      labels=c("Cerebellar \n Volume","Cerebral \n Volume","Ansiform \n Area Volume","Cerebral \n Surface Area","Absolute \n Gyrification","Folding \n Length", "Folding \n Number", "Fold \n Wavelength","Fold \nDepth"),
      panel=panel.smooth,
      cex.labels = 0.6,
      upper.panel=NULL)
graphics.off()

## Plot & save correlation matrix
Mlog <- cbind(pheno$CerebellarVol, pheno$CerebralVol, pheno$CrusVol_1_3, pheno$SurfaceArea, pheno$AbsGI,pheno$FoldingLength,pheno$FoldingNumber,pheno$Lambda,pheno$Delta)
Rlog <- cor(na.omit(Mlog))
rownames(Rlog) <- colnames(Rlog) <- c("Cerebellar \n Volume","Cerebral \n Volume","Ansiform \ Area","Cerebral \n Surface Area","Absolute \n Gyrification","Folding \n Length", "Folding \n Number", "Fold \n Wavelength","Fold \nDepth")
pdf('1.raw_descriptives/1.8.correlationsEllipse.pdf')
corrplot(R, method="ellipse", order="AOE")
dev.off()

# With original data order
pdf('1.raw_descriptives/1.8a.correlationsEllipse.pdf')
corrplot(Rlog, method="ellipse", order="original")
dev.off()

# With numerical correlation representation
pdf('1.raw_descriptives/1.9.correlationsNumerical.pdf')
corrplot(R, method="number", order="AOE")
dev.off()

# And with original data order
pdf('1.raw_descriptives/1.9a.correlationsNumerical.pdf')
corrplot(R, method="number", order="original",
         tl.cex = 0.6,
         type = "upper",
         tl.col = "black",
         tl.srt = 45,
         tl.offset = 1.0,
         tl.pos = "n",
         diag = FALSE)
dev.off()

## Create logged cerebellar and cerebral volume vectors.
Clog <- pheno$CerebellarVol
names(Clog) <- pheno$SpecimenID
Vlog <- pheno$CerebralVol
names(Vlog) <- pheno$SpecimenID
CRlog <- pheno$ROC
names(CRlog) <- pheno$SpecimenID
lillie.test(Clog) # not normal
lillie.test(Vlog) # not normal

## Select only those specimens that have complete data (based on ansiform area segmentation availability i.e., 13 species and 30 specimens)
pheno_omit <- pheno %>% dplyr::select("SpecimenID", "CerebellarVol", "CerebralVol", "CrusVol_1_3", "ROC") %>% na.omit()
Clog_omit <- pheno_omit$CerebellarVol
Vlog_omit <- pheno_omit$CerebralVol
AAlog_omit <- pheno_omit$CrusVol_1_3
CRlog_omit <- pheno_omit$ROC

## Save figures
# Cerebellum vs Cerebrum
pdf('1.raw_descriptives/1.10.cerebellar-cerebral.pdf')
loessCI(Vlog, Clog, "Cerebellar Volume versus Cerebral Volume", xlab = "Cerebral Volume", ylab = "Cerebellar Volume")
graphics.off()

pdf('1.raw_descriptives/1.10a.cerebellar-cerebral.pdf')
loessCI_text(Vlog, Clog, "Cerebellar Volume versus Cerebral Volume", xlab = "Cerebral Volume", ylab = "Cerebellar Volume")
graphics.off()

# Ansiform Area vs Cerebellum
pdf('1.raw_descriptives/1.11.crus-cerebellar.pdf')
loessCI(Clog_omit, AAlog_omit, "Ansiform Area Volume versus Cerebellar Volume", xlab = "Cerebellar Volume", ylab = "Ansiform Area Volume")
graphics.off()

pdf('1.raw_descriptives/1.11a.crus-cerebellar.pdf')
loessCI_text(Clog_omit, AAlog_omit, "Ansiform Area Volume versus Cerebellar Volume", xlab = "Cerebellar Volume", ylab = "Ansiform Area Volume")
graphics.off()

# Ansiform Area vs Rest of cerebellum
pdf('1.raw_descriptives/1.11b.crus-ROC.pdf')
loessCI(CRlog_omit, AAlog_omit, "Ansiform Area Volume versus Rest of Cerebellar Volume", xlab = "Rest of Cerebellar Volume", ylab = "Ansiform Area Volume")
graphics.off()

pdf('1.raw_descriptives/1.11c.crus-ROC.pdf')
loessCI_text(CRlog_omit, AAlog_omit, "Ansiform Area Volume versus Rest of Cerebellar Volume", xlab = "Rest of Cerebellar Volume", ylab = "Ansiform Area Volume")
graphics.off()

#Ansiform Area vs Cerebrum
pdf('1.raw_descriptives/1.12.crus-cerebral.pdf')
loessCI(Vlog_omit, AAlog_omit, "Ansiform Area Volume versus Cerebral Volume", xlab = "Cerebral Volume", ylab = "Ansiform Area Volume")
graphics.off()

pdf('1.raw_descriptives/1.12a.crus-cerebral.pdf')
loessCI_text(Vlog_omit, AAlog_omit, "Ansiform Area Volume versus Cerebral Volume", xlab = "Cerebral Volume", ylab = "Ansiform Area Volume")
graphics.off()
  

### ------------------------------- ###
# Species median data, log-transformed
### ------------------------------- ###    

## Perform LOESS between specific variables for logged phenotypes
Cmed <- phenoMedian$MedianCerebellum
names(Cmed) <- phenoMedian$WilsonReederName # No more SpecimenID column, since every species is represented by its median specimen (based on cerebellar volume).
Vmed <- phenoMedian$MedianCerebrum
names(Vmed) <- phenoMedian$WilsonReederName
CRmed <- phenoMedian$MedianROC
names(CRmed) <- phenoMedian$WilsonReederName
lillie.test(Clog) # still not normal
lillie.test(Vlog) # not normal

## Select only those specimens that have complete data (based on ansiform area segmentation availability; 13 species & 30 specimens)
phenoMedian_omit <- phenoMedian %>% dplyr::select("WilsonReederName", "MedianCerebellum", "MedianCerebrum", "MedianCrus","MedianROC") %>% na.omit()
Cmed_omit <- phenoMedian_omit$MedianCerebellum
Vmed_omit <- phenoMedian_omit$MedianCerebrum
AAmed_omit <- phenoMedian_omit$MedianCrus
CRmed_omit <- phenoMedian_omit$MedianROC

## Save figures
# Cerebellum vs Cerebrum
pdf('1.raw_descriptives/1.13.cerebellar-cerebral.pdf')
loessCI(Vmed, Cmed, "Cerebellar Volume versus Cerebral Volume", xlab = "Cerebral Volume", ylab = "Cerebellar Volume")
graphics.off()

pdf('1.raw_descriptives/1.13a.cerebellar-cerebral.pdf')
loessCI_text(Vmed, Cmed, "Cerebellar Volume versus Cerebral Volume", xlab = "Cerebral Volume", ylab = "Cerebellar Volume")
graphics.off()

# Ansiform Area vs Cerebellum
pdf('1.raw_descriptives/1.14.crus-cerebellar.pdf')
loessCI(Cmed_omit, AAmed_omit, "Ansiform Area Volume versus Cerebellar Volume", xlab = "Cerebellar Volume", ylab = "Ansiform Area Volume")
graphics.off()

pdf('1.raw_descriptives/1.14a.crus-cerebellar.pdf')
loessCI_text(Cmed_omit, AAmed_omit, "Ansiform Area Volume versus Cerebellar Volume", xlab = "Cerebellar Volume", ylab = "Ansiform Area Volume")
graphics.off()

# Ansiform Area vs Rest of Cerebellum
pdf('1.raw_descriptives/1.14b.crus-ROC.pdf')
loessCI(CRmed_omit, AAmed_omit, "Ansiform Area Volume versus Rest of Cerebellar Volume", xlab = "Rest of Cerebellar Volume", ylab = "Ansiform Area Volume")
graphics.off()

pdf('1.raw_descriptives/1.14c.crus-ROC.pdf')
loessCI_text(CRmed_omit, AAmed_omit, "Ansiform Area Volume versus Rest of Cerebellar Volume", xlab = "Rest of Cerebellar Volume", ylab = "Ansiform Area Volume")
graphics.off()

#Ansiform Area vs Cerebrum
pdf('1.raw_descriptives/1.15.crus-cerebral.pdf')
loessCI(Vmed_omit, AAmed_omit, "Ansiform Area Volume versus Cerebral Volume", xlab = "Cerebral Volume", ylab = "Ansiform Area Volume")
graphics.off()

pdf('1.raw_descriptives/1.15a.crus-cerebral.pdf')
loessCI_text(Vlog_omit, AAlog_omit, "Ansiform Area Volume versus Cerebral Volume", xlab = "Cerebral Volume", ylab = "Ansiform Area Volume")
graphics.off()


### -------------------------- ###
# Basis for table 1.
### -------------------------- ###     

library(sjPlot)
library(dplyr)
library(tidyr)

# We start by selecting the variables of interest from our full phenotype table.
# We will also rename all variables to be more descriptive.
pheno_table1 <- pheno %>% dplyr::select("WilsonReederName", "CerebellarVol", "CerebralVol","CrusVol_1_3", "CerebroCerebellar", "CerebellarCrus") %>% dplyr::rename(Species = WilsonReederName) %>% dplyr::rename('Cerebellar Volume' = CerebellarVol) %>% dplyr::rename('Cerebral Volume' = CerebralVol) %>% dplyr::rename('Ansiform Area Volume' = CrusVol_1_3) %>% dplyr::rename('Cerebellum-to-Cerebrum Ratio' = CerebroCerebellar) %>% dplyr::rename('Ansiform Area-to-Cerebellum Ratio' = CerebellarCrus)
# Same thing for the raw phenotypes
pheno_table1a <- pheno_raw %>% dplyr::select("WilsonReederName", "CerebellarVol", "CerebralVol","CrusVol_1_3", "CerebroCerebellar", "CerebellarCrus") %>% dplyr::rename(Species = WilsonReederName) %>% dplyr::rename('Cerebellar Volume' = CerebellarVol) %>% dplyr::rename('Cerebral Volume' = CerebralVol) %>% dplyr::rename('Ansiform Area Volume' = CrusVol_1_3) %>% dplyr::rename('Cerebellum-to-Cerebrum Ratio' = CerebroCerebellar) %>% dplyr::rename('Ansiform Area-to-Cerebellum Ratio' = CerebellarCrus)

#Take the median value for all the main traits of interest in one go.
pheno_table1_med <- pheno_table1 %>% group_by(Species) %>%  
  summarise_all(.funs = median)
#Same for the raw phenotypes
pheno_table1a_med <- pheno_table1a %>% group_by(Species) %>%  
  summarise_all(.funs = median)

# Plot the publication-ready table. (Supposedly; we performed touching-up in Microsoft Word and Inkscape after)
tab_df(pheno_table1a_med, file="1.raw_descriptives/1.16.phenotypes_raw_median_alphabetical.doc") # save this as a Word document

# Order the data based on the tree
tree <- force.ultrametric(read.nexus("./input/consensusTree_10kTrees_Primates_Version3_34.nex"), method="extend")
plotTree(tree) # If you want to visualize your tree to see if the data matches it.
orderphylogeny <- c("Galago_demidoff", "Loris_tardigradus", "Daubentonia_madagascariensis", "Varecia_variegata_variegata", "Lemur_catta", "Eulemur_mongoz", "Lepilemur_ruficaudatus", "Mirza_coquereli", "Microcebus_murinus", "Saimiri_sciureus", "Cebus_apella", "Cebus_capucinus", "Saguinus_oedipus", "Callithrix_penicillata", "Aotus_trivirgatus", "Lagothrix_lagotricha", "Ateles_paniscus", "Pongo_pygmaeus", "Pan_troglodytes_troglodytes", "Pan_paniscus", "Homo_sapiens", "Gorilla_gorilla_gorilla", "Gorilla_beringei_graueri", "Hylobates_lar", "Trachypithecus_(Trachypithecus)_germaini", "Semnopithecus_entellus", "Colobus_polykomos", "Macaca_mulatta", "Macaca_fascicularis", "Papio_hamadryas", "Lophocebus_albigena", "Cercocebus_atys", "Chlorocebus_sabaeus", "Cercopithecus_cephus_cephus")
# There are functions to do this as well; here it was done manually.
pheno_table1a_med <- pheno_table1a_med[match(orderphylogeny, pheno_table1a_med$Species),]

tab_df(pheno_table1a_med, file="1.raw_descriptives/1.16a.phenotypes_raw_median_tree.doc", alternate.rows = T, show.footnote = T, title = "Neuroanatomical measurements", footnote = "Neuroanatomical measurements collected in the current study. Species median measurements for previously reported cerebral volumes (Heuer et al., 2019) are reported alongside cerebellar and ansiform area volumes. Species median ratios between cerebellar and cerebral volume, and ansiform area and cerebellar volume are also given. Species are ordered by the phylogenetic tree and clades are colored to correspond to Heuer et al. 2019. For four species (Homo sapiens, Pan troglodytes, Macaca fascicularis, and Macaca mulatta) several specimens were available. For these species, median absolute deviations are also reported.")

### -------------------------- ###
#  Median abolute deviations.
### -------------------------- ###     

hom_sap <- pheno_table1a[pheno_table1a$Species == 'Homo_sapiens',] #10
pan_trogtrog <- pheno_table1a[pheno_table1a$Species == 'Pan_troglodytes_troglodytes',] #9
mac_fasc <- pheno_table1a[pheno_table1a$Species == 'Macaca_fascicularis',] #8
mac_mul <- pheno_table1a[pheno_table1a$Species == 'Macaca_mulatta',] #6
# Perhaps less informative: since only two specimens. Included for now. Were ultimated not mentioned in the manuscript.
ate_pan <- pheno_table1a[pheno_table1a$Species == 'Ateles_paniscus',] #2
sai_sci <- pheno_table1a[pheno_table1a$Species == 'Saimiri_sciureus',] #2

# Calculate the Median Absolute Deviations. Manually added to the median observation table for publication.
hs_mad_cerebellum <- mad(hom_sap$`Cerebellar Volume`)
hs_mad_cerebrum <- mad(hom_sap$`Cerebral Volume`)
hs_mad_ansiform <- mad(hom_sap$`Ansiform Area Volume`)
hs_mad_CC_ratio <- mad(hom_sap$`Cerebellum-to-Cerebrum Ratio`)
hs_mad_AAC_ratio <- mad(hom_sap$`Ansiform Area-to-Cerebellum Ratio`)

ptt_mad_cerebellum <- mad(pan_trogtrog$`Cerebellar Volume`)
ptt_mad_cerebrum <- mad(pan_trogtrog$`Cerebral Volume`)
ptt_mad_ansiform <- mad(pan_trogtrog$`Ansiform Area Volume`)
ptt_mad_CC_ratio <- mad(pan_trogtrog$`Cerebellum-to-Cerebrum Ratio`)
ptt_mad_AAC_ratio <- mad(pan_trogtrog$`Ansiform Area-to-Cerebellum Ratio`)

mf_mad_cerebellum <- mad(mac_fasc$`Cerebellar Volume`)
mf_mad_cerebrum <- mad(mac_fasc$`Cerebral Volume`)
mf_mad_CC_ratio <- mad(mac_fasc$`Cerebellum-to-Cerebrum Ratio`)

mm_mad_cerebellum <- mad(mac_mul$`Cerebellar Volume`)
mm_mad_cerebrum <- mad(mac_mul$`Cerebral Volume`)
mm_mad_CC_ratio <- mad(mac_mul$`Cerebellum-to-Cerebrum Ratio`)

#Less informative, low n.
ap_mad_cerebellum <- mad(ate_pan$`Cerebellar Volume`)
ap_mad_cerebrum <- mad(ate_pan$`Cerebral Volume`)
ap_mad_CC_ratio <- mad(ate_pan$`Cerebellum-to-Cerebrum Ratio`)

ss_mad_cerebellum <- mad(sai_sci$`Cerebellar Volume`)
ss_mad_cerebrum <- mad(sai_sci$`Cerebral Volume`)
ss_mad_CC_ratio <- mad(sai_sci$`Cerebellum-to-Cerebrum Ratio`)

sink("1.raw_descriptives/1.17.Median_absolute_deviations.txt")
print("Human: Cerebellar Volume - Cerebral Volume - Ansiform Area Volume - Cerebro-Cerebellar Ratio - Ansiform Area-Cerebellar Ratio")
hs_mad_cerebellum
hs_mad_cerebrum
hs_mad_ansiform
hs_mad_CC_ratio
hs_mad_AAC_ratio

print("Chimpanzee: Cerebellar Volume - Cerebral Volume - Ansiform Area Volume - Cerebro-Cerebellar Ratio - Ansiform Area-Cerebellar Ratio")
ptt_mad_cerebellum
ptt_mad_cerebrum
ptt_mad_ansiform
ptt_mad_CC_ratio
ptt_mad_AAC_ratio

print("Crab-eating macaque: Cerebellar Volume - Cerebral Volume - Cerebro-Cerebellar Ratio")
mf_mad_cerebellum
mf_mad_cerebrum
mf_mad_CC_ratio

print("Rhesus macaque: Cerebellar Volume - Cerebral Volume - Cerebro-Cerebellar Ratio")
mm_mad_cerebellum
mm_mad_cerebrum
mm_mad_CC_ratio

# Less informative?
print("Black spider monkey: Cerebellar Volume - Cerebral Volume - Cerebro-Cerebellar Ratio")
ap_mad_cerebellum
ap_mad_cerebrum
ap_mad_CC_ratio

print("Common squirrel monkey: Cerebellar Volume - Cerebral Volume - Cerebro-Cerebellar Ratio")
ss_mad_cerebellum
ss_mad_cerebrum
ss_mad_CC_ratio

sink()


### -------------------------- ###
# Normality tests for all species with several specimens.
### -------------------------- ###  


sink("1.raw_descriptives/1.18.shapiro-tests_multipleObservations.txt")
print("Homo sapiens")
lapply(hom_sap[,2:6], shapiro.test)
print("Pan troglodytes")
lapply(pan_trogtrog[,2:6], shapiro.test)
print("Macaca fascicularis")
lapply(mac_fasc[,c(2,3,5)], shapiro.test)
print("Macaca mulatta")
lapply(mac_mul[,c(2,3,5)], shapiro.test)
sink()

## For Macaca mulatta, cerebellar and cerebral volumes are not normally distributed.
## We create box plots for these species, specifically.

## To be able to label the outliers, we create a function to identify outliers ##
find_outlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}


### ---------------------------- ###
# Box plots, highlighting outliers
### ---------------------------- ###    

## Load packages
library(ggpubr)
library(viridis)
library(hrbrthemes)
#hrbrthemes::import_roboto_condensed()

# Create a viridis 6-color scale vector. This allows us to color-match the box plots which can only represent two species (those measuring ansiform area volumes).
sixcol <- viridis(6, alpha = 0.6, begin = 0, end = 1, direction = 1, option = "D")
names(sixcol) <- c('Black spider monkey', 'Central chimpanzee', 'Common squirrel monkey', 'Crab-eating macaque', 'Human', 'Rhesus monkey')

## Crab-eating macaque ##
mf_box_cerebell <- pheno[pheno$`English Name` %in% c('Crab-eating macaque'), ] %>%
  group_by(`English Name`) %>%
  mutate(outlier = ifelse(find_outlier(CerebellarVol), SpecimenID , NA)) %>%
  ggplot(aes(y= CerebellarVol, x=`English Name`, fill=`English Name`)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = "triangle", outlier.size = 3) +
  scale_fill_manual(values = sixcol) +
  geom_point(color="black", size=1.4, alpha=0.9) +
  theme_ipsum(base_size = 10) +
  theme(legend.position="top", legend.text = element_text(size = 10), legend.key.size = unit(1, 'cm'),
        plot.title = element_text(size=16, face = "bold", hjust = 0.5), axis.title.y = element_text(hjust=0.5, size = 8, margin = margin(r = 40)), axis.text.x = element_text(size = 10)
        )  +
  ggtitle("Cerebellar volume distribution") +
  labs(y="Cerebellar volumes (in mm3, log-transformed)", x="") +
  geom_text(aes(label=outlier), na.rm=TRUE, hjust=-.5, size = 2) ## Mark the outlier
mf_box_cerebell

mf_box_cerebr <- pheno[pheno$`English Name` %in% c('Crab-eating macaque'), ] %>%
  group_by(`English Name`) %>%
  mutate(outlier = ifelse(find_outlier(CerebellarVol), SpecimenID , NA)) %>%
  ggplot(aes(y= CerebralVol, x=`English Name`, fill=`English Name`)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = "triangle", outlier.size = 3) +
  scale_fill_manual(values = sixcol) +
  geom_point(color="black", size=1.4, alpha=0.9) +
  theme_ipsum(base_size = 10) +
  theme(legend.position="top", legend.text = element_text(size = 10), legend.key.size = unit(1, 'cm'),
        plot.title = element_text(size=16, face = "bold", hjust = 0.5), axis.title.y = element_text(hjust=0.5, size = 8, margin = margin(r = 40)), axis.text.x = element_text(size = 10)
  )  +
  ggtitle("Cerebral volume distribution") +
  labs(y="Cerebral volumes (in mm3, log-transformed)", x="") +
  geom_text(aes(label=outlier), na.rm=TRUE, hjust=-.5, size = 2) ## Mark the outlier
mf_box_cerebr

## Rhesus monkey ##
mm_box_cerebell <- pheno[pheno$`English Name` %in% c('Rhesus monkey'), ] %>%
  group_by(`English Name`) %>%
  mutate(outlier = ifelse(find_outlier(CerebellarVol), SpecimenID , NA)) %>%
  ggplot(aes(y= CerebellarVol, x=`English Name`, fill=`English Name`)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = "triangle", outlier.size = 3) +
  scale_fill_manual(values = sixcol) +
  geom_point(color="black", size=1.4, alpha=0.9) +
  theme_ipsum(base_size = 10) +
  theme(legend.position="top", legend.text = element_text(size = 10), legend.key.size = unit(1, 'cm'),
        plot.title = element_text(size=16, face = "bold", hjust = 0.5), axis.title.y = element_text(hjust=0.5, size = 8, margin = margin(r = 40)), axis.text.x = element_text(size = 10)
  )  +
  labs(y="Cerebellar volume (in mm3, log-transformed)", x="") +
  geom_text(aes(label=outlier), na.rm=TRUE, hjust=-.5, size = 2) ## Mark the outlier
mm_box_cerebell

mm_box_cerebr <- pheno[pheno$`English Name` %in% c('Rhesus monkey'), ] %>%
  group_by(`English Name`) %>%
  mutate(outlier = ifelse(find_outlier(CerebellarVol), SpecimenID , NA)) %>%
  ggplot(aes(y= CerebralVol, x=`English Name`, fill=`English Name`)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = "triangle", outlier.size = 3) +
  scale_fill_manual(values = sixcol) +
  geom_point(color="black", size=1.4, alpha=0.9) +
  theme_ipsum(base_size = 10) +
  theme(legend.position="top", legend.text = element_text(size = 10), legend.key.size = unit(1, 'cm'),
        plot.title = element_text(size=16, face = "bold", hjust = 0.5), axis.title.y = element_text(hjust=0.5, size = 8, margin = margin(r = 40)), axis.text.x = element_text(size = 10)
  )  +
  labs(y="Cerebral volume (in mm3, log-transformed)", x="") +
  geom_text(aes(label=outlier), na.rm=TRUE, hjust=-.5, size = 2) ## Mark the outlier
mm_box_cerebr

### ------------------------------------------------------------ ###
# Create a supplemental figure with boxplots, and outliers marked.
### ------------------------------------------------------------ ###    
require(grid)   # for the textGrob() function
tiff("./1.raw_descriptives/1.19.CE-Rhesus_outliers.tiff", width = 8, height = 6, units = 'in', res = 300)
figure_suppl1 <- ggarrange(mf_box_cerebell + rremove("ylab") + rremove("xlab"), mf_box_cerebr + rremove("ylab") + rremove("xlab"), mm_box_cerebell + rremove("ylab") + rremove("xlab"), mm_box_cerebr + rremove("ylab") + rremove("xlab"),
                           labels = c("a", "b", "c", "d"),
                           ncol = 2, nrow = 2, common.legend = F) # Common legend does not give color for Rhesus monkey.
figure_suppl1 <- annotate_figure(figure_suppl1, left = textGrob("Volume (log-transformed)", rot = 90, vjust = 0.35, hjust = 0.5, gp = gpar(cex = 1.0)))
figure_suppl1 
graphics.off()


## Cerebellar and cerebral outliers belong to: MacaqueRhesusMac_bdf8 & 	MacaqueCrabierMa_8f9e

### ---------------------------------------- ###
# Box plots for all intraspecific variability.
### ---------------------------------------- ###    

### ---------------------------------------- ###
# Absolute measurements.
### ---------------------------------------- ###    

## Cerebellar Volume
box_cerebell <- pheno[pheno$`English Name` %in% c('Human', 'Central chimpanzee', 'Crab-eating macaque', 'Rhesus monkey', 'Black spider monkey', 'Common squirrel monkey'), ] %>%
  group_by(`English Name`) %>%
  mutate(outlier = ifelse(find_outlier(CerebellarVol), SpecimenID , NA)) %>%
  ggplot(aes(y= CerebellarVol, x=`English Name`, fill=`English Name`)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = "triangle", outlier.size = 3) +
  scale_fill_manual(values = sixcol) +
  geom_point(color="black", size=1.4, alpha=0.9) +
  theme_ipsum(base_size = 14) +
  theme(legend.position="top", legend.text = element_text(), legend.key.size = unit(3, 'cm'),
        plot.title = element_text(size=30, face = "bold", hjust = 0.5), axis.title.y = element_text(hjust=0.5, size = 15, margin = margin(r = 40))
  ) +
  ggtitle("Cerebellar volumes distribution") +
  labs(y="Cerebellar volumes (in mm3, log-transformed)", x="") +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier
box_cerebell

## Cerebral Volume
box_cerebr <- pheno[pheno$`English Name` %in% c('Human', 'Central chimpanzee', 'Crab-eating macaque', 'Rhesus monkey', 'Black spider monkey', 'Common squirrel monkey'), ] %>%
  group_by(`English Name`) %>%
  mutate(outlier = ifelse(find_outlier(CerebellarVol), SpecimenID , NA)) %>%
  ggplot(aes(y= CerebralVol, x=`English Name`, fill=`English Name`)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = "triangle", outlier.size = 3) +
  scale_fill_manual(values = sixcol) +
  geom_point(color="black", size=1.4, alpha=0.9) +
  theme_ipsum(base_size = 14) +
  theme(legend.position="top", legend.text = element_text(), legend.key.size = unit(3, 'cm'),
        plot.title = element_text(size=30, face = "bold", hjust = 0.5), axis.title.y = element_text(hjust=0.5, size = 15, margin = margin(r = 40))
  ) +
  ggtitle("Cerebral volumes distribution") +
  labs(y="Cerebral volumes (in mm3, log-transformed)", x="") +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier
box_cerebr

## Ansiform Area Volume
box_aa <- pheno[pheno$`English Name` %in% c('Human', 'Central chimpanzee'), ] %>%
  group_by(`English Name`) %>%
  mutate(outlier = ifelse(find_outlier(CrusVol_1_3), SpecimenID , NA)) %>%
  ggplot(aes(y= CrusVol_1_3, x=`English Name`, fill=`English Name`)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = "triangle", outlier.size = 3) +
  scale_fill_manual(values = sixcol) +
  geom_point(color="black", size=1.4, alpha=0.9) +
  theme_ipsum(base_size = 14) +
  theme(legend.position="top", legend.text = element_text(), legend.key.size = unit(3, 'cm'),
        plot.title = element_text(size=30, face = "bold", hjust = 0.5), axis.title.y = element_text(hjust=0.5, size = 15, margin = margin(r = 40))
  ) +
  ggtitle("Ansiform area volumes distribution") +
  labs(y="Ansiform area volumes (in mm3, log-transformed)", x="") +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier
box_aa

## Combining cerebellar, cerebral, and ansiform area volumes in a single plot.
sixcol2 <- viridis(6, alpha = 0.6, begin = 0, end = 1, direction = 1, option = "D") %>% rep(2) 
sixcol2 <- c(sixcol2, "#41448799", '#7AD15199')   
names(sixcol2) <- c('Black spider monkey.Cerebellum', 'Central chimpanzee.Cerebellum', 'Common squirrel monkey.Cerebellum', 'Crab-eating macaque.Cerebellum', 'Human.Cerebellum', 'Rhesus monkey.Cerebellum', 'Black spider monkey.Cerebrum', 'Central chimpanzee.Cerebrum', 'Common squirrel monkey.Cerebrum', 'Crab-eating macaque.Cerebrum', 'Human.Cerebrum', 'Rhesus monkey.Cerebrum', 'Central chimpanzee.Ansiform area', 'Human.Ansiform area')

library(ggpattern)

tiff("./1.raw_descriptives/1.20. spread_absolute_measurements.tiff", width = 8, height = 10, units = 'in', res = 300)
absolute.measurements <- pheno[pheno$`English Name` %in% c('Human', 'Central chimpanzee', 'Crab-eating macaque', 'Rhesus monkey', 'Black spider monkey', 'Common squirrel monkey'), ] %>%
  group_by(`English Name`) %>% dplyr::select(`English Name`, CerebellarVol, CerebralVol, CrusVol_1_3) %>% dplyr::rename(Cerebellum = CerebellarVol, Cerebrum = CerebralVol, `Ansiform area` = CrusVol_1_3) %>%
  gather(Phenotype, Volume, Cerebellum:`Ansiform area`) %>% drop_na(Volume) %>% dplyr::filter(!Volume < 2.9) %>%
  mutate(outlier = ifelse(find_outlier(Volume), SpecimenID , NA)) %>%
  ggplot(aes(interaction(`English Name`, Phenotype), Volume)) + geom_boxplot(aes(fill = interaction(`English Name`, Phenotype)), position = "identity", outlier.colour = "red", outlier.shape = "triangle", outlier.size = 3) +
  scale_fill_manual(values = sixcol2) +
  geom_point(color="black", size=1.4, alpha=0.9) +
  theme_ipsum(base_size = 10) +
  theme(legend.position="off", legend.text = element_text(), legend.key.size = unit(2, 'cm'), # We turn the legend off and append legends from the individual measure plots.
        plot.title = element_text(size=20, face = "bold", hjust = 0.5), axis.title.y = element_text(hjust=0.5, size = 14, margin = margin(r = 40)), axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.title.x = element_text(hjust=0.5, size = 14, margin = margin(t = 30))
  ) +
  labs(y="Volume (log-transformed)", x="Species and measurement") +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=20, hjust=0.52)
absolute.measurements
graphics.off()

# Just to get a common legend for the 1.19 and 1.21 plots. A little hack-y but it works.
tiff("./1.raw_descriptives/1.20a.legend.tiff", width = 12, height = 6, units = 'in', res = 300)
box_cerebell
graphics.off()

### ---------------------------------------- ###
# Relative measurements.
### ---------------------------------------- ###  

## Cerebello-cerebral ratio
box_cc <- pheno[pheno$`English Name` %in% c('Human', 'Central chimpanzee', 'Crab-eating macaque', 'Rhesus monkey', 'Black spider monkey', 'Common squirrel monkey'), ] %>%
  group_by(`English Name`) %>%
  ggplot(aes(y= CerebroCerebellar, x=`English Name`, fill=`English Name`)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = "triangle", outlier.size = 3) +
  scale_fill_manual(values = sixcol) +
  geom_point(color="black", size=1.4, alpha=0.9) +
  theme_ipsum(base_size = 14) +
  theme(legend.position="top", legend.text = element_text(), legend.key.size = unit(3, 'cm'),
        plot.title = element_text(size=20, face = "bold", hjust = 0.5), axis.title.y = element_text(hjust=0.5, size = 20, margin = margin(r = 40)), axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle("Cerebello-cerebral ratio distribution") +
  labs(y="Cerebello-cerebral ratio (in percentages)", x="")
box_cc

## Ansiform Area Ratio with Cerebellum
box_aac <- pheno[pheno$`English Name` %in% c('Human', 'Central chimpanzee'), ] %>%
  group_by(`English Name`) %>%
  mutate(outlier = ifelse(find_outlier(CerebellarVol), SpecimenID , NA)) %>%
  ggplot(aes(y= CerebellarCrus, x=`English Name`, fill=`English Name`)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = "triangle", outlier.size = 3) +
  scale_fill_manual(values = sixcol) +
  geom_point(color="black", size=1.4, alpha=0.9) +
  theme_ipsum(base_size = 14) +
  theme(legend.position="top", legend.text = element_text(), legend.key.size = unit(3, 'cm'),
        plot.title = element_text(size=20, face = "bold", hjust = 0.5), axis.title.y = element_text(hjust=0.5, size = 20, margin = margin(r = 40)), axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle("Ansiform area ratio distribution") +
  labs(y="Ansiform area ratio (in percentages)", x="") +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier
box_aac


## Combining cerebellum/cerebrum and ansiform/cerebellar ratios.
sixcol3 <- viridis(6, alpha = 0.6, begin = 0, end = 1, direction = 1, option = "D")
sixcol3 <- c("#41448799", '#7AD15199', sixcol3)   
names(sixcol3) <- c('Central chimpanzee.AAC', 'Human.AAC', 'Black spider monkey.CC', 'Central chimpanzee.CC', 'Common squirrel monkey.CC', 'Crab-eating macaque.CC', 'Human.CC', 'Rhesus monkey.CC')

## Plot ratios together in one plot.
tiff("./1.raw_descriptives/1.21.spread_relative_measurements.tiff", width = 8, height = 10, units = 'in', res = 300)
relative.measurements <- pheno[pheno$`English Name` %in% c('Human', 'Central chimpanzee', 'Crab-eating macaque', 'Rhesus monkey', 'Black spider monkey', 'Common squirrel monkey'), ] %>%
  group_by(`English Name`) %>% dplyr::select(`English Name`, CerebroCerebellar, CerebellarCrus) %>% dplyr::rename(CC = CerebroCerebellar, AAC = CerebellarCrus) %>%
  gather(Phenotype, Ratio, CC:AAC) %>% drop_na(Ratio) %>% filter(!((`English Name` == 'Common squirrel monkey' & Phenotype == 'AAC') |
                                                                     (`English Name` == 'Crab-eating macaque' & Phenotype == 'AAC') |
                                                                     (`English Name` == 'Rhesus monkey' & Phenotype == 'AAC') )) %>%
  mutate(outlier = ifelse(find_outlier(Ratio), `English Name`, NA)) %>% 
  ggplot(aes(interaction(`English Name`, Phenotype), Ratio)) + geom_boxplot(aes(fill = interaction(`English Name`, Phenotype)), position = "identity", outlier.colour = "red", outlier.shape = "triangle", outlier.size = 3) +
  scale_fill_manual(values = sixcol3) +
  geom_point(color="black", size=1.4, alpha=0.9) +
  theme_ipsum(base_size = 10) +
  theme(legend.position="off", legend.text = element_text(), legend.key.size = unit(2, 'cm'),
        plot.title = element_text(size=20, face = "bold", hjust = 0.5), axis.title.y = element_text(hjust=0.5, size = 14, margin = margin(r = 40)), axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.title.x = element_text(hjust=0.5, size = 14, margin = margin(t = 30))
  ) +
  labs(y="Ratio (in percentages)", x="Species and ratio") +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=20, hjust=0.52)
relative.measurements
graphics.off()

### ---------------------------------------------------------- ###
# Create figure with full intraspecific variability box plots.
# NOT USED, Saved in case.
### ---------------------------------------------------------- ###

## Arrange the figures
library(gridExtra)
# Figure S1: Intraspecies variation with outliers.
S1_fig <- ggarrange(box_cerebell, 
                    box_cerebr, 
                    box_cc, 
                    arrangeGrob(box_aa, box_aac, ncol = 2),
                    labels = c("a", "b", "c", "d", "e"),
                    ncol = 2, nrow = 4, common.legend = T)
S1_fig


library(grid)
# Move to a new page
grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(box_cerebell, vp = define_region(row = 1, col = 1:2))   # Span over two columns
print(box_cerebr, vp = define_region(row = 2, col = 1:2))   # Span over two columns
print(box_cc, vp = define_region(row = 3, col = 1:2))  # Span over two columns
print(box_aa, vp = define_region(row = 4, col = 1))
print(box_aac, vp = define_region(row = 4, col = 1))



## Main text # SAVE MANUALLY; not used
figure_main1 <- ggarrange(box_cc, box_aac,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1, common.legend = T)
figure_main1

##  Supplemental # SAVE MANUALLY: not used
figure_suppl2 <- ggarrange(box_cerebell, box_cerebr, box_aa,
                          labels = c("A", "B", "C"),
                          ncol = 1, nrow = 3, common.legend = T) 
figure_suppl2



### --------------------------------------------------- ###
# Plot supplemental figure for ordered bar graphs per clade
### --------------------------------------------------- ###    

# Cerebellum
tiff("./1.raw_descriptives/1.22.barplot_cerebellum.tiff", width = 20, height = 5, units = 'in', res = 300)
pheno_split %>% 
  ggplot(aes(fct_reorder(SpecimenID,
                         CerebellarVol), 
             CerebellarVol))+
  geom_col(aes(col = Clade, fill = Clade)) +   theme_classic() +
  labs(x="Specimen ID", y = "Cerebellar volume") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_text(size = 0, margin = margin(t = 0)), axis.title.y = element_text(size = 14, margin = margin(r = 40))) 
graphics.off()

# Cerebrum
tiff("./1.raw_descriptives/1.23.barplot_cerebrum.tiff", width = 20, height = 5, units = 'in', res = 300)
pheno_split %>% 
  ggplot(aes(fct_reorder(SpecimenID,
                         CerebralVol), 
             CerebralVol))+
  geom_col(aes(col = Clade, fill = Clade)) +   theme_classic() +
  labs(x="Specimen ID", y = "Cerebral volume") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_text(size = 0, margin = margin(t = 0)), axis.title.y = element_text(size = 14, margin = margin(r = 40))) 
graphics.off()

# Ansiform Area
tiff("./1.raw_descriptives/1.24.barplot_ansiform.tiff", width =10, height = 5, units = 'in', res = 300)
pheno_split %>% drop_na(CrusVol_1_3) %>%
  ggplot(aes(fct_reorder(SpecimenID,
                         CrusVol_1_3), 
             CrusVol_1_3))+
  geom_col(aes(col = Clade, fill = Clade)) +   theme_classic() +
  labs(x="Specimen ID", y = "Ansiform area volume") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_text(size = 0, margin = margin(t = 0)), axis.title.y = element_text(size = 14, margin = margin(r = 40))) 
graphics.off()

# Ratio Cerebro/Cerebellar
tiff("./1.raw_descriptives/1.25.barplot_cerebroCerebellar.tiff", width = 20, height = 5, units = 'in', res = 300)
pheno_split %>% 
  ggplot(aes(fct_reorder(SpecimenID,
                         CerebroCerebellar), 
             CerebroCerebellar))+
  geom_col(aes(col = Clade, fill = Clade)) +   theme_classic() +
  labs(x="Specimen ID", y = "Ratio cerebellum/cerebrum") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_text(size = 0, margin = margin(t = 0)), axis.title.y = element_text(size = 14, margin = margin(r = 40))) 
graphics.off()

# Ratio Ansiform area/Cerebellum
tiff("./1.raw_descriptives/1.26.barplot_ansiformRatio.tiff", width =10, height = 5, units = 'in', res = 300)
pheno_split %>% drop_na(CrusVol_1_3) %>%
  ggplot(aes(fct_reorder(SpecimenID,
                         CerebellarCrus), 
             CerebellarCrus))+
  geom_col(aes(col = Clade, fill = Clade)) +   theme_classic() +
  labs(x="Specimen ID", y = "Ratio ansiform area/cerebellum") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_text(size = 0, margin = margin(t = 0)), axis.title.y = element_text(size = 14, margin = margin(r = 40))) 
graphics.off()


### --------------------------------------------------- ###
# Box plots per clade to visualize distribution.
# Not really used at the moment. Kept in case.
# Also, due to the obvious effect of phylogeny, as well as 
# body size, box plots of species median value have little
# value other than to describe data spread.
### --------------------------------------------------- ### 


# Cerebellum
pdf('1.raw_descriptives/1.27.boxplot_cerebellum.pdf')
bcerebellum <- phenoMedian %>%
  group_by(Clade) %>%
  mutate(outlier = ifelse(find_outlier(MedianCerebellum), WilsonReederName , NA)) %>%
  ggplot(aes(x = Clade, y = MedianCerebellum, fill = Clade)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Clade")) + 
  labs(title="Cerebellar Volume",x="Clade", y = "Cerebellar Volume") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier
bcerebellum
dev.off()

# Cerebrum
pdf('1.raw_descriptives/1.28.boxplot_cerebrum.pdf')
bcerebrum <- phenoMedian %>% 
  group_by(Clade) %>%
  mutate(outlier = ifelse(find_outlier(MedianCerebrum), WilsonReederName , NA)) %>%
  ggplot(aes(x = Clade, y = MedianCerebrum, fill = Clade)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Clade")) +
  labs(title="Cerebral Volume",x="Clade", y = "Cerebral Volume") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier
bcerebrum
dev.off()

# Ansiform Area
pdf('1.raw_descriptives/1.29.boxplot_ansiformarea.pdf')
bAA <- phenoMedian %>%
  ggplot(aes(x = Clade, y = MedianCrus, fill = Clade)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Clade")) +
  labs(title="Ansiform Area Volume",x="Clade", y = "Ansiform Area Volume") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1)) 
bAA
dev.off()

# Cerebellar/Cerebral Volume
pdf('1.raw_descriptives/1.30.boxplot_cerebrocerebellar.pdf')
bcc <- phenoMedian %>%
  group_by(Clade) %>%
  mutate(outlier = ifelse(find_outlier(MedianCerebroCerebellar), WilsonReederName , NA)) %>%
  ggplot(aes(x = Clade, y = MedianCerebroCerebellar, fill = Clade)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Clade")) + 
  labs(title="Cerebello/Cerebral Volume Ratio", x="Clade", y = "Cerebellar/ Cerebral Volume") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier
bcc
dev.off()

# Crus/Cerebellar Volume
pdf('1.raw_descriptives/1.31.boxplot_aacerebellar.pdf')
bAAc <- bind_rows(phenoMedian %>% filter(Clade == c("Hominoidea")), phenoMedian %>% filter(Clade == c("Cebidae")), phenoMedian %>% filter(Clade == c("Papionini"))) %>%
  group_by(Clade) %>% drop_na(MedianCerebellarCrus) %>%
  mutate(outlier = ifelse(find_outlier(MedianCerebellarCrus), WilsonReederName , NA)) %>%
  ggplot(aes(x = Clade, y = MedianCerebellarCrus, fill = Clade)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Clade")) + 
  labs(title="Ansiform Area/ Whole Cerebellar Volume Ratio",x="Clade", y = "Ansiform Area/ Cerebellar Volume") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier 
bAAc
dev.off()

### ------------------------------------- ###
# Create new phenotype data, exclude outliers.
### ------------------------------------- ###  

## We will perform regressions with data excluding outliers!

## Hence, we also make a new median data table! ##
phenoOutliers_raw <- pheno_raw %>% filter(!SpecimenID %in% c("MacaqueCrabierMa_8f9e", "MacaqueRhesusMac_bdf8"))
write_csv2(phenoOutliers_raw, "./0.cleanInput/Phenotypes_raw_outliersRM.csv")

## Log-transformed ##
phenoOutliers <- pheno %>% filter(!SpecimenID %in% c("MacaqueCrabierMa_8f9e", "MacaqueRhesusMac_bdf8"))
write_csv2(phenoOutliers, "./0.cleanInput/Phenotypes_outliersRM.csv")

## Median, raw ##
phenoMedianOutliers_raw <- phenoOutliers_raw %>%
  group_by(WilsonReederName) %>%
  dplyr::summarize(MedianCerebellum = median(CerebellarVol, na.rm=TRUE)) %>% left_join(phenoOutliers_raw %>%
                                                                                         group_by(WilsonReederName) %>%
                                                                                         dplyr::summarize(MedianCerebrum = median(CerebralVol, na.rm=TRUE)), by = "WilsonReederName") %>% left_join(phenoOutliers_raw %>%
                                                                                                                                                                                                      group_by(WilsonReederName) %>%
                                                                                                                                                                                                      dplyr::summarize(MedianCrus = median(CrusVol_1_3, na.rm=TRUE)), by = "WilsonReederName") %>% left_join(phenoOutliers_raw %>%
                                                                                                                                                                                                                                                                                                               group_by(WilsonReederName) %>%
                                                                                                                                                                                                                                                                                                               dplyr::summarize(MedianROC = median(ROC, na.rm=TRUE)), by = "WilsonReederName")

## Append species information
phenoMedianOutliers_raw <- phenoMedianOutliers_raw %>%
  left_join(lut, by="WilsonReederName") %>% dplyr::select(c("WilsonReederName","Hominoidea", "Clade", "MedianCerebellum", "MedianCerebrum", "MedianCrus", "MedianROC")) %>% distinct()

## We also add columns for the ratios of Cerebellar/Cerebral volume, and Crus/Cerebellar Volume, using only species medians.
# We have to first 'unlog' the values, so we can get percentages for the ratios.
phenoMedianOutliers_raw <- phenoMedianOutliers_raw %>%
  mutate(MedianCerebroCerebellar = ((MedianCerebellum)/(MedianCerebrum))*100) %>% mutate(MedianCerebellarCrus = ((MedianCrus)/(MedianCerebellum))*100) %>% mutate(MedianROCCrus = ((MedianCrus)/(MedianROC))*100)

write_csv2(phenoMedianOutliers_raw, "./0.cleanInput/PhenotypesMedian_raw_outliersRM.csv")

## Log-transform ##
phenoMedianOutliers <- phenoMedianOutliers_raw
phenoMedianOutliers$MedianCerebellum <- log10(phenoMedianOutliers$MedianCerebellum)
phenoMedianOutliers$MedianCerebrum <- log10(phenoMedianOutliers$MedianCerebrum)
phenoMedianOutliers$MedianCrus <- log10(phenoMedianOutliers$MedianCrus)
phenoMedianOutliers$MedianROC <- log10(phenoMedianOutliers$MedianROC)
# Cerebellum/Cerebral and Ansiform Area/Cerebellum ratio do not need to be log-transformed, as we turned them into percentages earlier (in same order of magnitude).

write_csv2(phenoMedianOutliers, "./0.cleanInput/PhenotypesMedian_outliersRM.csv")


### ------------------------------------------- ###
# Correlate variables not corrected for phylogeny 
# ('raw': means uncorrected in a phylogenetic sense; 
# Includes raw data AND log-transformed data).  
### ------------------------------------------- ###  

## Fit linear models using generalized least squares

## Load libraries
library(nlme) 
#library(lmtest)
library(ggpmisc)
library(psych)

### ---------------------------------------- ###
# Simple linear mode, check heteroscedasticity.
### ---------------------------------------- ###    

f.full.lm <- lm(CerebellarVol ~ CerebralVol, data=pheno) # initial model
f.apes.lm <- lm(CerebellarVol ~ CerebralVol, data=pheno, subset = Hominoidea == "Ape")
f.nonapes.lm <- lm(CerebellarVol ~ CerebralVol, data=pheno,subset = Hominoidea == "Nonape")
## Visualize
par(mfrow=c(2,2)) # init 4 charts in 1 panel
plot(f.full.lm)

par(mfrow=c(2,2))
plot(f.apes.lm)

par(mfrow=c(2,2))
plot(f.nonapes.lm)
## Test statistically with the Breusch-Pagan test
lmtest::bptest(f.full.lm)  # Not sign.
lmtest::bptest(f.apes.lm) # Sign.: Even more strongly influenced by several observations per species for Chimpanzees and Humans, who both have high values.
lmtest::bptest(f.nonapes.lm) # Not sign.

### ------------------------------------------------------------------------- ###
# Ordinary least-squares regression
# For species Median, excluding outlier Crab-eating macaque and Rhesus monkey.
# USE THIS, or lm in 5.phylogenetic regression?
# They should amount to the same (regressions using median values, that is).
# We can keep this as an alternative visualization.
### ------------------------------------------------------------------------- ###

## Cerebellum vs Cerebrum
## All species
RegressionFull_CC <- ggplot(data = phenoMedianOutliers, aes(x = MedianCerebrum, y = MedianCerebellum)) + 
  geom_point(aes(col=WilsonReederName)) +
  theme_classic() +
  ggtitle("Full data") +
  ylab("") +
  xlab("") +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))
RegressionFull_CC

## Apes and non-apes separately
Regression_Apes_CC <- ggplot(phenoMedianOutliers, aes(x = MedianCerebrum, y = MedianCerebellum, linetype = Hominoidea)) +
  geom_point(aes(col=WilsonReederName))  +
  theme_classic() +
  ggtitle("Apes versus nonapes") +
  ylab("") +
  xlab("") +
  geom_smooth(method = "lm", fill = NA) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  theme_classic() +
  scale_linetype_manual(values=c("twodash", "solid"))
Regression_Apes_CC

### -------------------------------------------------------------- ###
# Create OLS regression plots, supplemental figure, species medians.
### -------------------------------------------------------------- ###    

figure_suppl4 <- ggarrange(RegressionFull_CC, Regression_Apes_CC,
                          labels = c("A", "B"),
                          ncol = 2, nrow = 1, common.legend = T, legend = "bottom") +
  ggtitle("Cerebellar versus cerebral volume") # title does not show yet
figure_suppl4 <- annotate_figure(figure_suppl4, left = textGrob("Cerebellar volume (in mm3, log-transformed)", rot = 90, vjust = 0.5, hjust = 0.25, gp = gpar(cex = 1.3)))
figure_suppl4 <- annotate_figure(figure_suppl4, bottom = textGrob("Cerebral volume (in mm3, log-transformed)", vjust = -13, hjust = 0.5, gp = gpar(cex = 1.3)))
figure_suppl4 # SAVE MANUALLY

## As expected, the ape regression sees the largest effect of summarizing the data into median values, due to the higher number of Chimpanzee and Human observations.
f.full.lm.med <- lm(MedianCerebellum ~ MedianCerebrum, data=phenoMedianOutliers) # initial model
f.apes.lm.med <- lm(MedianCerebellum ~ MedianCerebrum, data=phenoMedianOutliers, subset = Hominoidea == "Ape")
f.nonapes.lm.med <- lm(MedianCerebellum ~ MedianCerebrum, data=phenoMedianOutliers,subset = Hominoidea == "Nonape")
## Test statistically with the Breusch-Pagan test
lmtest::bptest(f.full.lm.med)  # Not sign.
lmtest::bptest(f.apes.lm.med) # Not sign.
lmtest::bptest(f.nonapes.lm.med) # Not sign.

## Ansiform Area vs. Rest of Cerebellar volume
## All species
RegressionFull_AAC <- ggplot(data = phenoMedianOutliers, aes(x = MedianROC, y = MedianCrus)) + 
  geom_point(aes(col=WilsonReederName)) +
  theme_classic() +
  ggtitle("Full data") +
  ylab("") +
  xlab("") +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))
RegressionFull_AAC

## Apes and non-apes separately
Regression_Apes_AAC <- ggplot(phenoMedianOutliers, aes(x = MedianROC, y = MedianCrus, linetype = Hominoidea)) +
  geom_point(aes(col=WilsonReederName))  +
  theme_classic() +
  ggtitle("Apes versus nonapes") +
  ylab("") +
  xlab("") +
  geom_smooth(method = "lm", fill = NA) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  theme_classic() +
  scale_linetype_manual(values=c("twodash", "solid"))
Regression_Apes_AAC

figure_suppl5 <- ggarrange(RegressionFull_AAC, Regression_Apes_AAC,
                           labels = c("A", "B"),
                           ncol = 2, nrow = 1, common.legend = T, legend = "bottom") + 
  ggtitle("Ansiform area versus cerebellar volume") # title does not show yet
figure_suppl5 <- annotate_figure(figure_suppl5, left = textGrob("Ansiform area volume (in mm3, log-transformed)", rot = 90, vjust = 0.5, hjust = 0.25, gp = gpar(cex = 1.3)))
figure_suppl5 <- annotate_figure(figure_suppl5, bottom = textGrob(" Rest of Cerebellar volume (in mm3, log-transformed)", vjust = -13, hjust = 0.5, gp = gpar(cex = 1.3)))
figure_suppl5 # SAVE MANUALLY

## Ansiform Area vs. Cerebral Volume
## All species
RegressionFull_AACerebr <- ggplot(data = phenoMedianOutliers, aes(x = MedianCerebrum, y = MedianCrus)) + 
  geom_point(aes(col=WilsonReederName)) +
  theme_classic() +
  ggtitle("Full data") +
  ylab("") +
  xlab("") +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))
RegressionFull_AACerebr

## Apes and non-apes separately
Regression_Apes_AACerebr <- ggplot(phenoMedianOutliers, aes(x = MedianCerebrum, y = MedianCrus, linetype = Hominoidea)) +
  geom_point(aes(col=WilsonReederName))  +
  theme_classic() +
  ggtitle("Apes versus nonapes") +
  ylab("") +
  xlab("") +
  geom_smooth(method = "lm", fill = NA) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  theme_classic() +
  scale_linetype_manual(values=c("twodash", "solid"))
Regression_Apes_AACerebr

figure_suppl6 <- ggarrange(RegressionFull_AACerebr, Regression_Apes_AACerebr,
                           labels = c("A", "B"),
                           ncol = 2, nrow = 1, common.legend = T, legend = "bottom") + 
  ggtitle("Ansiform area versus cerebral volume")
figure_suppl6 <- annotate_figure(figure_suppl6, left = textGrob("Ansiform area volume (in mm3, log-transformed)", rot = 90, vjust = 0.5, hjust = 0.25, gp = gpar(cex = 1.3)))
figure_suppl6 <- annotate_figure(figure_suppl6, bottom = textGrob("Cerebral volume (in mm3, log-transformed)", vjust = -13, hjust = 0.5, gp = gpar(cex = 1.3)))
figure_suppl6 # SAVE MANUALLY

### --------------------------------------------------------------- ###
# R-values of the 3 regressions with Fisher's R-to-Z transformation.
### --------------------------------------------------------------- ###    

# Cerebellum vs Cerebrum
f.full.lm.med <- lm(MedianCerebellum ~ MedianCerebrum, data=phenoMedianOutliers)
f.apes.lm.med <- lm(MedianCerebellum ~ MedianCerebrum, data=phenoMedianOutliers, subset = Hominoidea == "Ape")
f.nonapes.lm.med <- lm(MedianCerebellum ~ MedianCerebrum, data=phenoMedianOutliers,subset = Hominoidea == "Nonape")
paired.r(summary(f.full.lm.med)$r.squared, summary(f.apes.lm.med)$r.squared, summary(f.nonapes.lm.med)$r.squared, n = 34)

# Ansiform area vs Rest of Cerebellum
f.full.lm.med.aaROC <- lm(MedianCrus ~ MedianROC, data=phenoMedianOutliers)
f.apes.lm.med.aaROC <- lm(MedianCrus ~ MedianROC, data=phenoMedianOutliers, subset = Hominoidea == "Ape")
f.nonapes.lm.med.aaROC <- lm(MedianCrus ~ MedianROC, data=phenoMedianOutliers,subset = Hominoidea == "Nonape")
paired.r(summary(f.full.lm.med.aaROC)$r.squared, summary(f.apes.lm.med.aaROC)$r.squared, summary(f.nonapes.lm.med.aaROC)$r.squared, n = 13)

# Ansiform area vs Cerebrum
f.full.lm.med.aacerebr <- lm(MedianCrus ~ MedianCerebrum, data=phenoMedianOutliers) 
f.apes.lm.med.aacerebr <- lm(MedianCrus ~ MedianCerebrum, data=phenoMedianOutliers, subset = Hominoidea == "Ape")
f.nonapes.lm.med.aacerebr <- lm(MedianCrus ~ MedianCerebrum, data=phenoMedianOutliers,subset = Hominoidea == "Nonape")
paired.r(summary(f.full.lm.med.aacerebr)$r.squared, summary(f.apes.lm.med.aacerebr)$r.squared, summary(f.nonapes.lm.med.aacerebr)$r.squared, n = 13)


### -------------------------- ###
### END of script
### -------------------------- ###   


