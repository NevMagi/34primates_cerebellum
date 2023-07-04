############################################
### Distribution of scaling in brain volumes
### and all phenotypes of interest (incl. ratios)
### per sex (and also age for humans). Although
### these analyses are definitely liited by 
### the small sample size, they may show whether
### in humans and chimpanzees, ratios are sensitive
### to sex or age. 
############################################

library(tidyverse)
library(magrittr)
library(ggplot2)
library(rstatix)
library(corrplot)


### SET UP ###
getwd()
setwd("~/Documents/PhD/EvolutionPrimates/RevisionCode")

outputdir <- "7.Sex_age"
dir.create(outputdir)

sex <- read_csv2("./input/Phenotypes_split_BrainVolume.csv") %>% dplyr::rename(species = WilsonReederName) %>% dplyr::select(species, CerebellarVol, CerebralVol, CrusVol_1_3, ROC, CerebroCerebellar, CerebellarCrus, ROCCrus, BrainVol, Sex, Age) %>% drop_na(Sex) %>% as.data.frame()
as.factor(sex$Sex)
#sex2 <- sex %>% select(species, Sex, CerebellarVol, CerebralVol, CrusVol_1_3, BrainVol)
sex$CerebellarVol <- 10^(sex$CerebellarVol) # Unlog
sex$CerebralVol <- 10^(sex$CerebralVol) # Unlog
sex$CrusVol_1_3 <- 10^(sex$CrusVol_1_3) # Unlog
sex$BrainVol <- 10^(sex$BrainVol) # Unlog
sex$CerebellarNormalized <- sex$CerebellarVol/sex$BrainVol # Normalize
sex$CerebralNormalized <- sex$CerebralVol/sex$BrainVol # Normalize
sex$CrusNormalized <- sex$CrusVol_1_3/sex$BrainVol # Normalize
sex$CerebellarVol <- log10(sex$CerebellarVol) # Reassign
sex$CerebralVol <- log10(sex$CerebralVol) # Reassign
sex$CrusVol_1_3 <- log10(sex$CrusVol_1_3) # Reassign
sex$BrainVol <- log10(sex$BrainVol) # Reassign

#sex <- sex[c(1,10,11,9,5,2,3,4,7,6,8,12,13,14)]

## Long formatting
homosapiens <- sex %>% filter(species== c("Homo_sapiens"))
homosapiens <- homosapiens %>% select(-species)
homosapiens$Sex <- as.factor(homosapiens$Sex)
homosapiens <- homosapiens %>%
  as_tibble()
homosapiens.long <- homosapiens %>%
  pivot_longer(-Sex, names_to = "variables", values_to = "value")
homosapiens.long$variables <- as.factor(homosapiens.long$variables)

## Grouped t-test
stat.test.hs <- homosapiens.long %>%
  group_by(variables) %>%
  t_test(value ~ Sex) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.hs

# Chimpanzee
pantroglodytes <- sex %>% filter(species== c("Pan_troglodytes_troglodytes"))
pantroglodytes <- pantroglodytes %>% select(-species, -Age)
pantroglodytes$Sex <- as.factor(pantroglodytes$Sex)
pantroglodytes <- pantroglodytes %>%
  as_tibble()
pantroglodytes.long <- pantroglodytes %>%
  pivot_longer(-Sex, names_to = "variables", values_to = "value")
pantroglodytes.long$variables <- as.factor(pantroglodytes.long$variables)

## Grouped t-test
stat.test.ptt <- pantroglodytes.long %>%
  group_by(variables) %>%
  t_test(value ~ Sex) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.ptt

## To be able to label the outliers in individual trait plots, we create a function to identify them ##
find_outlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}


#####################################
### Distribution sex

### + T-tests between the sexes. From inspection
### of the data it already becomes clear that
### population distribution is significantly
### larger than the sample size is able to capture.
### 
#####################################

## Cerebellar Volume
pdf('7.Sex_age/7.1.boxplot_cerebellar-hs.pdf')
hs_cerebelllarvol <- homosapiens %>%
  group_by(Sex) %>%
  mutate(outlier = ifelse(find_outlier(CerebellarVol), CerebellarVol , NA)) %>%
  ggplot(aes(x = Sex, y = CerebellarVol, fill = Sex)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Sex")) + 
  labs(title="Cerebellar Volume Humans",x="Sex", y = "Cerebellar Volume") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier 
hs_cerebelllarvol
dev.off()

pdf('7.Sex_age/7.1a.boxplot_cerebellar-ptt.pdf')
ptt_cerebelllarvol <- pantroglodytes %>%
  group_by(Sex) %>%
  mutate(outlier = ifelse(find_outlier(CerebellarVol), CerebellarVol , NA)) %>%
  ggplot(aes(x = Sex, y = CerebellarVol, fill = Sex)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Sex")) + 
  labs(title="Cerebellar Volume Chimpanzees",x="Sex", y = "Cerebellar Volume") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier 
ptt_cerebelllarvol
dev.off()

## Cerebral Volume
pdf('7.Sex_age/7.2.boxplot_cerebral-hs.pdf')
hs_cerebral <- homosapiens %>%
  group_by(Sex) %>%
  mutate(outlier = ifelse(find_outlier(CerebralVol), CerebralVol , NA)) %>%
  ggplot(aes(x = Sex, y = CerebralVol, fill = Sex)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Sex")) + 
  labs(title="Cerebral Volume Humans",x="Sex", y = "Cerebral Volume") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier 
hs_cerebral
dev.off()

pdf('7.Sex_age/7.2a.boxplot_cerebral-ptt.pdf')
ptt_cerebral <- pantroglodytes %>%
  group_by(Sex) %>%
  mutate(outlier = ifelse(find_outlier(CerebralVol), CerebralVol , NA)) %>%
  ggplot(aes(x = Sex, y = CerebralVol, fill = Sex)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Sex")) + 
  labs(title="Cerebral Volume Chimpanzees",x="Sex", y = "Cerebral Volume") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier 
ptt_cerebral
dev.off()


## Brain Volume
pdf('7.Sex_age/7.3.boxplot_brain-hs.pdf')
hs_brain <- homosapiens %>%
  group_by(Sex) %>%
  mutate(outlier = ifelse(find_outlier(BrainVol), BrainVol , NA)) %>%
  ggplot(aes(x = Sex, y = BrainVol, fill = Sex)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Sex")) + 
  labs(title="Brain Volume Humans",x="Sex", y = "Brain Volume") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier 
hs_brain
dev.off()

pdf('7.Sex_age/7.3a.boxplot_brain-ptt.pdf')
ptt_brain <- pantroglodytes %>%
  group_by(Sex) %>%
  mutate(outlier = ifelse(find_outlier(BrainVol), BrainVol , NA)) %>%
  ggplot(aes(x = Sex, y = BrainVol, fill = Sex)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Sex")) + 
  labs(title="Brain Volume Chimpanzees",x="Sex", y = "Brain Volume") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier 
ptt_brain
dev.off()

## Cerebro-cerebellar ratio
pdf('7.Sex_age/7.4.boxplot_cerebrocerebellar-hs.pdf')
hs_cc <- homosapiens %>%
  group_by(Sex) %>%
  mutate(outlier = ifelse(find_outlier(CerebroCerebellar),CerebroCerebellar , NA)) %>%
  ggplot(aes(x = Sex, y = CerebroCerebellar, fill = Sex)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Sex")) + 
  labs(title="Ratio Cerebellum/Cerebrum Humans",x="Sex", y = "Ratio Cerebellum/Cerebrum") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier 
hs_cc
dev.off()

pdf('7.Sex_age/7.4a.boxplot_cerebrocerebellar-ptt.pdf')
ptt_cc <- pantroglodytes %>%
  group_by(Sex) %>%
  mutate(outlier = ifelse(find_outlier(CerebroCerebellar),CerebroCerebellar , NA)) %>%
  ggplot(aes(x = Sex, y = CerebroCerebellar, fill = Sex)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Sex")) + 
  labs(title="Ratio Cerebellum/Cerebrum Chimpanzees",x="Sex", y = "Ratio Cerebellum/Cerebrum") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier 
ptt_cc
dev.off()

## Ansiform area-cerebellar ratio
pdf('7.Sex_age/7.5.boxplot_cerebellarcrus-hs.pdf')
hs_aac <- homosapiens %>%
  group_by(Sex) %>%
  mutate(outlier = ifelse(find_outlier(CerebellarCrus),CerebellarCrus , NA)) %>%
  ggplot(aes(x = Sex, y = CerebellarCrus, fill = Sex)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Sex")) + 
  labs(title="Ratio Ansiform Area/Cerebellum Humans",x="Sex", y = "Ratio Ansiform Area/Cerebellum") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier 
hs_aac
dev.off()

pdf('7.Sex_age/7.5a.boxplot_cerebellarcrus-ptt.pdf')
ptt_aac <- pantroglodytes %>%
  group_by(Sex) %>%
  mutate(outlier = ifelse(find_outlier(CerebellarCrus),CerebellarCrus , NA)) %>%
  ggplot(aes(x = Sex, y = CerebellarCrus, fill = Sex)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(outlier.shape = "triangle", outlier.size = 2, outlier.color = "red") +
  guides(fill = guide_legend(title = "Sex")) + 
  labs(title="Ratio Ansiform Area/Cerebellum Chimpanzees",x="Sex", y = "Ansiform Area/Cerebrum") + theme_classic() +
  theme(axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(label=outlier), na.rm=TRUE, vjust=2) ## Mark the outlier 
ptt_aac
dev.off()


## HUMAN PLOTS

# T-tests (likely quite questionable due to the sample size; but shows the difference between relative and absolute traits)
# Volumes
cerebell.t <- t.test(CerebellarVol ~ Sex, data=homosapiens) # Sign.
cerebral.t <- t.test(CerebralVol ~ Sex, data=homosapiens)  # Sign.
crus.t <- t.test(CrusVol_1_3 ~ Sex, data=homosapiens) # Sign.
ROC.t <- t.test(ROC ~ Sex, data=homosapiens) # Sign.
BrainVol.t <- t.test(BrainVol ~ Sex, data=homosapiens) # Sign.
# Ratios
CC.t.test <- t.test(CerebroCerebellar ~ Sex, data=homosapiens) # Non-sign
AAC.t.test <- t.test(CerebellarCrus ~ Sex, data=homosapiens) # Non-sign
AAROC.t.test <- t.test(ROCCrus ~ Sex, data=homosapiens) # Non-sign.

## Facet-wrapped plot
plot.hs <- ggboxplot(
  homosapiens.long, x = "Sex", y = "value",
  fill = "Sex", palette = "RdPu", legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) +
  facet_wrap(~factor(variables, levels =c('Age', 'CerebellarVol',  'CerebellarCrus','CerebellarNormalized', 'BrainVol','CerebralVol','CerebroCerebellar', 'CerebralNormalized','CrusVol_1_3','ROC', 'ROCCrus','CrusNormalized')), scales = "free_y", )
plot.hs

# Add statistical test p-values
stat.test.hs <- stat.test.hs %>% add_xy_position(x = "Sex", scales = c("free_y"))

pdf('7.Sex_age/7.6.boxplots_humans_significance.pdf')
plot.hs + stat_pvalue_manual(stat.test.hs, label = "p.adj.signif", color = "red", size = 2.5)
graphics.off()


## Rename variables to be pretty for publication.
variable_labels <- c(
  Age = "Age in years",
  CerebellarVol = "Cerebellar volume",
  CerebellarCrus = "Ansiform area-to-cerebellum ratio",
  CerebellarNormalized = "Normalized cerebellar volume",
  BrainVol = "Brain volume",
  CerebralVol = "Cerebral volume",
  CerebroCerebellar = "Cerebellum-to-cerebrum ratio",
  CerebralNormalized = "Normalized cerebral volume",
  CrusVol_1_3 = "Ansiform area volume",
  ROC = "Rest of cerebellum volume",
  ROCCrus = "Ansiform area-to-rest of cerebellar volume ratio",
  CrusNormalized = "Normalized ansiform area volume"
)


## Redo as a violin plot.
pdf('7.Sex_age/7.6a.violins_humans_significance.pdf', width = 16, height=8)
plot.hs.violin <- ggviolin(
  homosapiens.long, x = "Sex", y = "value",
  fill = "Sex", palette = "RdPu", legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) +
  facet_wrap(~factor(variables, levels =c('Age', 'CerebellarVol',  'CerebellarCrus','CerebellarNormalized', 'BrainVol','CerebralVol','CerebroCerebellar', 'CerebralNormalized','CrusVol_1_3','ROC', 'ROCCrus','CrusNormalized')), scales = "free_y", labeller = as_labeller(variable_labels))
plot.hs.violin + 
  geom_boxplot(width=0.1) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  stat_pvalue_manual(stat.test.hs, label = "p.adj.signif", color = "red") +
  ylab("") + xlab("")
graphics.off()


## CHIMPANZEE PLOTS
# T-tests (likely quite questionable due to the sample size; but shows the difference between relative and absolute traits)
# Volumes
cerebell.t.ptt <- t.test(CerebellarVol ~ Sex, data=pantroglodytes) # Sign.
cerebral.t.ptt <- t.test(CerebralVol ~ Sex, data=pantroglodytes)  # Sign.
crus.t.ptt <- t.test(CrusVol_1_3 ~ Sex, data=pantroglodytes) # Sign.
ROC.t.ptt <- t.test(ROC ~ Sex, data=pantroglodytes) # Sign.
BrainVol.t.ptt <- t.test(BrainVol ~ Sex, data=pantroglodytes) # Sign.
# Ratios
CC.t.test.ptt <- t.test(CerebroCerebellar ~ Sex, data=pantroglodytes) # Non-sign
AAC.t.test.pt <- t.test(CerebellarCrus ~ Sex, data=pantroglodytes) # Non-sign
AAROC.t.test.ptt <- t.test(ROCCrus ~ Sex, data=pantroglodytes) # Non-sign.


## Add an artifical age plot, which is empty because data does not exist.
## This will make the plot more visually comparable to the human plot
pantroglodytes2 <- sex %>% filter(species== c("Pan_troglodytes_troglodytes"))
pantroglodytes2 <- pantroglodytes2 %>% select(-species)
pantroglodytes2$Sex <- as.factor(pantroglodytes$Sex)
pantroglodytes2[is.na(pantroglodytes2)] <-  runif(sum(is.na(pantroglodytes2$Age)), min = 18, max = 25)
pantroglodytes2 <- pantroglodytes2 %>%
  as_tibble()
pantroglodytes2.long <- pantroglodytes2 %>%
  pivot_longer(-Sex, names_to = "variables", values_to = "value")


# Add statistical test p-values
stat.test.ptt <- stat.test.ptt %>% add_xy_position(x = "Sex", scales = c("free_y"))


## Facet-wrapped plot
plot.ptt <- ggboxplot(
  pantroglodytes.long, x = "Sex", y = "value",
  fill = "Sex", palette = "Purples", legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) +
  facet_wrap(~factor(variables, levels =c('Age', 'CerebellarVol',  'CerebellarCrus','CerebellarNormalized', 'BrainVol','CerebralVol','CerebroCerebellar', 'CerebralNormalized','CrusVol_1_3','ROC', 'ROCCrus','CrusNormalized')), scales = "free_y")
plot.ptt
# Add statistical test p-values
stat.test.ptt <- stat.test.ptt %>% add_xy_position(x = "Sex", scales = c("free_y"))
pdf('7.Sex_age/7.7.boxplots_chimpanzees_significance.pdf')
plot.ptt + stat_pvalue_manual(stat.test.ptt, label = "p.adj.signif", color = "red", size = 2.5) 
graphics.off()


## Redo as a violin plot.
pdf('7.Sex_age/7.7a.violins_chimpanzees_significance.pdf', width = 16, height=8)
plot.ptt.violin <- ggviolin(
  pantroglodytes2.long, x = "Sex", y = "value",
  fill = "Sex", palette = "Purples", legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) +
  facet_wrap(~factor(variables, levels =c('Age', 'CerebellarVol',  'CerebellarCrus','CerebellarNormalized', 'BrainVol','CerebralVol','CerebroCerebellar', 'CerebralNormalized','CrusVol_1_3','ROC', 'ROCCrus','CrusNormalized')), scales = "free_y", labeller = as_labeller(variable_labels))
plot.ptt.violin + 
  geom_boxplot(width=0.1) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  stat_pvalue_manual(stat.test.ptt, label = "p.adj.signif", color = "red") +
  ylab("") + xlab("")
graphics.off()




#####################################
### AGE
### Just for humans, check if there is a 
### systematic relationship between trait 
### values and age.
#####################################
age <- sex %>% drop_na(Age) %>% select (-species, -Sex)
age <- age %>% 
  as.matrix %>%
  cor 


# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat) 
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
## matrix of the p-value of the correlation
p.mat <- cor.mtest(age)
head(p.mat[, 1:9])


## To display both the correlation and significance level, we have to change corrplot slightly. Otherwise both overlap.
## https://stackoverflow.com/questions/63227830/r-corrplot-plot-correlation-coefficients-along-with-significance-stars
## This have to be done manually upon loading of the package, so just running the whole script will make the figure look less than ideal.
corrplot(age, type="upper", tl.col = "black", tl.cex = 0.6, tl.srt = 45, method = "number",
         p.mat = p.mat, sig.level = c(.001, .01, .05), insig = "label_sig", diag = F, order = "original",
         col = colorRampPalette(c("darkblue", "white", "darkred"))(100)) ### (uncorrected)

pdf('7.Sex_age/7.8.humans_correlation_with_age.pdf')
corrplot(age, type="upper", tl.col = "black", tl.cex = 0.6, tl.srt = 45, method = "number",
         p.mat = p.mat,  sig.level = c(.001/36, .01/36, .05/36), insig = "label_sig",, diag = F, order = "original",
         col = colorRampPalette(c("darkblue", "white", "darkred"))(100)) ## Maximally conservative (/n)
graphics.off()

#################################
#### END OF SCRIPT
#################################

