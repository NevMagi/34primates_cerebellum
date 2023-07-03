## Set working directory
getwd()
setwd("~/Documents/PhD/EvolutionPrimates/Final_code/")

## Data input
AAmass <- read_csv2("./input/AA_corrMASS.csv")
AAmass$lograt <- log10(AAmass$ratio)
AAmass$logmass <- log10(AAmass$bodymass)

CCmass <- read_csv2("./0.cleanInput/PhenotypesMedian.csv") %>% 
  dplyr::rename(mass = `Body mass species mean`) %>%
  filter(!is.na(mass))
CCmass$CCrat <- log10(CCmass$MedianCerebroCerebellar) 

### -------------------------------------------------------- ###
### Corrected for phylogeny.
### -------------------------------------------------------- ### 

# Mod: Ansiform area: cerebellum ratio regressed on body mass
tiff("./FOLDER/5.23.BM-PGLS-AAratio-bodymass-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod5 <- gls(lograt ~ logmass, data=AAmass, correlation=bm13)
plot(mod5)
plot(AAmass$lograt ~ AAmass$logmass, xlab = "Body mass (log-transformed)", ylab = "Ratio ansiform area/ cerebellum (log-transformed)")
abline(mod5)
abline(a=1.165705, b=1.0, col= "blue") #MANUAL INPUT!! Is this correct?? Don't think so. B is incorrect (negative logged values).
# Confidence intervals
pGLS.mod5.ci <-gls.ci(AAmass$lograt, AAmass$logmass,vcv(tree13))
lines(pGLS.mod5.ci$CI.plot$X,pGLS.mod5.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod5.ci$CI.plot$X,pGLS.mod5.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod5.pi <-gls.pi(AAmass$lograt, AAmass$logmass,vcv(tree13), 1)
lines(pGLS.mod5.pi$PI.plot$X,pGLS.mod5.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod5.pi$PI.plot$X,pGLS.mod5.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()


## For the cerebellum-to-cerebrum ratio, which we expected to not be correlated with body mass ##
tiff("./FOLDER/5.24.BM-PGLS-CCratio-bodymass-ci-pi.tiff", width = 8, height = 6, units = 'in', res = 300)
mod6 <- gls(CCrat ~ mass, data=CCmass, correlation=bm)
plot(mod6)
plot(CCmass$CCrat ~ CCmass$mass, xlab = "Body mass (log-transformed)", ylab = "Ratio cerebellum/cerebrum (log-transformed)")
abline(mod6)
abline(a=1.195468, b=1.0, col= "blue") #MANUAL INPUT!! Is this correct?? Don't think so. B is incorrect (negative logged values).
# Confidence intervals
pGLS.mod6.ci <-gls.ci(CCmass$CCrat, CCmass$mass,vcv(tree))
lines(pGLS.mod6.ci$CI.plot$X,pGLS.mod6.ci$CI.plot$Lower2.5,lty=2)
lines(pGLS.mod6.ci$CI.plot$X,pGLS.mod6.ci$CI.plot$Upper2.5,lty=2)
# Prediction intervals
pGLS.mod6.pi <-gls.pi(CCmass$CCrat, CCmass$mass,vcv(tree), 1)
lines(pGLS.mod6.pi$PI.plot$X,pGLS.mod6.pi$PI.plot$Lower2.5,lty=2, col = "red")
lines(pGLS.mod6.pi$PI.plot$X,pGLS.mod6.pi$PI.plot$Upper2.5,lty=2, col = "red")
graphics.off()



## Apes and non-apes
CCmass.ape <- CCmass %>% dplyr::filter(Hominoidea == "Ape")
CCmass.nonape <- CCmass %>% dplyr::filter(Hominoidea == "Nonape")
AAmass.ape <- AAmass %>% dplyr::filter(Hominoidea == "Ape")
AAmass.nonape <- AAmass %>% dplyr::filter(Hominoidea == "Nonape")

reg.mod1.MASS.ape <- lm(MedianCerebroCerebellar ~ mass, data = CCmass.ape)
reg.mod1.MASS.nonape <- lm(MedianCerebroCerebellar ~ mass, data = CCmass.nonape)
reg.mod2.MASS.ape <- lm(MedianCerebellarCrus ~ mass, data = AAmass.ape)
reg.mod2.MASS.nonape <- lm(MedianCerebellarCrus ~ mass, data = AAmass.nonape)

summary(reg.mod1.MASS.ape)
summary(reg.mod1.MASS.nonape)
summary(reg.mod2.MASS.ape)
summary(reg.mod2.MASS.nonape)

tiff("./5.phylogenetic_regressions/5.24a.regularRegression-CCmass-apes.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplot(CCmass, aes(y=MedianCerebroCerebellar, x=mass, colour = Hominoidea)) + 
  geom_point(na.rm=T, aes(color=CCmass$Hominoidea)) +
  geom_smooth(method="lm", na.rm=T, se=F, formula=y~x, aes(color=CCmass$Hominoidea)) +
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

tiff("./5.phylogenetic_regressions/5.25a.regularRegression-AAmass-apes.tiff", width = 8, height = 8, units = 'in', res = 300)
ggplot(AAmass, aes(y=MedianCerebellarCrus, x=mass, colour = Hominoidea)) + 
  geom_point(na.rm=T, aes(color=AAmass$Hominoidea)) +
  geom_smooth(method="lm", na.rm=T, se=F, formula=y~x, aes(color=AAmass$Hominoidea)) +
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