# Max Brown 28.6.19

# tidying up the Euphrasiaspeciesscript.R script

# libraries

library(data.table)
library(ggplot2)
library(MCMCglmm)
library(lme4)
library(VCVglmm)

##### Are there parasite-host interactions in Euphrasia? #####
##### Part 1: Data tidying #####
expt2dat <- fread("/Users/mbrown/OneDrive - University of Edinburgh/EUPHRASIA EXPERIMENT 2 HOSTS AND SPECIES/Experiment2firstflowering1.csv")
expt2key <- fread("/Users/mbrown/OneDrive - University of Edinburgh/EUPHRASIA EXPERIMENT 2 HOSTS AND SPECIES/Experiment2Species.csv")

colnames(expt2dat)[1] <- "Unique_ID"

expt2dat.2 <- expt2key[,-c(7:9)][expt2dat[,-c(15:21)], on = "Unique_ID"]

# how many vigursii ended up flowering?
table(expt2dat.2$Host_code[expt2dat.2$`Flower colour` == "V" & expt2dat.2$Euphrasia_sp == "Euphrasia vigursii"])
# new column of Euphrasia species
expt2dat.2[, `:=`(Euphrasia_sp2 = ifelse(`Flower colour` != "V" & Euphrasia_sp == "Euphrasia vigursii", 
                                         "Euphrasia tetraquetra", Euphrasia_sp))]
# new column of population
# all data
expt2dat.2[Euphrasia_sp2 == "Euphrasia tetraquetra"]$Population <- "T1761"

# and turn it into a factor
expt2dat.2$Euphrasia_sp2 <- as.factor(expt2dat.2.test$Euphrasia_sp2)



##### Plot 1: Visualise quickly what is going on! #####
expt2dat.2$Height <- as.numeric(expt2dat.2$Height)
expt2dat.2 <- expt2dat.2[!is.na(expt2dat.2$Height)]
expt2dat.2 <- expt2dat.2[!is.na(expt2dat.2$Host_code)]

# for experiment A (from here-on out). Replace "A" with "B" for experiment B
ggplot(expt2dat.2[grepl("A", Unique_ID)], aes(x=Host_code, y = Height))+
  geom_point()+
  facet_wrap(~Euphrasia_sp2)+
  theme_bw()+
  stat_summary(fun.y = mean, col="red", geom = "point")+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))+
  ylab(label = "Height (mm)")+
  xlab(label = "Host species")

##### Part 2: Add data on reproductive nodes at the end of the season #####

# read in nodes data
rnodes <- fread("/Users/mbrown/OneDrive - University of Edinburgh/EUPHRASIA EXPERIMENT 2 HOSTS AND SPECIES/REPRODUCTIVENODES.csv")
# merge with data generated at first flowering. Merged so that UniqueID's that survived to end of season are preserved 
# (not all that first flowered...)
rnodes2 <- expt2dat.2[rnodes, on = "Unique_ID"]
rnodes2$Euphrasia_sp2 <- factor(rnodes2$Euphrasia_sp2)


##### Part 3: Model of reproductive nodes as a function of Euphrasia species and look for interaction #####

prior.manysp <- list(R=list(V=diag(1), nu=0.002), 
                     G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                            G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                            G3=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))


rnodes2$Host_given_Euphrasia <- interaction(rnodes2$Euphrasia_sp2, rnodes2$Host_code)
rnodes2$Host_given_population <- interaction(rnodes2$Population, rnodes2$Host_code)

rnodes3 <- rnodes2[complete.cases(Euphrasia_sp2)]

manysp.2<-MCMCglmm(Reproductive_nodes ~ Euphrasia_sp2,
                   random = ~ Host_code + Host_given_Euphrasia + Host_given_population,
                   family = "poisson",
                   prior=prior.manysp,
                   data = rnodes3,
                   nitt = 13000*8,
                   burnin = 3000*8,
                   thin = 10*8,
                   pr=T,
                   verbose = T)

# correlation in host effects
posterior.mode(manysp.2$VCV[,"Host_code"]/(manysp.2$VCV[,"Host_code"]+manysp.2$VCV[,"Host_given_population"]))

write.csv(x = summary(manysp.2)$solutions,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Host_parasite_interaction/Model_solutions.csv")

write.csv(x = data.table(HOST = MCMCReppois(mod = manysp.2, y = "Host_code"),
                         HOST_GIVEN_EUPHRASIA = MCMCReppois(mod = manysp.2, y = "Host_given_Euphrasia"),
                         HOST_GIVEN_POPULATION = MCMCReppois(mod = manysp.2, y = "Host_given_population"), keep.rownames = TRUE),
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Host_parasite_interaction/Variance_Components.csv")

# visualise the variance components
VCVdensity(manysp.2)+xlim(0, 1)
# so there are interactions here, let's dig further

# so this is great, we have the posterior modes for population:host interactions and host
write.csv(x = Solapply(manysp.2)[order(Grouped_Value)][Group %in% c("Host_code", "Host_given_population")][, .(Variable = Variable,
                                                                                                               Grouped_Value = exp(Grouped_Value))],
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Host_parasite_interaction/Posterior_Modes.csv")
# significance of random effects
# need to add Obs
rnodes3$Obs <- as.factor(1:nrow(rnodes3)) # maybe do not add...

# full model
manysp.2LR1 <- glmer(Reproductive_nodes ~ Euphrasia_sp2 + (1 | Host_code) + (1 | Host_given_Euphrasia) + (1 | Host_given_population),
                     family = "poisson", data = rnodes3)
# is Host_code significant?
manysp.2LR2 <- glmer(Reproductive_nodes ~ Euphrasia_sp2 + (1 | Host_given_Euphrasia) + (1 | Host_given_population),
                     family = "poisson", data = rnodes3)
# is Host_given_Euphrasia significant?
manysp.2LR3 <- glmer(Reproductive_nodes ~ Euphrasia_sp2 + (1 | Host_code) + (1 | Host_given_population),
                     family = "poisson", data = rnodes3)
# is Host_given_population significant?
manysp.2LR4 <- glmer(Reproductive_nodes ~ Euphrasia_sp2 + (1 | Host_code) +  (1 | Host_given_Euphrasia),
                     family = "poisson", data = rnodes3)

write.csv(
  x = data.table(
    # Host_code is significant
    HOST_CODE = anova(manysp.2LR1, manysp.2LR2),
    # Host_given_code_Euphrasia is not significant
    HOST_GIVEN_EUPHRASIA = anova(manysp.2LR1, manysp.2LR3),
    # Host_given_population is highly significant
    HOST_GIVEN_POPULATION = anova(manysp.2LR1, manysp.2LR4), keep.rownames = TRUE
  ), file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Host_parasite_interaction/LRTs_of_models.csv"
)


##### Plot 2: Posterior Modes of the interaction model #####
# and we can plot it easily thanks to Solapply and some dplyr magic.
Solapply(manysp.2, HPDinterval)[order(`Posterior Mode`)][Group %in% c("Host_code", "Host_given_population")] %>%
  ggplot(aes(x = reorder(Variable, `Posterior Mode`), y = `Posterior Mode`))+
  geom_point()+
  geom_errorbar(aes(ymin = lowerHPD, ymax = upperHPD))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

##### Plot 3: Raw data for the manuscript, means and SE's #####

rnodes3[, .(mean = mean(Reproductive_nodes),
            sem = sd(Reproductive_nodes)/sqrt(.N),
            N = .N), by = c("Euphrasia_sp2", "Host_code", "Population")] %>% #[Population != "M1767"]
  
  ggplot(aes(x = reorder(Host_code, mean) , y = mean))+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem, group=Euphrasia_sp2), position = position_dodge(width = 0.9), width=0.4)+
  geom_point(aes(colour = Euphrasia_sp2), position = position_dodge(width = 0.9), size=3)+
  facet_wrap(~Population, scales = "free_y")+
  theme_bw()+ theme(strip.text.x = element_text(size=20),
                    strip.background = element_rect(colour="white", fill="white"),
                    axis.line.x = element_line(colour = "black"),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.text.x = element_text(angle = 60, hjust = 1),
                    axis.title.x.bottom = element_text(size = 20),
                    axis.title.y.left = element_text(size = 20),
                    legend.title = element_text(size = 20))+
  xlab(label = "Host Species")+
  ylab(label = "Mean reproductive nodes at end of season")+
  scale_colour_discrete(name = "Euphrasia species")

# just Euphrasia vigursii and tetraquetra, ready for model comparison

rnodes3[, .(mean = mean(Reproductive_nodes),
            sem = sd(Reproductive_nodes)/sqrt(.N),
            N = .N), by = c("Euphrasia_sp2", "Host_code", "Population")][Euphrasia_sp2 %in% c("Euphrasia tetraquetra", "Euphrasia vigursii")] %>%
  
  ggplot(aes(x = reorder(Host_code, mean) , y = mean))+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem, group=Euphrasia_sp2), position = position_dodge(width = 0.9), width=0.4)+
  geom_point(aes(colour = Euphrasia_sp2), position = position_dodge(width = 0.9), size=3)+
  facet_wrap(~Population)+
  theme_bw()+ theme(strip.text.x = element_text(size=20),
                    strip.background = element_rect(colour="white", fill="white"),
                    axis.line.x = element_line(colour = "black"),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.text.x = element_text(angle = 60, hjust = 1),
                    axis.title.x.bottom = element_text(size = 20),
                    axis.title.y.left = element_text(size = 20),
                    legend.title = element_text(size = 20))+
  xlab(label = "Host Species")+
  ylab(label = "Mean reproductive nodes at end of season")+
  scale_colour_discrete(name = "Euphrasia species")


##### Part 4: Model two populations from the same location for differences #####

# subset data

rnodes4 <- rnodes3[Euphrasia_sp2 %in% c("Euphrasia tetraquetra", "Euphrasia vigursii")]

prior.manysp.3 <- list(R=list(V=diag(1), nu=0.002), 
                       G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                              G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

rnodes4$Euphrasia_sp2 <- factor(rnodes4$Euphrasia_sp2)
rnodes4$Host_code <- factor(rnodes4$Host_code)
rnodes4$Host_given_Euphrasia <- factor(rnodes4$Host_given_Euphrasia)
rnodes4$Host_given_population <- factor(rnodes4$Host_given_population)

manysp.3<-MCMCglmm(Reproductive_nodes ~ Euphrasia_sp2,
                   random = ~ Host_code + Host_given_Euphrasia,
                   family = "poisson",
                   prior=prior.manysp.3,
                   data = rnodes4,
                   nitt = 13000*10,
                   burnin = 3000*10,
                   thin = 10*10,
                   pr=T,
                   verbose = T)

write.csv(x = summary(manysp.3)$solutions,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Tet_vs_Vig/Model_Solutions.csv")

write.csv(x = data.table(HOST = MCMCReppois(mod = manysp.3, y = "Host_code"),
                         HOST_GIVEN_EUPHRASIA = MCMCReppois(mod = manysp.3, y = "Host_given_Euphrasia"), keep.rownames = TRUE),
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Tet_vs_Vig/Variance_Components.csv")

VCVdensity(manysp.3)+xlim(0, 1)

# posterior modes?

write.csv(x = Solapply(manysp.3)[order(Grouped_Value)][, .(Variable = Variable,
                                                           Grouped_Value = exp(Grouped_Value))],
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Tet_vs_Vig/Posterior_Modes.csv")

# significance of random effects
# need to add Obs
rnodes3$Obs <- as.factor(1:nrow(rnodes3)) # maybe do not add...

# full model
manysp.3LR1 <- glmer(Reproductive_nodes ~ 1 + (1 | Host_code) + (1 | Host_given_Euphrasia),
                     family = "poisson", data = rnodes4)
# is Host_code significant?
manysp.3LR2 <- glmer(Reproductive_nodes ~ 1 + (1 | Host_given_Euphrasia),
                     family = "poisson", data = rnodes4)
# is Host_given_Euphrasia significant?
manysp.3LR3 <- glmer(Reproductive_nodes ~ 1 + (1 | Host_code),
                     family = "poisson", data = rnodes4)

# Host_code is significant
anova(manysp.3LR1, manysp.3LR2)
# Host_given_Euphrasia is not significant
anova(manysp.3LR1, manysp.3LR3)

write.csv(
  x = data.table(
    # Host_code is significant
    HOST_CODE = anova(manysp.3LR1, manysp.3LR2),
    # Host_given_code_Euphrasia is not significant
    HOST_GIVEN_EUPHRASIA = anova(manysp.3LR1, manysp.3LR3),
    # Host_given_population is not significant (yay)
    HOST_GIVEN_POPULATION = anova(manysp.3LR1, manysp.3LR4), keep.rownames = TRUE
  ), file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Tet_vs_Vig/LRTs_of_models.csv"
)

##### Plot 4: Two population model plot of posterior modes #####

Solapply(manysp.3, HPDinterval)[order(`Posterior Mode`)][Group %in% c("Host_code", "Host_given_population")] %>%
  ggplot(aes(x = reorder(Variable, `Posterior Mode`), y = `Posterior Mode`))+
  geom_point()+
  geom_errorbar(aes(ymin = lowerHPD, ymax = upperHPD))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

##### Part 5: Differences in E.anglica and E. micrantha performance #####

# subset data for micrantha

rnodes_micrantha <- rnodes3[Euphrasia_sp2 %in% c("Euphrasia micrantha")]

prior.manysp.4 <- list(R=list(V=diag(1), nu=0.002), 
                       G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                              G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

rnodes_micrantha$Euphrasia_sp2 <- factor(rnodes_micrantha$Euphrasia_sp2)
rnodes_micrantha$Host_code <- factor(rnodes_micrantha$Host_code)
rnodes_micrantha$Host_given_Euphrasia <- factor(rnodes_micrantha$Host_given_Euphrasia)
rnodes_micrantha$Host_given_population <- factor(rnodes_micrantha$Host_given_population)

manysp.4<-MCMCglmm(Reproductive_nodes ~ Population,
                   random = ~ Host_code + Host_given_population,
                   family = "poisson",
                   prior=prior.manysp.4,
                   data = rnodes_micrantha,
                   nitt = 13000*10,
                   burnin = 3000*10,
                   thin = 10*10,
                   pr=T,
                   verbose = T)
# looks like Host more important thank host given pop


write.csv(x = summary(manysp.4)$solutions,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Micrantha/Model_Solutions.csv")

write.csv(x = specify_decimal(data.table(HOST = MCMCReppois(mod = manysp.4, y = "Host_code"),
                                         HOST_GIVEN_POPULATION = MCMCReppois(mod = manysp.4, y = "Host_given_population")), 6),
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Micrantha/Variance_Components.csv")

VCVdensity(manysp.4)+xlim(0, 1)

# posterior modes?

write.csv(x = Solapply(manysp.4)[order(Grouped_Value)][, .(Variable = Variable,
                                                           Grouped_Value = exp(Grouped_Value))],
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Micrantha/Posterior_Modes.csv")

# significance of random effects
# need to add Obs
rnodes_micrantha$Obs <- as.factor(1:nrow(rnodes_micrantha)) # maybe do not add...

# full model
manysp.4LR1 <- glmer(Reproductive_nodes ~ Population + (1 | Host_code) + (1 | Host_given_population) + (1 | Obs),
                     family = "poisson", data = rnodes_micrantha)
# is Host_code significant?
manysp.4LR2 <- glmer(Reproductive_nodes ~ Population + (1 | Host_given_population) + (1 | Obs),
                     family = "poisson", data = rnodes_micrantha)
# is Host_given_population significant?
manysp.4LR3 <- glmer(Reproductive_nodes ~ Population + (1 | Host_code) + (1 | Obs),
                     family = "poisson", data = rnodes_micrantha)

# Host_code is significant
anova(manysp.4LR1, manysp.4LR2)
# Host_given_code_Euphrasia is not significant
anova(manysp.4LR1, manysp.4LR3)

write.csv(
  x = data.table(
    # Host_code is not significant
    HOST_CODE = anova(manysp.4LR1, manysp.4LR2),
    # Host_given_population is significant
    HOST_GIVEN_POPULATION = anova(manysp.4LR1, manysp.4LR3), keep.rownames = TRUE
  ), file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Micrantha/LRTs_of_models.csv"
)

# subset data for anglica and repeat - NOTE NOT YET FINISHED
rnodes_anglica <- rnodes3[Euphrasia_sp2 %in% c("Euphrasia anglica")]

prior.manysp.5 <- list(R=list(V=diag(1), nu=0.002), 
                       G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                              G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

rnodes_anglica$Euphrasia_sp2 <- factor(rnodes_anglica$Euphrasia_sp2)
rnodes_anglica$Host_code <- factor(rnodes_anglica$Host_code)
rnodes_anglica$Host_given_Euphrasia <- factor(rnodes_anglica$Host_given_Euphrasia)
rnodes_anglica$Host_given_population <- factor(rnodes_anglica$Host_given_population)

manysp.5<-MCMCglmm(Reproductive_nodes ~ Population,
                   random = ~ Host_code + Host_given_population,
                   family = "poisson",
                   prior=prior.manysp.5,
                   data = rnodes_anglica,
                   nitt = 13000*10,
                   burnin = 3000*10,
                   thin = 10*10,
                   pr=T,
                   verbose = T)
# looks like Host more important than host given pop


write.csv(x = summary(manysp.5)$solutions,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Micrantha/Model_Solutions.csv")

write.csv(x = specify_decimal(data.table(HOST = MCMCReppois(mod = manysp.5, y = "Host_code"),
                                         HOST_GIVEN_POPULATION = MCMCReppois(mod = manysp.5, y = "Host_given_population")), 6),
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Micrantha/Variance_Components.csv")

VCVdensity(manysp.4)+xlim(0, 1)

# posterior modes?

write.csv(x = Solapply(manysp.4)[order(Grouped_Value)][, .(Variable = Variable,
                                                           Grouped_Value = exp(Grouped_Value))],
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Micrantha/Posterior_Modes.csv")

# significance of random effects
# need to add Obs
rnodes_micrantha$Obs <- as.factor(1:nrow(rnodes_micrantha)) # maybe do not add...

# full model
manysp.4LR1 <- glmer(Reproductive_nodes ~ Population + (1 | Host_code) + (1 | Host_given_population) + (1 | Obs),
                     family = "poisson", data = rnodes_micrantha)
# is Host_code significant?
manysp.4LR2 <- glmer(Reproductive_nodes ~ Population + (1 | Host_given_population) + (1 | Obs),
                     family = "poisson", data = rnodes_micrantha)
# is Host_given_population significant?
manysp.4LR3 <- glmer(Reproductive_nodes ~ Population + (1 | Host_code) + (1 | Obs),
                     family = "poisson", data = rnodes_micrantha)

# Host_code is significant
anova(manysp.4LR1, manysp.4LR2)
# Host_given_code_Euphrasia is not significant
anova(manysp.4LR1, manysp.4LR3)

write.csv(
  x = data.table(
    # Host_code is not significant
    HOST_CODE = anova(manysp.4LR1, manysp.4LR2),
    # Host_given_population is significant
    HOST_GIVEN_POPULATION = anova(manysp.4LR1, manysp.4LR3), keep.rownames = TRUE
  ), file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Multiple_Euphrasia_sp/Micrantha/LRTs_of_models.csv"
)
