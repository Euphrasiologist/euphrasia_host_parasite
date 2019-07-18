## Are there parasite-host interactions in Euphrasia? ##
# Created: 28.6.19 by Max Brown #
# Major update: 17.7.19 # 

## Libraries needed ##

library(data.table)
library(ggplot2)
library(MCMCglmm)
library(lme4)
library(VCVglmm)

# VCVglmm function not exporting properly... so rewrite here
Solapply <- function (model, FUN = posterior_mode, ...) {
  if (attributes(model)$class != "MCMCglmm") {
    stop("Object not of class MCMCglmm")
  }
  if (dim(model$Sol)[2] <= model$Fixed$nfl) {
    stop("Re-run the model with parameter option pr = TRUE")
  }
  posterior_mode <- function(x, adjust = 0.1, ...) {
    find.mode <- function(x, adjust, ...) {
      dx <- density(x, adjust = adjust, ...)
      dx$x[which.max(dx$y)]
    }
    apply(as.matrix(x), 2, find.mode, adjust = adjust, ...)
  }
  FUN <- match.fun(FUN)
  if (identical(FUN, HPDinterval)) {
    temp <- setDT(as.data.frame(HPDinterval(model$Sol)), 
                  keep.rownames = TRUE)[-c(1:model$Fixed$nfl)]
    temp$posterior.mode <- as.data.frame(apply(model$Sol[, 
                                                         -c(1:model$Fixed$nfl)], 2, posterior_mode))[, 1]
    colnames(temp) <- c("Variable", "lowerHPD", "upperHPD", 
                        "Posterior Mode")
    temp$Group <- gsub("\\..*", "", temp$Variable)
    return(temp)
  }
  temp <- setDT(as.data.frame(apply(model$Sol[, -c(1:model$Fixed$nfl)], 
                                    2, FUN)), keep.rownames = TRUE)
  colnames(temp) <- c("Variable", "Grouped_Value")
  temp$Group <- gsub("\\..*", "", temp$Variable)
  return(temp)
}


##### Part 1: Data tidying #####
expt2dat <- fread("./Data/Many_species/Experiment2firstflowering1.csv")
expt2key <- fread("./Data/Many_species/Experiment2Species.csv")

# change column name
colnames(expt2dat)[1] <- "Unique_ID"

expt2dat.2 <- expt2key[,-c("Notes", "Euphrasia dead", "Host dead")][expt2dat[,-c("Aborted node",
                                                                                 "Leaf circularity index",
                                                                                 "Capsule length",
                                                                                 "Capsule width",
                                                                                 "Calyx length",
                                                                                 "Calyx ratio",
                                                                                 "Notes")], on = "Unique_ID"]

# how many vigursii ended up flowering?
table(expt2dat.2$Host_code[expt2dat.2$`Flower colour` == "V" & expt2dat.2$Euphrasia_sp == "Euphrasia vigursii"])
# new column of Euphrasia species, including tetraquetra
expt2dat.2[, `:=`(Euphrasia_sp2 = ifelse(`Flower colour` != "V" & Euphrasia_sp == "Euphrasia vigursii", 
                                         "Euphrasia tetraquetra", Euphrasia_sp))]
# give Euphrasia tetraquetra its own population level.
expt2dat.2[Euphrasia_sp2 == "Euphrasia tetraquetra"]$Population <- "T1761"
# and turn species into a factor
expt2dat.2$Euphrasia_sp2 <- as.factor(expt2dat.2$Euphrasia_sp2)



##### Plot 1: Visualise quickly what is going on! #####
# make sure height is numeric
expt2dat.2$Height <- as.numeric(expt2dat.2$Height)
# remove NA's introduced by coercion
expt2dat.2 <- expt2dat.2[!is.na(expt2dat.2$Height)]
expt2dat.2 <- expt2dat.2[!is.na(expt2dat.2$Host_code)]

# for experiment A (from here-on out). Replace "A" with "B" for experiment B
# for a quick idea of how species are holistically reacting to each host.
plot_1.1 <- ggplot(expt2dat.2[grepl("A", Unique_ID)], aes(x=Host_code, y = Height))+
  geom_jitter(width = 0.25)+
  facet_wrap(~Euphrasia_sp2)+
  theme_bw()+
  stat_summary(fun.y = mean, col="red", geom = "point")+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))+
  ylab(label = "Height (mm)")+
  xlab(label = "Host species")

ggsave(filename = "./Figures/Many_species/quick_visualisation", plot = plot_1.1, 
       device = "pdf", width = 6, height = 5, units = "in")

##### Part 2: Add data on reproductive nodes at the end of the season #####

# read in nodes data
rnodes <- fread("./Data/Many_species/REPRODUCTIVENODES.csv")
# merge with data generated at first flowering. Merged so that UniqueID's that survived to end of season are preserved 
# (not all that first flowered...)
rnodes2 <- expt2dat.2[rnodes, on = "Unique_ID"]
rnodes2$Euphrasia_sp2 <- factor(rnodes2$Euphrasia_sp2)


##### Part 3: Model of reproductive nodes as a function of Euphrasia species and look for interaction #####

# prior for interaction model.
prior.manysp <- list(R=list(V=diag(1), nu=0.002), 
                     G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                            G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                            G3=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

# create the interaction factors
# between species of Euphrasia and host
rnodes2$Host_given_Euphrasia <- interaction(rnodes2$Euphrasia_sp2, rnodes2$Host_code)
# and between population of Euphrasia and host
rnodes2$Host_given_population <- interaction(rnodes2$Population, rnodes2$Host_code)

# remove any NA's
rnodes3 <- rnodes2[complete.cases(Euphrasia_sp2)]

manysp.2<-MCMCglmm(Reproductive_nodes ~ Euphrasia_sp2,
                   random = ~ Host_code + Host_given_Euphrasia + Host_given_population,
                   family = "poisson",
                   prior=prior.manysp,
                   data = rnodes3,
                   nitt = 13000*8,
                   burnin = 3000*8,
                   thin = 10*8,
                   pr=TRUE,
                   verbose = TRUE)

# correlation in host effects
posterior.mode(manysp.2$VCV[,"Host_code"]/(manysp.2$VCV[,"Host_code"]+manysp.2$VCV[,"Host_given_population"]))

write.csv(x = summary(manysp.2)$solutions,
          file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/Model_solutions.csv")

write.csv(x = data.table(HOST = MCMCReppois(mod = manysp.2, y = "Host_code"),
                         HOST_GIVEN_EUPHRASIA = MCMCReppois(mod = manysp.2, y = "Host_given_Euphrasia"),
                         HOST_GIVEN_POPULATION = MCMCReppois(mod = manysp.2, y = "Host_given_population"), keep.rownames = TRUE),
          file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/Variance_Components.csv")

# visualise the variance components
VCVdensity(manysp.2)+xlim(0, 1)

# so there are interactions here, let's dig further
# so this is great, we have the posterior modes for population:host interactions and host
write.csv(x = Solapply(manysp.2)[order(Grouped_Value)][Group %in% c("Host_code", "Host_given_population")][, .(Variable = Variable,
                                                                                                               Grouped_Value = exp(Grouped_Value))],
          file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/Posterior_Modes.csv")
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
  ), file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/LRTs_of_models.csv"
)


##### Plot 2: Posterior Modes of the interaction model #####

plot_2.1 <- Solapply(manysp.2, coda::HPDinterval)[order(`Posterior Mode`)][Group %in% c("Host_code", "Host_given_population")] %>%
  ggplot(aes(x = reorder(Variable, `Posterior Mode`), y = `Posterior Mode`))+
  geom_point()+
  geom_errorbar(aes(ymin = lowerHPD, ymax = upperHPD))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(filename = "./Figures/Many_species/posterior_interaction_modes_mod_1", plot = plot_2.1, 
       device = "pdf", width = 21, height = 6, units = "in")

##### Plot 3: Raw data for the manuscript, means and SE's #####

plot_3.1<- rnodes3[, .(mean = mean(Reproductive_nodes),
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

ggsave(filename = "./Figures/Many_species/population_cum_nodes", plot = plot_3.1, 
       device = "pdf", width = 10, height = 6, units = "in")

# just Euphrasia vigursii and tetraquetra, ready for model comparison

plot_3.2 <- rnodes3[, .(mean = mean(Reproductive_nodes),
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

ggsave(filename = "./Figures/Many_species/population_cum_nodes_vig_tet", plot = plot_3.2, 
       device = "pdf", width = 10, height = 6, units = "in")

##### Part 4: Model two populations from the same location for differences #####

# subset data to investigate only tetraquetra and vigursii
rnodes4 <- rnodes3[Euphrasia_sp2 %in% c("Euphrasia tetraquetra", "Euphrasia vigursii")]

prior.manysp.3 <- list(R=list(V=diag(1), nu=0.002), 
                       G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                              G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

# remove unwanted factor levels
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
                   pr=TRUE,
                   verbose = TRUE)

write.csv(x = summary(manysp.3)$solutions,
          file = "./Data/Many_species/Model_outputs/Tet_vs_Vig/Model_Solutions.csv")

write.csv(x = data.table(HOST = MCMCReppois(mod = manysp.3, y = "Host_code"),
                         HOST_GIVEN_EUPHRASIA = MCMCReppois(mod = manysp.3, y = "Host_given_Euphrasia"), keep.rownames = TRUE),
          file = "./Data/Many_species/Model_outputs/Tet_vs_Vig/Variance_Components.csv")

VCVdensity(manysp.3)+xlim(0, 1)

# posterior modes?

write.csv(x = Solapply(manysp.3)[order(Grouped_Value)][, .(Variable = Variable,
                                                           Grouped_Value = exp(Grouped_Value))],
          file = "./Data/Many_species/Model_outputs/Tet_vs_Vig/Posterior_Modes.csv")

# significance of random effects
length(unique(rnodes4$Obs))

# full model
manysp.3LR1 <- glmer(Reproductive_nodes ~ 1 + (1 | Host_code) + (1 | Host_given_Euphrasia) + (1 | Obs),
                     family = "poisson", data = rnodes4)
# is Host_code significant?
manysp.3LR2 <- glmer(Reproductive_nodes ~ 1 + (1 | Host_given_Euphrasia) + (1 | Obs),
                     family = "poisson", data = rnodes4)
# is Host_given_Euphrasia significant?
manysp.3LR3 <- glmer(Reproductive_nodes ~ 1 + (1 | Host_code) + (1 | Obs),
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
    HOST_GIVEN_EUPHRASIA = anova(manysp.3LR1, manysp.3LR3), keep.rownames = TRUE
  ), file = "./Data/Many_species/Model_outputs/Tet_vs_Vig/LRTs_of_models.csv"
)

##### Plot 4: Two population model plot of posterior modes #####

plot_4.1 <- Solapply(manysp.3, HPDinterval)[order(`Posterior Mode`)][Group %in% c("Host_code", "Host_given_Euphrasia")] %>%
  ggplot(aes(x = reorder(Variable, `Posterior Mode`), y = `Posterior Mode`))+
  geom_point()+
  geom_errorbar(aes(ymin = lowerHPD, ymax = upperHPD))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  xlab(label = "Factor")

ggsave(filename = "./Figures/Many_species/vig_tet_interation_plot", plot = plot_4.1, 
       device = "pdf", width = 10, height = 6, units = "in")

##### Part 5: Differences in E.anglica and E. micrantha performance #####

# subset data for micrantha
rnodes_micrantha <- rnodes3[Euphrasia_sp2 %in% c("Euphrasia micrantha")]

# prior for micrantha model
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
                   pr=TRUE,
                   verbose = TRUE)

# looks like Host more important thank host given pop
write.csv(x = summary(manysp.4)$solutions,
          file = "./Data/Many_species/Model_outputs/Micrantha/Model_Solutions.csv")

write.csv(x = specify_decimal(data.table(HOST = MCMCReppois(mod = manysp.4, y = "Host_code"),
                                         HOST_GIVEN_POPULATION = MCMCReppois(mod = manysp.4, y = "Host_given_population")), 6),
          file = "./Data/Many_species/Model_outputs/Micrantha/Variance_Components.csv")

VCVdensity(manysp.4)+xlim(0, 1)

# posterior modes?

write.csv(x = Solapply(manysp.4)[order(Grouped_Value)][, .(Variable = Variable,
                                                           Grouped_Value = exp(Grouped_Value))],
          file = "./Data/Many_species/Model_outputs/Micrantha/Posterior_Modes.csv")

# significance of random effects
# need to add Obs
rnodes_micrantha$Obs <- as.factor(1:nrow(rnodes_micrantha)) 

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
  ), file = "./Data/Many_species/Model_outputs/Micrantha/LRTs_of_models.csv"
)

# subset data for anglica and repeat
rnodes_anglica <- rnodes3[Euphrasia_sp2 %in% c("Euphrasia anglica")]

# prior for anglica
prior.manysp.5 <- list(R=list(V=diag(1), nu=0.002), 
                       G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                              G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

rnodes_anglica$Euphrasia_sp2 <- factor(rnodes_anglica$Euphrasia_sp2)
rnodes_anglica$Host_code <- factor(rnodes_anglica$Host_code)
rnodes_anglica$Host_given_Euphrasia <- factor(rnodes_anglica$Host_given_Euphrasia)
rnodes_anglica$Host_given_population <- factor(rnodes_anglica$Host_given_population)

# 
manysp.5<-MCMCglmm(Reproductive_nodes ~ Population,
                   random = ~ Host_code + Host_given_population,
                   family = "poisson",
                   prior=prior.manysp.5,
                   data = rnodes_anglica,
                   nitt = 13000*10,
                   burnin = 3000*10,
                   thin = 10*10,
                   pr=TRUE,
                   verbose = TRUE)

# solutions
write.csv(x = summary(manysp.5)$solutions,
          file = "./Data/Many_species/Model_outputs/Anglica/Model_Solutions.csv")

write.csv(x = specify_decimal(data.table(HOST = MCMCReppois(mod = manysp.5, y = "Host_code"),
                                         HOST_GIVEN_POPULATION = MCMCReppois(mod = manysp.5, y = "Host_given_population")), 6),
          file = "./Data/Many_species/Model_outputs/Anglica/Variance_Components.csv")

VCVdensity(manysp.4)+xlim(0, 1)

# posterior modes?

write.csv(x = Solapply(manysp.5)[order(Grouped_Value)][, .(Variable = Variable,
                                                           Grouped_Value = exp(Grouped_Value))],
          file = "./Data/Many_species/Model_outputs/Anglica/Posterior_Modes.csv")

# significance of random effects
# need to add Obs
rnodes_anglica$Obs <- as.factor(1:nrow(rnodes_anglica)) 

# full model
manysp.5LR1 <- glmer(Reproductive_nodes ~ Population + (1 | Host_code) + (1 | Host_given_population) + (1 | Obs),
                     family = "poisson", data = rnodes_anglica)
# is Host_code significant?
manysp.5LR2 <- glmer(Reproductive_nodes ~ Population + (1 | Host_given_population) + (1 | Obs),
                     family = "poisson", data = rnodes_anglica)
# is Host_given_population significant?
manysp.5LR3 <- glmer(Reproductive_nodes ~ Population + (1 | Host_code) + (1 | Obs),
                     family = "poisson", data = rnodes_anglica)

# Host_code is not significant
anova(manysp.5LR1, manysp.5LR2)
# Host_given_code_Euphrasia is not significant
anova(manysp.5LR1, manysp.5LR3)

write.csv(
  x = data.table(
    # Host_code is not significant
    HOST_CODE = anova(manysp.5LR1, manysp.5LR2),
    # Host_given_population is significant
    HOST_GIVEN_POPULATION = anova(manysp.5LR1, manysp.5LR3), keep.rownames = TRUE
  ), file = "./Data/Many_species/Model_outputs/Anglica/LRTs_of_models.csv"
)

##### Plot 5: E. anglica and E. micrantha interactions #####

# Euphrasia micrantha
plot_5.1 <- Solapply(manysp.4, HPDinterval)[order(`Posterior Mode`)][Group %in% c("Host_code", "Host_given_population")] %>%
  ggplot(aes(x = reorder(Variable, `Posterior Mode`), y = `Posterior Mode`))+
  geom_point()+
  geom_errorbar(aes(ymin = lowerHPD, ymax = upperHPD))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  xlab(label = "Factor")

ggsave(filename = "./Figures/Many_species/micrantha_interation", plot = plot_5.1, 
       device = "pdf", width = 10, height = 6, units = "in")

# Euphrasia anglica
plot_5.2 <- Solapply(manysp.5, HPDinterval)[order(`Posterior Mode`)][Group %in% c("Host_code", "Host_given_population")] %>%
  ggplot(aes(x = reorder(Variable, `Posterior Mode`), y = `Posterior Mode`))+
  geom_point()+
  geom_errorbar(aes(ymin = lowerHPD, ymax = upperHPD))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  xlab(label = "Factor")

ggsave(filename = "./Figures/Many_species/anglica_interation", plot = plot_5.2, 
       device = "pdf", width = 10, height = 6, units = "in")
