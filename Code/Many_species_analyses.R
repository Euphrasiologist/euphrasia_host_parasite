## Are there parasite-host interactions in Euphrasia? ##
# Created: 28.6.19 by Max Brown #
# Major update: 17.7.19 # 

## Libraries needed ##

library(data.table)
library(ggplot2)
library(MCMCglmm)
library(lme4)
library(VCVglmm)
library(dplyr)

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

# "not" in
`%!in%` <- function(x,y)!('%in%'(x,y))

# little function to convert poisson slopes and intercepts.
ablines <- function(intercept, slope, x1, x2){
  x <- seq(x1,x2,0.01)
  res <- exp(intercept + slope*x)
  res
}

# euphrasia anglica is actually only from one population...
# so maybe sink 


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
#write.csv(x = expt2dat.2, file = "./Data/Many_species/expt2dat.2.csv")
# how many vigursii ended up flowering?
table(expt2dat.2$Host_code[expt2dat.2$`Flower colour` == "V" & expt2dat.2$Euphrasia_sp == "Euphrasia vigursii"])
# new column of Euphrasia species, including tetraquetra
expt2dat.2[, `:=`(Euphrasia_sp2 = ifelse(`Flower colour` != "V" & Euphrasia_sp == "Euphrasia vigursii", 
                                         "Euphrasia tetraquetra", Euphrasia_sp))]
# give Euphrasia tetraquetra its own population level.
expt2dat.2[Euphrasia_sp2 == "Euphrasia tetraquetra"]$Population <- "T1761"
# and turn species into a factor
expt2dat.2$Euphrasia_sp2 <- as.factor(expt2dat.2$Euphrasia_sp2)
# only one population of Euphrasia anglica
unique(expt2dat.2$Population)
expt2dat.2[, Population2 := ifelse(test = Population %in% c("A1766", "A1763"),
       yes = "A1766",
       no = Population)]

expt2dat.2[is.na(Population2)]

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
rnodes2$Host_given_population <- NULL

# remove any NA's
rnodes3 <- rnodes2[complete.cases(Euphrasia_sp2)]

# subset only the micrantha population I used -> M1767
# non significant, but positive effect of host co-occurrence (increases above intercept by ~3 nodes.)

manysp.2<-MCMCglmm(Reproductive_nodes ~ Euphrasia_sp2,
                   random = ~ Host_code + Population + Host_given_Euphrasia,
                   family = "poisson",
                   prior=prior.manysp,
                   data = rnodes3,
                   nitt = 13000*8,
                   burnin = 3000*8,
                   thin = 10*8,
                   pr=TRUE,
                   verbose = TRUE)
summary(manysp.2)

# correlation in host effects
posterior.mode(manysp.2$VCV[,"Host_code"]/(manysp.2$VCV[,"Host_code"]+manysp.2$VCV[,"Host_given_Euphrasia"]))
HPDinterval(manysp.2$VCV[,"Host_code"]/(manysp.2$VCV[,"Host_code"]+manysp.2$VCV[,"Host_given_Euphrasia"]))
# and the inverse for completeness
posterior.mode(manysp.2$VCV[,"Host_given_Euphrasia"]/(manysp.2$VCV[,"Host_code"]+manysp.2$VCV[,"Host_given_Euphrasia"]))
HPDinterval(manysp.2$VCV[,"Host_given_Euphrasia"]/(manysp.2$VCV[,"Host_code"]+manysp.2$VCV[,"Host_given_Euphrasia"]))


write.csv(x = summary(manysp.2)$solutions,
          file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/Model_solutions.csv")

write.csv(x = data.table(HOST = MCMCReppois(mod = manysp.2, y = "Host_code"),
                         HOST_GIVEN_EUPHRASIA = MCMCReppois(mod = manysp.2, y = "Host_given_Euphrasia"),
                         POPULATION = MCMCReppois(mod = manysp.2, y = "Population"),
                         UNITS = MCMCReppois(mod = manysp.2, y = "units"), keep.rownames = TRUE),
          file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/Variance_Components.csv")

# visualise the variance components
VCVdensity(manysp.2)+xlim(0, 1)

# so there are interactions here, let's dig further
# so this is great, we have the posterior modes for population:host interactions and host
write.csv(x = Solapply(manysp.2)[order(Grouped_Value)][Group %in% c("Host_code", "Host_given_Euphrasia")][, .(Variable = Variable,
                                                                                                               Grouped_Value = exp(Grouped_Value))],
          file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/Posterior_Modes.csv")
# significance of random effects
# need to add Obs
rnodes3$Obs <- as.factor(1:nrow(rnodes3)) # maybe do not add...

# full model
manysp.2LR1 <- glmer(Reproductive_nodes ~ Euphrasia_sp2 + (1 | Host_code) + (1 | Host_given_Euphrasia) + (1 | Population) + (1|Obs),
                     family = "poisson", data = rnodes3)
# is Host_code significant?
manysp.2LR2 <- glmer(Reproductive_nodes ~ Euphrasia_sp2 + (1 | Host_given_Euphrasia) + (1 | Population) + (1 | Obs),
                     family = "poisson", data = rnodes3)
# is Host_given_Euphrasia significant?
manysp.2LR3 <- glmer(Reproductive_nodes ~ Euphrasia_sp2 + (1 | Host_code) + (1 | Population) + (1 | Obs),
                     family = "poisson", data = rnodes3)
# is Population significant?
manysp.2LR4 <- glmer(Reproductive_nodes ~ Euphrasia_sp2 + (1 | Host_code) +  (1 | Host_given_Euphrasia) + (1 | Obs),
                     family = "poisson", data = rnodes3)

write.csv(
  x = data.table(
    # Host_code is significant
    HOST_CODE = anova(manysp.2LR1, manysp.2LR2),
    # Host_given_code_Euphrasia is significant
    HOST_GIVEN_EUPHRASIA = anova(manysp.2LR1, manysp.2LR3),
    # Population is highly significant
    HOST_GIVEN_POPULATION = anova(manysp.2LR1, manysp.2LR4), keep.rownames = TRUE
  ), file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/LRTs_of_models.csv"
)




##### Plot 2: Posterior Modes of the interaction model #####

plot_2.1 <- Solapply(manysp.2, coda::HPDinterval)[order(`Posterior Mode`)][Group %in% c("Host_code", "Host_given_Euphrasia")] %>%
  ggplot(aes(x = reorder(Variable, `Posterior Mode`), y = `Posterior Mode`))+
  geom_point()+
  geom_errorbar(aes(ymin = lowerHPD, ymax = upperHPD))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(filename = "./Figures/Many_species/posterior_interaction_modes_mod_1", plot = plot_2.1, 
       device = "pdf", width = 21, height = 6, units = "in")

m1 <- data.table(Density =(manysp.2$VCV[,c(1)])/rowSums(manysp.2$VCV),
                 Trait = as.factor(rep("Host contribution", 1000)))
m2 <- data.table(Density =(manysp.2$VCV[,c(2)])/rowSums(manysp.2$VCV),
                 Trait = as.factor(rep("Population", 1000)))
m3 <- data.table(Density =(manysp.2$VCV[,c(3)])/rowSums(manysp.2$VCV),
                 Trait = as.factor(rep("Host:Euphrasia", 1000)))
m4 <- data.table(Density =(manysp.2$VCV[,c(4)])/rowSums(manysp.2$VCV),
                 Trait = as.factor(rep("Residual", 1000)))
final <- rbind(m1,m2,m3,m4)
final$Density <- as.numeric(final$Density)

final$Trait <- factor(x = final$Trait, levels = rev(c("Host:Euphrasia",
                                                  "Population",
                                                  "Host contribution",
                                                  "Residual")))

plot_2.2<-ggplot(final, aes(x = Density, y = Trait))+
  theme_bw()+ggridges::geom_density_ridges()+geom_vline(xintercept = 0, col = "red", lty = 2, size = 1)+
  scale_x_continuous(limits = c(-0.2,1))+
  theme(strip.text.x = element_text(size=20),
        strip.background = element_rect(colour="white", fill="white"),
        axis.line.x = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x.bottom = element_text(size = 20),
        axis.title.y.left = element_text(size = 20),
        legend.title = element_text(size = 20))

ggsave(filename = "./Figures/Many_species/posterior_interaction_dist", plot = plot_2.2, 
       device = "pdf", width = 10, height = 6, units = "in")

##### Plot 3: Raw data for the manuscript, means and SE's #####

plot_3.1<- rnodes3[, .(mean = mean(Reproductive_nodes, na.rm = TRUE),
            sem = sd(Reproductive_nodes, na.rm = TRUE)/sqrt(.N),
            N = .N), by = c("Euphrasia_sp2", "Host_code","Population2")] %>% #[Population != "M1767"]
  
  ggplot(aes(x = reorder(Host_code, mean) , y = mean))+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem, group=Euphrasia_sp2), position = position_dodge(width = 0.9), width=0.4)+
  geom_point(aes(colour = Euphrasia_sp2), position = position_dodge(width = 0.9), size=3)+
  facet_wrap(~Population2, scales = "free_y")+
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
       device = "png", width = 10, height = 6, units = "in")

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





### REMOVE THESE ANALYSES BELOW ###
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

##### Part 6 onwards deprecated, to be deleted #####
##### Part 6: Host association and interactions #####

# import co-occurrence data

cooccurtet <- fread("./Data/Many_species/Cooccur/Cooccurrence_dataTet.csv")
cooccurvig <- fread("./Data/Many_species/Cooccur/Cooccurrence_dataVig.csv")
cooccurang <- fread("./Data/Many_species/Cooccur/Cooccurrence_dataAng.csv")
cooccurmic <- fread("./Data/Many_species/Cooccur/Cooccurrence_dataMic.csv")

mymerge <- function(x,y) merge(x,y,all=TRUE)

cooccur <- Reduce(mymerge,list(cooccurtet,cooccurvig,cooccurang,cooccurmic))

cooccur[is.na(cooccur)] <- 0

write.csv(x = cooccur, file = "Data/Many_species/Cooccur/Cooccur_all.csv")

cooccurPCA <- prcomp(x = t(cooccur[!"Nodes",-"Species"]), scale. = TRUE, center = TRUE)
cooccurPCA2 <- cooccurPCA$x
cooccurPCA2 <- as.data.frame(cooccurPCA2)
cooccurPCA2$Species <- rownames(cooccurPCA2)

cooccurPCA2$Species <- gsub(pattern = ".[0123456789]+", replacement = "", x = cooccurPCA2$Species)

ggplot(cooccurPCA2, aes(x = PC1, y = PC2))+
  geom_point(aes(colour = Species))+
  theme_bw()

# for adding to model.
# tetraquetra
tet <- data.table(Host = cooccurtet[-1,"Species"], 
                  Rel_abundance = apply(X = cooccurtet[-1,-1], MARGIN = 1, FUN = function(x) sum(x)/length(x)),
                  Diversity = mean(apply(X = cooccurtet[-1,-1], MARGIN = 2, FUN = function(x) sum(x))),
                  Euphrasia_sp2 = "Euphrasia tetraquetra")
vig <- data.table(Host = cooccurvig[-1,"Species"], 
                  Rel_abundance = apply(X = cooccurvig[-1,-1], MARGIN = 1, FUN = function(x) sum(x)/length(x)),
                  Diversity = mean(apply(X = cooccurvig[-1,-1], MARGIN = 1, FUN = function(x) sum(x)), na.rm = TRUE),
                  Euphrasia_sp2 = "Euphrasia vigursii")
ang <- data.table(Host = cooccurang[-1,"Species"],
                  Rel_abundance = apply(X = cooccurang[-1,-1], MARGIN = 1, FUN = function(x) sum(x)/length(x)),
                  Diversity = mean(apply(X = cooccurang[-1,-1], MARGIN = 1, FUN = function(x) sum(x))),
                  Euphrasia_sp2 = "Euphrasia anglica")
mic <- data.table(Host = cooccurmic[-1,"Species"],
                  Rel_abundance = apply(X = cooccurmic[-1,-1], MARGIN = 1, FUN = function(x) sum(x)/length(x)),
                  Diversity = mean(apply(X = cooccurmic[-1,-1], MARGIN = 1, FUN = function(x) sum(x))),
                  Euphrasia_sp2 = "Euphrasia micrantha")

relabund <- rbind(tet,vig,ang,mic)

# add in host species names properly
mnames <-data.table(Host.Species = c("Holcus_lanatus", "Plantago_lanceolata", "Agrostis_curtisii",
                                     "Lolium_perenne", "Festuca_ovina", "Lotus_corniculatus", "Plantago_maritima",
                                     "Hypericum_pulchrum", "Origanum_vulgare", "Ulex_gallii", "Deschampsia_flexuosa",
                                     "Veronica_chamaedrys", "Calluna_vulgaris"),
                    Host_code = unique(rnodes3$Host_code))
# names merged
rnodes4 <- rnodes3[mnames, on= "Host_code"]
# relative abundance merged
rnodes5 <- relabund[rnodes4, on = c("Host.Species", "Euphrasia_sp2")]
# if relative abundance NA, replace with zero
rnodes5$Rel_abundance[is.na(rnodes5$Rel_abundance)] <- 0



prior.manysp.6 <- list(R=list(V=diag(1), nu=0.002), 
                     G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                            #G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                            G3=list(V=diag(2), nu=1, alpha.mu=rep(0,2), alpha.V=diag(2)*1000)))


manysp.6 <-MCMCglmm(Reproductive_nodes ~ Rel_abundance + Euphrasia_sp2,
                   random = ~ Host_code + us(1 + Rel_abundance):Host_given_population,
                   family = "poisson",
                   prior=prior.manysp.6,
                   data = rnodes5[Population2 %!in% c("M1768", "M1769")],
                   nitt = 13000*8,
                   burnin = 3000*8,
                   thin = 10*8,
                   pr=TRUE,
                   verbose = TRUE,
                   saveX = TRUE,
                   saveZ = TRUE)

summary(manysp.6)

MCMCReppois(mod = manysp.6, y = "Host_code")

MCMCReppois(mod = manysp.6, y = "(Intercept):(Intercept).Host_given_population")
MCMCReppois(mod = manysp.6, y = "Rel_abundance:Rel_abundance.Host_given_population")

# visualise the vcv

VCVdensity(manysp.6)+xlim(0,4)

# visualise the model
manysp.6_vis <-Solapply(manysp.6)

manysp.6_vis$Variable <- c(Solapply(manysp.6)$Variable[1:13], gsub(x = Solapply(manysp.6)$Variable, pattern = "^.*?\\.", replacement = "")[14:95])

manysp.6_vis[Group == "(Intercept)"][order(-Grouped_Value)]

# plot the fitted lines
res <- list()
for(i in 1:41){
  res[[i]] <- ablines(intercept = manysp.6_vis$Grouped_Value[manysp.6_vis$Group == "(Intercept)"][i], 
                      slope = manysp.6_vis$Grouped_Value[manysp.6_vis$Group == "Rel_abundance"][i],
                      0,1)
}
names(res) <- manysp.6_vis$Variable[manysp.6_vis$Group == "Rel_abundance"]

res2 <- setDT(as.data.frame(t(data.frame(matrix(unlist(res), nrow=length(res), byrow=T, dimnames = list(names(res)))))))
res2$Rel_abundance <- seq(0,1,0.01)
res3 <- melt(res2, id.vars = "Rel_abundance")

ggplot(rnodes5, aes(x = Rel_abundance, y = (Reproductive_nodes)))+
  geom_jitter()+
  facet_wrap(~Host_code)
  #geom_line(data = res3, aes(x =Rel_abundance, y = value, group = variable))


# but anyway there are only really two large interactions in the slope

manysp.6_vis[Group %in% c("Rel_abundance"),][order(-Grouped_Value)]

##### Plot 6: Alternative visualisations for relative co-occurrence interactions #####

plot6.1 <- rnodes5[, .(mean = mean(Reproductive_nodes),
            sem = sd(Reproductive_nodes)/sqrt(.N),
            N = .N), by = c("Euphrasia_sp2", "Host_code", "Rel_abundance" ,"Population2")] %>% #[Population != "M1767"]
  
  ggplot(aes(x = reorder(Host_code, mean) , y = mean))+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem, group=Euphrasia_sp2), position = position_dodge(width = 0.9), width=0.4)+
  geom_point(aes(alpha = Rel_abundance), position = position_dodge(width = 0.9), size=3)+
  facet_wrap(~Population2, scales = "free_y")+
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

ggsave(filename = "./Figures/Many_species/population_cum_nodes_COOCCUR", plot = plot6.1, 
       device = "pdf", width = 10, height = 6, units = "in")


# next one

# let's loop through this vector to create a new factor
vec <- vector(length = length(rnodes5$Rel_abundance))

for(i in 1:length(rnodes5$Rel_abundance)){
  
  if(rnodes5$Rel_abundance[i] < 0.25){
    vec[i] <- "0-0.25"
  } else 
    if(rnodes5$Rel_abundance[i] >= 0.25 & rnodes5$Rel_abundance[i] < 0.5){
      vec[i] <- "0.25-0.5"
    } else 
      if(rnodes5$Rel_abundance[i] >= 0.5 & rnodes5$Rel_abundance[i] < 0.75){
        vec[i] <- "0.5-0.75"
      } else
        if(rnodes5$Rel_abundance[i] >= 0.75){
          vec[i] <- "0.75-1.0"
        }
}
rnodes5$Rel_abundance2 <- as.factor(vec)

# resulting plot
plot6.2 <- ggplot(rnodes5, aes(x = Host.Species, y = Reproductive_nodes, group = Euphrasia_sp2))+
  geom_point(aes(colour = Rel_abundance2), position = position_dodge(width = 0.75))+
  geom_boxplot(alpha = 0.2, aes(fill = Euphrasia_sp2))+
  facet_grid(~Host.Species, scales = "free_x")+
  theme_bw()

ggsave(filename = "./Figures/Many_species/population_cum_nodes_alternative_COOCCUR", plot = plot6.2, 
       device = "pdf", width = 17, height = 6, units = "in")  

##### Part 7: Analysis of reproductive nodes in the wild #####

cooccur <- read.csv("Data/Many_species/Cooccur/Cooccur_all.csv")
cooccur <- setDT(cooccur[,-1])
colnms <- colnames(cooccur)
tcooccur <- transpose(l = cooccur[,-1])
colnames(tcooccur) <- as.character(cooccur$Species)
tcooccur$Euphrasia <- colnms[-1]

# get rid of the numbers at the end
tcooccur[, Euphrasia := gsub(x = tcooccur$Euphrasia, pattern = ".[0123456789]+", replacement = "")]
# add a space between E. and specific epithet
tcooccur[, Euphrasia := gsub("([A-Z]\\.)([a-z])", "\\1 \\2", tcooccur$Euphrasia)]

tcooccur <- tcooccur[Nodes > 0]

tcooccur2 <- melt.data.table(tcooccur, id.vars = c("Nodes", "Euphrasia"))

# boxplots

ggplot(tcooccur2[value > 0], aes(x = Euphrasia, y = Nodes))+
  geom_jitter()+
  facet_wrap(~variable)+
  geom_boxplot(alpha = 0)

ggplot(tcooccur2[value > 0], aes(x = Euphrasia, y = Nodes))+
  geom_jitter(aes(colour = variable == "Lotus_corniculatus"))+
  geom_boxplot(alpha = 0)

ggplot(tcooccur2[value > 0], aes(x = Euphrasia, y = Nodes))+
  geom_jitter(aes(colour = variable == "Plantago_lanceolata"))+
  geom_boxplot(alpha = 0)


prior.nodes <- list(R=list(V=diag(1), nu=0.002), 
                    G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

nodes <- MCMCglmm(Nodes ~ Euphrasia,
                  random = ~Euphrasia:variable, 
                  data = tcooccur2[value > 0],
                  family = "poisson",
                  nitt = 13000*8,
                  burnin = 3000*8,
                  thin = 10*8,
                  pr=TRUE)

summary(nodes)
VCVdensity(mod = nodes)


# What about diversity?
tcooccur[, .(Diversity = apply(tcooccur[, -c("Euphrasia", "Nodes")], 1, sum),
             Euphrasia = Euphrasia, 
             Nodes = Nodes)] %>%
  ggplot(aes(x = Diversity, y = Nodes))+
  geom_point()+
  geom_smooth(method = "glm", se = FALSE,
              method.args = list(family = "poisson"))
