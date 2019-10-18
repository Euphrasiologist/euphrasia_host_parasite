# Data and analyses for the paper: "Host identity not functional group determines the performance of a parasitic plant"
# Section 1: Many host analyses
# Created: 7.8.18 by Max Brown
# © Max Brown 
# Major update: 11.7.19 # 

# Libraries needed

library(ape)
library(MCMCglmm)
library(data.table)
library(VCVglmm)
library(ggplot2) 
library(ggtree)
library(ggridges)
library(dplyr) 
library(phangorn) 
library(Hmisc)
library(ggrepel)

# functions used

# manipulating logits and probs
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
probit2prob <- function(probit){
  pnorm(probit)
}


# notes 

# Format is: simple stats for each dataset. The model then follows.

# The model is summarised, using summary() and the table is used for the supplementary material.
# Interesting effect sizes and directions can be eyeballed from the model.

# Variance explained by a) host species & phylogeny b) phylogeny (i.e. phylogenetic signal) are estimated.

# Wald tests of fixed effects for significance of functional group and annual/perennial.

# The phylogeny is constrained using the APG4 reference tree and contains 45 taxa.


##### Part 1: Prepare tree for models at first flowering #####

# constraint tree output
constraint <- read.tree(text = "(Agrostis_capillaris:0.0103868960,Lagurus_ovatus:0.0137334981,((((Cynosurus_cristatus:0.0119501076,Festuca_rubra:0.0143927581)56:0.0014297595,Holcus_lanatus:0.0116621244)61:0.0018419113,Phleum_pratense:0.0106048521)100:0.0080082346,((Zea_mays:0.0376747406,((((Allium_ursinum:0.0696170091,Galanthus_nivalis:0.0548052238)amaryllidaceae/100:0.0011547212,Hyacinthoides_non_scripta:0.0509074406)100:0.0322710961,Dactylorhiza_purpurella:0.0989415971)asparagales/100:0.0128940050,((((((Sorbus_aucuparia:0.0314958348,Fragaria_vesca:0.0582563355)rosaceae/100:0.0349788544,((((Ononis_spinosa:0.0220743766,(Lathyrus_japonicus:0.0294240865,Trifolium_pratense:0.0339786483)100:0.0021655535)100:0.0183607112,Lotus_corniculatus:0.0850694608)100:0.0000019136,Ulex_europaeus:0.0718966448)100:0.0000027395,Vicia_cracca:0.0479406736)fabaceae/100:0.0785149642)100:0.0148394884,(Arabidopsis_thaliana:0.1250612607,Helianthemum_nummularium:0.1469119452)malvales_to_brassicales/100:0.0195592609)100:0.0150299043,((((((Meum_athamanticum:0.0221268430,Anthriscus_sylvestris:0.0221983048)apiaceae/100:0.0502387307,Centranthus_ruber:0.1162127616)100:0.0022623680,((Centaurea_nigra:0.0172220154,(Leucanthemum_vulgare:0.0331192970,Tragopogon_pratensis:0.0182447955)100:0.0000023339)100:0.0048658069,Senecio_vulgaris:0.0209331617)asteraceae/100:0.0535848492)100:0.0118251564,(((Thymus_polytrichus:0.0579324847,Mimulus_guttatus:0.0360518925)100:0.0082342727,Plantago_lanceolata:0.0899883709)100:0.0336066248,(Galium_aparine:0.0100601356,Galium_verum:0.0102328988)galium/100:0.1193345682)100:0.0185788636)100:0.0088343295,Erica_tetralix:0.1160061555)ericales_to_asterales/100:0.0079753034,(((Chenopodium_album:0.0332635734,Chenopodium_bonus_henricus:0.0105604686)chenopodium/100:0.0238607464,(Silene_dioica:0.0014311631,Silene_latifolia:0.0000024274)silene/100:0.0808030168)100:0.0364786972,Rumex_acetosella:0.1210989875)caryophyllales/100:0.0287067285)100:0.0088664052)100:0.0236664505,Papaver_rhoeas:0.1015330721)eudicots/100:0.0318216294,(Pinus_sylvestris:0.2367123693,(Equisetum_arvense:0.3878874112,(Cystopteris_fragilis:0.1233674528,Pteridium_aquilinum:0.1016667919)100:1.0445751830)100:0.2401868160)seedplants/100:0.1314408258)poales_to_asterales/100:0.0454413537)100:0.1494656671)100:0.0242458362,Hordeum_vulgare:0.0237302454)93:0.0060960686)poaceae/100:0.0096602096);")

# remove the underscores
constraint$tip.label <- gsub("_", " ", constraint$tip.label) 
# correct the tip labels
constraint$tip.label[10] <- "Hyacinthoides non-scripta"
constraint$tip.label[36] <- "Chenopodium bonus-henricus"
constraint$tip.label[43] <- "Cystopteris dickieana"

# Remove tips that do not occur in the data.
constraint.1<-drop.tip(constraint, tip = c("Dactylorhiza purpurella", "Thymus polytrichus", "Pteridium aquilinum"))
constraint.1$node.label <- NULL
# root tree @ Cystopteris.
constraint.1 <- root(phy = constraint.1, outgroup = "Cystopteris dickieana", resolve.root = TRUE)
# edge cannot be zero, so make it tiny
constraint.1$edge.length[1] <- 1e-10

# comparison of chronoMPL, ultrametric and the difference it makes
# the second is the more realistic and not different from 
# the mean(diag(vcv(tree)))*VCV`animal` method

# create oject for MCMCglmm model
AinvULT <- ultm(constraint.1, MPL = FALSE, inverseA = TRUE)



##### Part 2: Prepare tree for models at end of season growth #####

# load in tree
constraint.2 <- read.tree(text = "(Agrostis_capillaris:0.0103868960,Lagurus_ovatus:0.0137334981,((((Cynosurus_cristatus:0.0119501076,Festuca_rubra:0.0143927581)56:0.0014297595,Holcus_lanatus:0.0116621244)61:0.0018419113,Phleum_pratense:0.0106048521)100:0.0080082346,((Zea_mays:0.0376747406,((((Allium_ursinum:0.0696170091,Galanthus_nivalis:0.0548052238)amaryllidaceae/100:0.0011547212,Hyacinthoides_non_scripta:0.0509074406)100:0.0322710961,Dactylorhiza_purpurella:0.0989415971)asparagales/100:0.0128940050,((((((Sorbus_aucuparia:0.0314958348,Fragaria_vesca:0.0582563355)rosaceae/100:0.0349788544,((((Ononis_spinosa:0.0220743766,(Lathyrus_japonicus:0.0294240865,Trifolium_pratense:0.0339786483)100:0.0021655535)100:0.0183607112,Lotus_corniculatus:0.0850694608)100:0.0000019136,Ulex_europaeus:0.0718966448)100:0.0000027395,Vicia_cracca:0.0479406736)fabaceae/100:0.0785149642)100:0.0148394884,(Arabidopsis_thaliana:0.1250612607,Helianthemum_nummularium:0.1469119452)malvales_to_brassicales/100:0.0195592609)100:0.0150299043,((((((Meum_athamanticum:0.0221268430,Anthriscus_sylvestris:0.0221983048)apiaceae/100:0.0502387307,Centranthus_ruber:0.1162127616)100:0.0022623680,((Centaurea_nigra:0.0172220154,(Leucanthemum_vulgare:0.0331192970,Tragopogon_pratensis:0.0182447955)100:0.0000023339)100:0.0048658069,Senecio_vulgaris:0.0209331617)asteraceae/100:0.0535848492)100:0.0118251564,(((Thymus_polytrichus:0.0579324847,Mimulus_guttatus:0.0360518925)100:0.0082342727,Plantago_lanceolata:0.0899883709)100:0.0336066248,(Galium_aparine:0.0100601356,Galium_verum:0.0102328988)galium/100:0.1193345682)100:0.0185788636)100:0.0088343295,Erica_tetralix:0.1160061555)ericales_to_asterales/100:0.0079753034,(((Chenopodium_album:0.0332635734,Chenopodium_bonus_henricus:0.0105604686)chenopodium/100:0.0238607464,(Silene_dioica:0.0014311631,Silene_latifolia:0.0000024274)silene/100:0.0808030168)100:0.0364786972,Rumex_acetosella:0.1210989875)caryophyllales/100:0.0287067285)100:0.0088664052)100:0.0236664505,Papaver_rhoeas:0.1015330721)eudicots/100:0.0318216294,(Pinus_sylvestris:0.2367123693,(Equisetum_arvense:0.3878874112,(Cystopteris_fragilis:0.1233674528,Pteridium_aquilinum:0.1016667919)100:1.0445751830)100:0.2401868160)seedplants/100:0.1314408258)poales_to_asterales/100:0.0454413537)100:0.1494656671)100:0.0242458362,Hordeum_vulgare:0.0237302454)93:0.0060960686)poaceae/100:0.0096602096);")
# remove underscores
constraint.2$tip.label <- gsub("_", " ", constraint.2$tip.label) 
# correct tip label names
constraint.2$tip.label[10] <- "Hyacinthoides non-scripta"
constraint.2$tip.label[36] <- "Chenopodium bonus-henricus"
constraint.2$tip.label[43] <- "Cystopteris dickieana"
constraint.2$node.label <- NULL


# root tree
constraint.2 <- root(phy = constraint.2, outgroup = "Pteridium aquilinum", resolve.root = TRUE)
# make sure there are no zero values
constraint.2$edge.length[1] <- 1e-10

# create the object for MCMCglmm model 
AinvULT2 <- ultm(constraint.2, MPL = FALSE, inverseA = TRUE)

##### Part 3: Models at first flowering #####

# load data
floweringsofar<- fread("./Data/Many_hosts/FloweringsofarV2.csv")
# Add underscores to spaces
colnames(floweringsofar) <- gsub(pattern = " ", replacement = "_", x = colnames(floweringsofar))

# some data preparation first
# make columns factors
names_factors <- c("AnnPer", "Functional_group", "CSR", "Family", "Host_Species", "ID_No", "Unique_ID")
for (col in names_factors) set(floweringsofar, j=col, value=as.factor(floweringsofar[[col]]))
# days to flower and corolla length change to numeric
names_numeric <- c("Days_since_germination")
for (col in names_numeric) set(floweringsofar, j=col, value=as.numeric(floweringsofar[[col]]))

# relevel so perennial is the baseline
floweringsofar$AnnPer <- relevel(floweringsofar$AnnPer, ref = "Per")
# relevel so grass is baseline
floweringsofar$Functional_group <- relevel(floweringsofar$Functional_group, ref = "Grass")

# change name of Cystopteris dickieana
levels(floweringsofar$Host_Species)[10] <- "Cystopteris dickieana"
# remove no host
floweringsofar <- floweringsofar[Host_Species != "No host"]
floweringsofar$Host_Species <- factor(floweringsofar$Host_Species)

# create 'animal' the phylogenetic random effect variable
floweringsofar$animal <- floweringsofar$Host_Species


## DAYS TO FLOWER ##

prior.2 <- list(R=list(V=diag(1), nu=0.002), 
                G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                       G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

# fixed effects estimates for both annual/perennial and functional group
# the random effects specification allows different variances for each functional group
# but no covariances

# add transplant date?
floweringsofar2 <- allgrowth3[,.(Unique_ID = as.factor(Unique_ID), Transplant.Date)][floweringsofar, on = "Unique_ID"]
floweringsofar2$Transplant.Date <- floweringsofar2$Transplant.Date

# fold difference across all Euphrasia plants in days to flower
max(floweringsofar2$Days_since_germination, na.rm = TRUE)/min(floweringsofar2$Days_since_germination, na.rm = TRUE)

mcmcfix2<- MCMCglmm(Days_since_germination ~ AnnPer + Functional_group + Transplant.Date, 
                    random = ~ Host_Species + animal, 
                    ginverse = list(animal=AinvULT),
                    data = floweringsofar2,
                    prior = prior.2,
                    family= "poisson",
                    nitt = 13000*10,
                    burnin = 3000*10,
                    thin = 10*10,
                    pr=TRUE,
                    verbose = TRUE)
# summary output
summary(mcmcfix2)

mcmc.fix2.plot <- function(){
  x <- seq(0,68,1)
  y <- exp(summary(mcmcfix2)$solutions[1,1] + summary(mcmcfix3)$solutions[7,1]*x)
  plot(x = x, y = y, type = "l", xlab = "", ylab = "")
}
mcmc.fix2.plot()

write.csv(x = specify_decimal(summary(mcmcfix2)$solutions, 4),
          file = "./Data/Many_hosts/Model_outputs/Days_to_flower/Days_solutions.csv")

# wald test of annual/perennial
write.csv(x = aod::wald.test(cov(mcmcfix2$Sol[,2, drop=F]), colMeans(mcmcfix2$Sol[,2, drop=F]), Terms=1)$result,
          file = "./Data/Many_hosts/Model_outputs/Days_to_flower/AnnPer_Wald_Test.csv")
# wald test of functional group
write.csv(x = aod::wald.test(cov(mcmcfix2$Sol[,3:6, drop=F]), colMeans(mcmcfix2$Sol[,3:6, drop=F]), Terms=1:4)$result,
          file = "./Data/Many_hosts/Model_outputs/Days_to_flower/Functional_Group_Wald_Test.csv")

# joint phylogenetic distribution
write.csv(list(posterior.mode(as.mcmc(rowSums(mcmcfix2$VCV[,c(1,2)])/rowSums(mcmcfix2$VCV))),
               HPDinterval(as.mcmc(rowSums(mcmcfix2$VCV[,c(1,2)])/rowSums(mcmcfix2$VCV)))), 
          file = "./Data/Many_hosts/Model_outputs/Days_to_flower/Joint_Phylogeny_Variance.csv")

# and the phylogenetic component

write.csv(list(posterior.mode(as.mcmc((mcmcfix2$VCV[,c(2)])/rowSums(mcmcfix2$VCV[,c(1,2)]))),
               HPDinterval(as.mcmc((mcmcfix2$VCV[,c(2)])/rowSums(mcmcfix2$VCV[,c(1,2)])))),
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Days_to_flower/Phylogeny_Variance.csv")



##### Part 4: Survival analysis (Event History Analysis) ######

# May == 1, June == 2, July == 3, August == 4, September == 5
survivaldata<- fread("./Data/Many_hosts/Survivalanalysis3.csv")
survivaldata
# add time as a column in the data, i.e. at which time points was each individual Euphrasia alive
survivaldata <- survivaldata[, Time := .(Time = 1:.N), by="Unique_ID"][Time < 6,]
# a table of individuals alive (0) or dead (1)
table(survivaldata$Time, survivaldata$y)

# Average proability of death on each host (Is this average probability??) (amended...)
survivaldata[, .(Mean = mean((y-1)*-1 )), by = c("Name", "Time")][,.(Mean = mean(Mean)), by = "Name"][order(-Mean)]
survivaldata[Name == "Zea mays", .(mean((y-1)*-1 )), by = "Time"]


# numbers surviving at each time point
survivaldata[, .(N = .N), by = c("Time")]
# Average probability of death at eg time point 3
survivaldata[, .(Mean = mean(y)), by = c("Name", "Time")][Time == 3,]
# 506 plants alive at time 4
table(survivaldata[Time == 4]$y)
# which ones were they?
table(survivaldata[Time == 4]$Name, survivaldata[Time == 4]$y)
# equivalent but more easily manipulated
# Galium aparine and Plantago lanceolata had most survival, 28 individuals on each
survivaldata[Time == 3, .("0" = sum(as.integer(y == 0)),
                          "1" = sum(as.integer(y != 0))) ,by = "Name"][order(-`1`)]


# some data tidying
tmp <- as.POSIXlt(survivaldata$Latest.Germ.Date, format = "%d/%m/%y")
survivaldata$Latest.Germ.Date<-tmp$yday

tmp <- as.POSIXlt(survivaldata$Transplant.Date, format = "%d/%m/%y")
survivaldata$Transplant.Date<-tmp$yday

# model 
survivaldata[, animal := as.factor(Name)]
levels(survivaldata$animal)[10] <- "Cystopteris dickieana"

# add annual perennial, functional group, and Nodes
allgrowth <- fread("./Data/Many_hosts/Allgrowthmeasurements1.csv")
colnames(allgrowth) <- gsub(pattern = " ", replacement = "_", x = colnames(allgrowth))

# select the columns needed
to_add <- allgrowth[, c("Unique_ID", "AnnPer", "Functional_group", "C.R.Nodes.F", "Family")]

# merge datasets
survivaldata.2 <- survivaldata[to_add, on = "Unique_ID"][Name != "No host"][Time < 5]
# change name, time, annper and functional group, animal to factor
names_factors2 <- c("Name", "Time", "AnnPer", "Functional_group", "animal")
for (col in names_factors2) set(survivaldata.2, j=col, value=as.factor(survivaldata.2[[col]]))

# make sure the levels are levelled rightly
survivaldata.2$AnnPer <- relevel(survivaldata.2$AnnPer, ref = "Per")
survivaldata.2$Functional_group <- relevel(survivaldata.2$Functional_group, ref = "Grass")

# probability of survival
survivaldata.2$y <- (survivaldata.2$y -1)*-1 
# make transplant date start at day 97 (=0)
survivaldata.2$Transplant.Date <- survivaldata.2$Transplant.Date -97

# prior for the event history analysis
prior.eha <-list(R=list(V=diag(1), fix=1), 
                 G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                        G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))


eha.1 <- MCMCglmm(y ~ 1 + Time+AnnPer+ Transplant.Date + Functional_group,
                  random = ~Name + animal,
                  ginverse = list(animal=AinvULT2),
                  data = survivaldata.2,
                  prior=prior.eha,
                  family = "threshold",
                  nitt = 13000*10,
                  burnin = 3000*10,
                  thin = 10*10,
                  verbose = TRUE,
                  pr=TRUE)

summary(eha.1)
# print table to put into the manuscript:
write.csv(x = specify_decimal(summary(eha.1)$solutions, 4), 
          file = "./Data/Many_hosts/Model_outputs/Survival/EHA_solutions.csv")

# functional group wald test
write.csv(x = aod::wald.test(cov(eha.1$Sol[,7:10, drop=F]), colMeans(eha.1$Sol[,7:10, drop=F]), Terms=1:4)$result$chi2,
          file = "./Data/Many_hosts/Model_outputs/Survival/Functional_group_Wald_Test.csv")
# annual perennial
write.csv(x = aod::wald.test(cov(eha.1$Sol[,5, drop=F]), colMeans(eha.1$Sol[,5, drop=F]), Terms=1)$result$chi2,
          file = "./Data/Many_hosts/Model_outputs/Survival/AnnPer_Wald_Test.csv")
# covariance of functional group. Save
write.csv(x = cov(eha.1$Sol[,7:10, drop=F]), 
          file = "./Data/Many_hosts/Model_outputs/Survival/Var_Covar_Functional_group.csv")

# joint phylogenetic distribution of variance
write.csv(x = list(posterior.mode(as.mcmc(rowSums(eha.1$VCV[,c(1,2)])/rowSums(eha.1$VCV))),
                   HPDinterval(as.mcmc(rowSums(eha.1$VCV[,c(1,2)])/rowSums(eha.1$VCV)))),
          file = "./Data/Many_hosts/Model_outputs/Survival/Joint_Phylogeny_Variance.csv")

# phylogenetic signal
write.csv(x = list(posterior.mode(eha.1$VCV[,c(2)]/rowSums(cbind(eha.1$VCV[,c(2)], 
                                                                 eha.1$VCV[,c(1)]))),
                   HPDinterval(eha.1$VCV[,c(2)]/rowSums(cbind(eha.1$VCV[,c(2)], 
                                                              eha.1$VCV[,c(1)])))),
          file = "./Data/Many_hosts/Model_outputs/Survival/Phylogenetic_Signal.csv")

# Visualise the fixed effects
MCMCfixplot(eha.1)
# Visualise the variances
VCVdensity(eha.1)

##### Part 5: End of season reproductive nodes #####

allgrowth<- read.csv("/Users/mbrown/OneDrive - University of Edinburgh/Euphrasia Experiment 1 HOSTS/GROWTH EXPT DATA/Allgrowthmeasurements1.csv", na.strings = "-")
# need to merge to get transplant date
allgrowth2<-allgrowth[, c("Unique_ID", "ID_No", "AnnPer", "Functional_group", "CSR", "Family", "Name", "C.R.Nodes.F")]
allgrowth3 <- allgrowth2[unique(survivaldata.2[,c("Unique_ID", "Transplant.Date")]), on = "Unique_ID"][Name != "No host"]
# make sure data is of right type
names_factors3 <- c("Name", "AnnPer", "Functional_group", "animal")
for (col in names_factors3) set(allgrowth3, j=col, value=as.factor(allgrowth3[[col]]))
# check tree and animal match
setdiff(constraint.2$tip.label, unique(allgrowth3$Name))
# there's always one...
levels(allgrowth3$Name)[10] <- "Cystopteris dickieana"

# add animal for phylogenetic effects
allgrowth3[, animal := Name]

# some summary stats

allgrowth3[, .(Mean.Nodes = mean(C.R.Nodes.F)), by = "Name"][order(-Mean.Nodes)]
# adjust transplant date
allgrowth3$Transplant.Date <- allgrowth3$Transplant.Date -97

# relevel factors
allgrowth3$AnnPer <- relevel(x = allgrowth3$AnnPer, ref = "Per")
allgrowth3$Functional_group <- relevel(x = allgrowth3$Functional_group, ref = "Grass")


# model #
prior.3<-list(R=list(V=diag(1), nu=0.002),
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                     G2=list(V=diag(1), n=1, alpha.mu=rep(0, 1),alpha.V=diag(1)*1000)))

mcmcfix3<-MCMCglmm(C.R.Nodes.F ~ AnnPer + Functional_group + Transplant.Date,
                   random = ~Name + animal,
                   ginverse = list(animal=AinvULT2),
                   family="poisson", 
                   prior=prior.3, 
                   data=allgrowth3,
                   nitt=13000*20,
                   thin=10*20,
                   burnin=3000*20,
                   verbose = TRUE,
                   pr=TRUE)
summary(mcmcfix3)

mcmc.fix3.plot <- function(){
  x <- seq(0,68,1)
  y <- exp(summary(mcmcfix3)$solutions[1,1] + summary(mcmcfix3)$solutions[4,1] + summary(mcmcfix3)$solutions[7,1]*x)
  plot(x = x, y = y, type = "l", xlab = "", ylab = "")
  y2 <- exp(summary(mcmcfix3)$solutions[1,1] + summary(mcmcfix3)$solutions[6,1] + summary(mcmcfix3)$solutions[7,1]*x)
  lines(x = x, y = y2, type = "l")
}
mcmc.fix3.plot()

# AnnPer
write.csv(x = aod::wald.test(cov(mcmcfix3$Sol[,2, drop=F]), colMeans(mcmcfix3$Sol[,2, drop=F]), Terms=1)$result$chi2,
          file = "./Data/Many_hosts/Model_outputs/Reproductive_nodes_end/AnnPer_Wald_Test.csv")
# Functional group
write.csv(x = aod::wald.test(cov(mcmcfix3$Sol[,3:6, drop=F]), colMeans(mcmcfix3$Sol[,3:6, drop=F]), Terms=1:4)$result$chi2,
          file = "./Data/Many_hosts/Model_outputs/Reproductive_nodes_end/Functional_Group_Wald_Test.csv")

# joint phylogenetic distribution
write.csv(list(posterior.mode(as.mcmc(rowSums(mcmcfix3$VCV[,c(1,2)])/rowSums(mcmcfix3$VCV))),
               HPDinterval(as.mcmc(rowSums(mcmcfix3$VCV[,c(1,2)])/rowSums(mcmcfix3$VCV)))),
          file = "./Data/Many_hosts/Model_outputs/Reproductive_nodes_end/Joint_Phylogenetic_Variance.csv")

# and the phylogenetic component
# 74.5%
write.csv(list(posterior.mode(mcmcfix3$VCV[,c(2)]/rowSums(cbind(mcmcfix3$VCV[,c(2)], 
                                                                mcmcfix3$VCV[,c(1)]))),
               # 14.8% - 96.1%
               HPDinterval(mcmcfix3$VCV[,c(2)]/rowSums(cbind(mcmcfix3$VCV[,c(2)], 
                                                             mcmcfix3$VCV[,c(1)]))) ),
          file = "./Data/Many_hosts/Model_outputs/Reproductive_nodes_end/Phylogenetic_Variance.csv")

##### Part 6: Reproductive nodes over time #####

allgrowth4 <- allgrowth[, .(Unique_ID, ID_No, AnnPer, Functional_group, CSR, Family, Name, R.Nodes.1, R.Nodes.2, R.Nodes.3, R.Nodes.4, R.Nodes.5)]

allgrowth5 <- melt.data.table(data = allgrowth4, id.vars = c("Unique_ID", "ID_No", "AnnPer", "Functional_group", "CSR", "Family", "Name"), 
                variable.name = "Time", value.name = "Nodes")
allgrowth6 <- allgrowth5[, Time := gsub(pattern = "R.Nodes.", replacement = "", x = Time)][Time %in% c("2","3", "4")][Name != "No host"]

names_factors4 <- c("Name", "AnnPer", "Functional_group", "animal", "Time")
for (col in names_factors4) set(allgrowth6, j=col, value=as.factor(allgrowth6[[col]]))

# check tip labels
setdiff(unique(allgrowth6$Name), constraint.2$tip.label)
levels(allgrowth6$Name)[10] <- "Cystopteris dickieana"

allgrowth6$AnnPer <- relevel(x = allgrowth6$AnnPer, ref = "Per")
allgrowth6$Functional_group <- relevel(x = allgrowth6$Functional_group, ref = "Grass")

# add animal
allgrowth6[, animal := Name]
# add transplant date
allgrowth7 <- allgrowth3[,.(Unique_ID = as.integer(Unique_ID), Transplant.Date)][allgrowth6, on = "Unique_ID"]

prior.4<-list(R=list(V=diag(3), nu=0.002), 
              G=list(G1=list(V=diag(3), nu=3, alpha.mu=rep(0,3), alpha.V=diag(3)*1000),
                     G2= list(V=diag(1), nu=3, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))
# run for 13,000 X 20 iterations...
mcmcfix4<-MCMCglmm(Nodes~Time*AnnPer+Functional_group + Transplant.Date,
                   random = ~ us(Time):Name + animal,
                   ginverse = list(animal=AinvULT2),
                   rcov=~idh(Time):units,
                   family="poisson", 
                   prior=prior.4, 
                   data=allgrowth7,
                   nitt=13000*20,
                   thin=10*20,
                   burnin=3000*20,
                   verbose = TRUE,
                   pr=TRUE)
summary(mcmcfix4)
# summary output
write.csv(x = summary(mcmcfix4)$solutions, 
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Reproduction_over_time/OT_Solutions.csv")
# functional group
write.csv(x = aod::wald.test(cov(mcmcfix4$Sol[,5:8, drop=F]), colMeans(mcmcfix4$Sol[,5:8, drop=F]), Terms=1:4)$result$chi2,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Reproduction_over_time/Functional_Group_Wald_Test.csv")
# life history host: time point
write.csv(x = aod::wald.test(cov(mcmcfix4$Sol[,9:10, drop=F]), colMeans(mcmcfix4$Sol[,9:10, drop=F]), Terms=1:2)$result$chi2,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Reproduction_over_time/AnnPerInt_Wald_Test.csv")
# life history of host
write.csv(x = aod::wald.test(cov(mcmcfix4$Sol[,4, drop=F]), colMeans(mcmcfix4$Sol[,4, drop=F]), Terms=1)$result$chi2,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Reproduction_over_time/AnnPer_Wald_Test.csv")

# expected number of nodes at time point 1 on perennial 

exp(3.27173) # = 26 times the expected number of nodes at time point 4 compared to time point 2
exp(0.50560) # 1.6 times the expected number of nodes compared to a perennial grass host at time 2
# interaction at time point four reduces expected number of nodes by 24 times
exp(3.27173+-2.38964)


##### Plot 1: Event history analysis #####

# make time numeric for smooth interpolation
survivaldata.2$Time <- as.numeric(survivaldata.2$Time)

# this plot takes three largest families, two coincide with their own functional group
plot_1 <- ggplot(survivaldata.2[Family %in% c("Poaceae", "Fabaceae")], 
                 aes(x=Time, y=y,group=Name))+
  geom_point(alpha=0.1, position = position_jitter(width = 0.2, height = 0.2))+
  facet_wrap(~Family, scales = "free_x")+
  geom_smooth(aes(x=Time, y=y, col = Name), 
              method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE, size = 1)+
  ylim(c(0,1))+
  theme_bw()+
  ylab(label = "Probability of death")+
  labs(colour = "Species")+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20))


ggsave(filename = "./Figures/Many_hosts/two_families_survival", plot = plot_1, 
       device = "pdf", width = 6, height = 5, units = "in")

plot_1_group <- ggplot(survivaldata.2[Family %in% c("Poaceae", "Fabaceae")], 
                 aes(x=Time, y=y,group=Family))+
  geom_point(alpha=0.1, position = position_jitter(width = 0.2, height = 0.2))+
  #facet_wrap(~Family, scales = "free_x")+
  geom_smooth(aes(x=Time, y=y, col = Family), 
              method = "glm", 
              method.args = list(family = "binomial"), 
              se = TRUE, size = 1)+
  ylim(c(0,1))+
  theme_bw()+
  ylab(label = "Probability of death")+
  labs(colour = "Species")+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20))

ggsave(filename = "./Figures/Many_hosts/two_families_survival_joined", plot = plot_1_group, 
       device = "pdf", width = 6, height = 5, units = "in")

##### Plot 2: Reproductive nodes at end of season with phylogeny #####

# create the ggtree object
(p<-ggtree(phytools::force.ultrametric(constraint.2))+
    geom_tiplab(size = 4)+ theme(axis.line = element_blank())+
    geom_cladelabel(node=64, label="Legumes", align=TRUE, angle = 270, offset =1, hjust = 0.5, offset.text = 0.05, barsize = 1.5)+
    geom_cladelabel(node=53, label="Grasses", align=TRUE, angle = 270, offset =1, hjust = 0.5, offset.text = 0.05, barsize = 1.5))
# data to go with the tree
phylodat<-allgrowth3[,c("Name","C.R.Nodes.F", "Functional_group", "Family")]
phylodat$group <- phylodat$Name

# group the data by host species
phylodat1.1 <- phylodat[, .(Mean.Cum = mean(na.omit(C.R.Nodes.F)),
                            SD.Cum = sd(na.omit(C.R.Nodes.F))/sqrt(.N),
                            N = .N,
                            Functional.group = (Functional_group),
                            group = group,
                            Family = (Family)), by = "Name"]
phylodat1.1 <- unique(phylodat1.1)

# for making the boxplots
(finalplot<-facet_plot(p + xlim_tree(1.1), panel = "Cumulative Reproductive Nodes",
                       data = phylodat1.1,
                       mapping = aes(x=Mean.Cum,
                                     xmin = Mean.Cum-SD.Cum,
                                     xmax = Mean.Cum +SD.Cum,
                                     group = group),
                       geom = geom_errorbarh)+
    theme_tree2())
# overlay the red observational points
plot_2 <- facet_plot(finalplot, panel = "Cumulative Reproductive Nodes",
           data = phylodat1.1,
           mapping = aes(x=Mean.Cum, group = group),
           geom = geom_point, col = "red", size = 3)+
  theme_bw()+
  theme(strip.text.x = element_text(size=20, face = "bold", vjust = 1),
        strip.background = element_rect(colour="white", fill="white"),
        axis.line.x = element_line(colour = "black"),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank())+
  geom_hilight(node=53, fill="darkgreen", alpha=.2, extendto = 1.5)+
  geom_hilight(node=64, fill="steelblue", alpha=.2, extendto = 1.5)+
  geom_hilight(node=31, fill="tomato", alpha=.2, extendto = 1.5)+
  geom_hline(yintercept = c(8.5, 16.5, 21.5, 27.5), lty =2, alpha = 0.3)

ggsave(filename = "./Figures/Many_hosts/Cum_reprod_nodes_phylo", plot = plot_2, 
       device = "pdf", width = 9, height = 7, units = "in")


##### Plot 3: Effect of annual/perennial over time (host life history) #####
# reproductive nodes over time (for functional group and annual or perennial!)
plot_3 <- allgrowth5[, .(Nodes = mean(Nodes), uCI = mean(Nodes)+sd(Nodes)/sqrt(.N), lCI = mean(Nodes)-sd(Nodes)/sqrt(.N)), by = c("AnnPer", "Time")] %>%
  ggplot(aes(x = Time, y= Nodes))+ 
  facet_wrap(~AnnPer)+
  geom_point()+
  geom_errorbar(aes(ymax = uCI, ymin =lCI), width = 0.05)+
  geom_line(aes(group = AnnPer))+
  theme_bw()+
  theme(strip.background = element_blank(), strip.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

ggsave(filename = "./Figures/Many_hosts/ann_per_time", plot = plot_3, 
       device = "pdf", width = 6, height = 5, units = "in")

##### Plot 4: What trajectory does an 'average' host follow over time? ######

# create the average host data
avghost<-allgrowth5[, .(Nodes = mean(Nodes), uCI = mean(Nodes)+sd(Nodes)/sqrt(.N), lCI = mean(Nodes)-sd(Nodes)/sqrt(.N)), by = c("Time")]
avghost[, c("Name", "Family", "Time") := list("Average Host", "Average", 1:5)]

# bind it to performance over time
overtime <- allgrowth5[, .(Family = Family, Nodes = mean(Nodes), uCI = mean(Nodes)+sd(Nodes)/sqrt(.N), lCI = mean(Nodes)-sd(Nodes)/sqrt(.N)), by = c("Time", "Name")]
overtime <- rbind(overtime, avghost)
# get levels for the plot
first11 <- allgrowth3[, .(Mean.Nodes = mean(C.R.Nodes.F)), by = "Name"][order(-Mean.Nodes)]$Name[1:11]
overtime <- overtime[Name %in% c(as.character(first11), "Average Host")]
# set the levels
overtime$Name <- factor(overtime$Name, levels = c(as.character(first11), "Average Host"))
set(overtime, j="Time", value=as.numeric(overtime[["Time"]]))

plot_4 <- overtime %>%
  ggplot(aes(x = Time, y= Nodes, group = Name))+ 
  facet_wrap(~Name)+
  geom_errorbar(aes(ymax = uCI, ymin =lCI), width = 0.05)+
  geom_line(aes(colour = Family))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20))

ggsave(filename = "./Figures/Many_hosts/average_host", plot = plot_4, 
       device = "pdf", width = 12.5, height = 7.5, units = "in")

##### Plot 5: Show the phylogenetic signal posterior distributions #####

  # ridge plot for the phylogenetic constribution to the four main traits. 
  d1 <- data.table(Density =rowSums(mcmcfix2$VCV[,c(1,2)])/rowSums(mcmcfix2$VCV),
                   Trait = as.factor(rep("Days to Flower (Poisson)", 1000)))
  s1 <- data.table(Density =rowSums(eha.1$VCV[,c(1,2)])/rowSums(eha.1$VCV),
                   Trait = as.factor(rep("Survival (Threshold, probit)", 1000)))
  r1 <- data.table(Density =rowSums(mcmcfix3$VCV[,c(1,2)])/rowSums(mcmcfix3$VCV),
                   Trait = as.factor(rep("Total Reproductive Output (Poisson)", 1000)))
  final <- rbind(d1, s1, r1)
  final$Density <- as.numeric(final$Density)
  
  final$Trait <- factor(x = final$Trait, levels = c("Days to Flower (Poisson)",
                                                    "Total Reproductive Output (Poisson)",
                                                    "Survival (Threshold, probit)"))
  
  plot_5 <- ggplot(final, aes(x = Density, y = Trait))+
    theme_bw()+geom_density_ridges()+geom_vline(xintercept = 0, col = "red", lty = 2, size = 1)+
    scale_x_continuous(limits = c(-0.2,1))

ggsave(filename = "./Figures/Many_hosts/joint_variance_explained", plot = plot_5, 
       device = "pdf", width = 6, height = 5, units = "in")

##### DUSTBIN #####

## COROLLA LENGTH ##
prior.1 <- list(R=list(V=diag(1), nu=0.002), 
                G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                       G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

# fixed effects estimates for both annual/perennial and functional group
mcmcfix1 <- MCMCglmm(Corolla_Length ~ AnnPer + Functional_group, 
                     random = ~ Host_Species + animal, 
                     ginverse = list(animal = AinvULT),
                     data = floweringsofar,
                     prior = prior.1,
                     family= "gaussian",
                     nitt = 13000*10,
                     burnin = 3000*10,
                     thin = 10*10,
                     pr=TRUE,
                     verbose = TRUE)


# summary output
summary(mcmcfix1)

write.csv(x = specify_decimal(summary(mcmcfix1)$solutions, 4),
          file = "./Data/Many_hosts/Model_outputs/Corolla_length/Corolla_solutions.csv")

# joint phylogenetic distribution
# 66.5%
write.csv(
  list(posterior.mode(as.mcmc(rowSums(mcmcfix1$VCV[,c(1,2)])/rowSums(mcmcfix1$VCV))),
       # 18.5% - 80.4%
       HPDinterval(as.mcmc(rowSums(mcmcfix1$VCV[,c(1,2)])/rowSums(mcmcfix1$VCV)))), 
  file = "./Data/Many_hosts/Model_outputs/Corolla_length/Joint_Phylogeny_Variance.csv")

# and the phylogenetic component
# 60.2%
write.csv(
  list(posterior.mode(as.mcmc((mcmcfix1$VCV[,c(2)])/rowSums(mcmcfix1$VCV))),
       # 9.41% - 80.6%
       HPDinterval(as.mcmc((mcmcfix1$VCV[,c(2)])/rowSums(mcmcfix1$VCV)))), 
  file = "./Data/Many_hosts/Model_outputs/Corolla_length/Phylogeny_Variance.csv")


# host variance (v. low)
posterior.mode(as.mcmc((mcmcfix1$VCV[,c(1)])/rowSums(mcmcfix1$VCV)))
HPDinterval(as.mcmc((mcmcfix1$VCV[,c(1)])/rowSums(mcmcfix1$VCV)))


# wald test of annual/perennial
write.csv(x = aod::wald.test(cov(mcmcfix1$Sol[,2, drop=F]), colMeans(mcmcfix1$Sol[,2, drop=F]), Terms=1)$result, 
          file = "./Data/Many_hosts/Model_outputs/Corolla_length/AnnPer_Wald_Test.csv")
# wald test of functional group
write.csv(x = aod::wald.test(cov(mcmcfix1$Sol[,3:6, drop=F]), colMeans(mcmcfix1$Sol[,3:6, drop=F]), Terms=1:4)$result, 
          file = "./Data/Many_hosts/Model_outputs/Corolla_length/Functional_Group_Wald_Test.csv")
