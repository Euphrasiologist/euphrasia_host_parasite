# Data and analyses for the paper: "Host identity not functional group determines the performance of a parasitic plant"
# Section 1: Many host analyses
# Created: 7.8.18 by Max Brown
# © Max Brown 
# Major update: 11.7.19 # 

# Libraries needed

library(ape)
library(MCMCglmm)
library(data.table)
#devtools::install_github("Euphrasiologist/VCVglmm")
library(VCVglmm) ## not on CRAN
library(ggplot2) 
library(ggtree)
library(ggridges)
library(dplyr) 
library(phangorn) 
library(Hmisc)
library(ggrepel)
library(ggnewscale)

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

VCVapply <- function(model, FUN = coda::HPDinterval){
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
  
  fun <- match.fun(FUN)
  
  if (identical(fun, HPDinterval)) {
    temp <- setDT(as.data.frame(HPDinterval(model$VCV)), 
                  keep.rownames = TRUE)
    temp$posterior.mode <- as.data.frame(apply(model$VCV, 2, posterior_mode))[, 1]
    colnames(temp) <- c("Variable", "lowerHPD", "upperHPD", 
                        "Posterior Mode")
    temp <- temp[, .(Variable,`Posterior Mode`, lowerHPD, upperHPD)]
    return(temp)
  }
  temp <- setDT(as.data.frame(apply(model$VCV, 2, fun)), keep.rownames = TRUE)
  colnames(temp) <- c("Variable", name <- paste(deparse(substitute(FUN))))
  name
  return(temp)
}

# notes 

# Format is: simple stats for each dataset. The model then follows.

# The model is summarised, using summary() and the table is used for the supplementary material.
# Interesting effect sizes and directions can be eyeballed from the model.

# Variance explained by a) host species & phylogeny b) phylogeny (i.e. phylogenetic signal) are estimated.

# Wald tests of fixed effects for significance of functional group and annual/perennial.

# The phylogeny is constrained using the APG4 reference tree and contains 45 taxa.

# time points are: May == 1, June == 2, July == 3, August == 4, September == 5

##### Part 1: Prepare tree for models at first flowering #####

# constraint tree output
constraint <- read.tree(text = "(Agrostis_capillaris:0.0103868960,Lagurus_ovatus:0.0137334981,((((Cynosurus_cristatus:0.0119501076,Festuca_rubra:0.0143927581)56:0.0014297595,Holcus_lanatus:0.0116621244)61:0.0018419113,Phleum_pratense:0.0106048521)100:0.0080082346,((Zea_mays:0.0376747406,((((Allium_ursinum:0.0696170091,Galanthus_nivalis:0.0548052238)amaryllidaceae/100:0.0011547212,Hyacinthoides_non_scripta:0.0509074406)100:0.0322710961,Dactylorhiza_purpurella:0.0989415971)asparagales/100:0.0128940050,((((((Sorbus_aucuparia:0.0314958348,Fragaria_vesca:0.0582563355)rosaceae/100:0.0349788544,((((Ononis_spinosa:0.0220743766,(Lathyrus_japonicus:0.0294240865,Trifolium_pratense:0.0339786483)100:0.0021655535)100:0.0183607112,Lotus_corniculatus:0.0850694608)100:0.0000019136,Ulex_europaeus:0.0718966448)100:0.0000027395,Vicia_cracca:0.0479406736)fabaceae/100:0.0785149642)100:0.0148394884,(Arabidopsis_thaliana:0.1250612607,Helianthemum_nummularium:0.1469119452)malvales_to_brassicales/100:0.0195592609)100:0.0150299043,((((((Meum_athamanticum:0.0221268430,Anthriscus_sylvestris:0.0221983048)apiaceae/100:0.0502387307,Centranthus_ruber:0.1162127616)100:0.0022623680,((Centaurea_nigra:0.0172220154,(Leucanthemum_vulgare:0.0331192970,Tragopogon_pratensis:0.0182447955)100:0.0000023339)100:0.0048658069,Senecio_vulgaris:0.0209331617)asteraceae/100:0.0535848492)100:0.0118251564,(((Thymus_polytrichus:0.0579324847,Mimulus_guttatus:0.0360518925)100:0.0082342727,Plantago_lanceolata:0.0899883709)100:0.0336066248,(Galium_aparine:0.0100601356,Galium_verum:0.0102328988)galium/100:0.1193345682)100:0.0185788636)100:0.0088343295,Erica_tetralix:0.1160061555)ericales_to_asterales/100:0.0079753034,(((Chenopodium_album:0.0332635734,Chenopodium_bonus_henricus:0.0105604686)chenopodium/100:0.0238607464,(Silene_dioica:0.0014311631,Silene_latifolia:0.0000024274)silene/100:0.0808030168)100:0.0364786972,Rumex_acetosella:0.1210989875)caryophyllales/100:0.0287067285)100:0.0088664052)100:0.0236664505,Papaver_rhoeas:0.1015330721)eudicots/100:0.0318216294,(Pinus_sylvestris:0.2367123693,(Equisetum_arvense:0.3878874112,(Cystopteris_fragilis:0.1233674528,Pteridium_aquilinum:0.1016667919)100:1.0445751830)100:0.2401868160)seedplants/100:0.1314408258)poales_to_asterales/100:0.0454413537)100:0.1494656671)100:0.0242458362,Hordeum_vulgare:0.0237302454)93:0.0060960686)poaceae/100:0.0096602096);")
bootstrap_constraint_tree <- read.tree(text = "(Agrostis_capillaris:0.0107710360,Lagurus_ovatus:0.0141798036,((((Cynosurus_cristatus:0.0123446665,Festuca_rubra:0.0149004779)55:0.0014738722,Holcus_lanatus:0.0120445040)61:0.0018973896,Phleum_pratense:0.0109832005)100:0.0082685225,((Zea_mays:0.0387534621,((((Allium_ursinum:0.0711800196,Galanthus_nivalis:0.0563932560)100:0.0012469692,Hyacinthoides_non_scripta:0.0521906528)100:0.0324660729,Dactylorhiza_purpurella:0.1012756953)100:0.0132311314,((((((Sorbus_aucuparia:0.0320477420,Fragaria_vesca:0.0597722657)100:0.0355510987,((((Ononis_spinosa:0.0225239509,(Lathyrus_japonicus:0.0303051384,Trifolium_pratense:0.0348177570)100:0.0024167048)100:0.0195230095,Lotus_corniculatus:0.0867822377)100:0.0000025677,Ulex_europaeus:0.0727823582)100:0.0000028302,Vicia_cracca:0.0505748224)100:0.0797846620)100:0.0152630847,(Arabidopsis_thaliana:0.1279453101,Helianthemum_nummularium:0.1509497337)100:0.0190292662)100:0.0150228594,((((((Meum_athamanticum:0.0230119844,Anthriscus_sylvestris:0.0227254472)100:0.0512075747,Centranthus_ruber:0.1187558648)100:0.0022990675,((Centaurea_nigra:0.0177791451,(Leucanthemum_vulgare:0.0341583391,Tragopogon_pratensis:0.0189617581)100:0.0000023355)100:0.0050033436,Senecio_vulgaris:0.0214150336)100:0.0547642617)100:0.0122871060,(((Thymus_polytrichus:0.0594431597,Mimulus_guttatus:0.0370743890)100:0.0083251593,Plantago_lanceolata:0.0918113826)100:0.0343244347,(Galium_aparine:0.0103710922,Galium_verum:0.0105322854)100:0.1214619832)100:0.0191267155)100:0.0093641745,Erica_tetralix:0.1184830873)100:0.0079251519,(((Chenopodium_album:0.0343216889,Chenopodium_bonus_henricus:0.0110239352)100:0.0247768734,(Silene_dioica:0.0014916328,Silene_latifolia:0.0000021549)100:0.0826029828)100:0.0373898514,Rumex_acetosella:0.1235274096)100:0.0292642418)100:0.0094711485)100:0.0238596466,Papaver_rhoeas:0.1030491281)100:0.0331534728,(Pinus_sylvestris:0.2391077191,(Equisetum_arvense:0.3946018682,(Cystopteris_fragilis:0.1247748389,Pteridium_aquilinum:0.1041459185)100:1.0852017555)100:0.2437925802)100:0.1327500423)100:0.0457150606)100:0.1511567568)100:0.0248840705,Hordeum_vulgare:0.0244810139)96:0.0062754228)100:0.0099786179);")

# remove the underscores
constraint$tip.label <- gsub("_", " ", constraint$tip.label) 
bootstrap_constraint_tree$tip.label <- gsub("_", " ", bootstrap_constraint_tree$tip.label) 
# correct the tip labels
constraint$tip.label[10] <- "Hyacinthoides non-scripta"
constraint$tip.label[36] <- "Chenopodium bonus-henricus"
constraint$tip.label[43] <- "Cystopteris dickieana"

bootstrap_constraint_tree$tip.label[10] <- "Hyacinthoides non-scripta"
bootstrap_constraint_tree$tip.label[36] <- "Chenopodium bonus-henricus"
bootstrap_constraint_tree$tip.label[43] <- "Cystopteris dickieana"

# Remove tips that do not occur in the data.
constraint.1 <- drop.tip(constraint, tip = c("Dactylorhiza purpurella", "Thymus polytrichus", "Pteridium aquilinum"))
bootstrap_constraint_tree.1 <- drop.tip(bootstrap_constraint_tree, tip = c("Dactylorhiza purpurella", "Thymus polytrichus", "Pteridium aquilinum"))
constraint.1$node.label <- NULL

# do not make node labels null
bootstrap_constraint_tree.1

# root tree @ Cystopteris.
constraint.1 <- root(phy = constraint.1, outgroup = "Cystopteris dickieana", resolve.root = TRUE)
bootstrap_constraint_tree.1 <- root(phy = bootstrap_constraint_tree.1, outgroup = "Cystopteris dickieana", resolve.root = TRUE)
# edge cannot be zero, so make it tiny
constraint.1$edge.length[1] <- 1e-10
bootstrap_constraint_tree.1$edge.length[1] <- 1e-10

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

##### Part 3: Model at first flowering #####

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
# days to flower means and errors
flmeans <- floweringsofar2[, .(mean = as.numeric(specify_decimal(mean(Days_since_germination, na.rm = TRUE), 4)),
                    sem = specify_decimal(sd(Days_since_germination, na.rm = TRUE)/sqrt(.N), 4)), by = .(Host_Species)][order(-mean)]
write.csv(x = flmeans, file = "./Data/Many_hosts/Model_outputs/Days_to_flower/means_sds_days.csv")

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
# lower CI
exp(4.6196626-0.2703376)
# upper CI
exp(4.6196626+0.0043234)

# wald test of functional group
write.csv(x = aod::wald.test(cov(mcmcfix2$Sol[,3:6, drop=F]), colMeans(mcmcfix2$Sol[,3:6, drop=F]), Terms=1:4)$result,
          file = "./Data/Many_hosts/Model_outputs/Days_to_flower/Functional_Group_Wald_Test.csv")

# effect of phylogeny
write.csv(list(posterior.mode(as.mcmc((mcmcfix2$VCV[,c(2)])/rowSums(mcmcfix2$VCV))),
               HPDinterval(as.mcmc((mcmcfix2$VCV[,c(2)])/rowSums(mcmcfix2$VCV)))), 
          file = "./Data/Many_hosts/Model_outputs/Days_to_flower/Phylogeny_Variance.csv")

# joint phylogenetic distribution
write.csv(list(posterior.mode(as.mcmc(rowSums(mcmcfix2$VCV[,c(1,2)])/rowSums(mcmcfix2$VCV))),
               HPDinterval(as.mcmc(rowSums(mcmcfix2$VCV[,c(1,2)])/rowSums(mcmcfix2$VCV)))), 
          file = "./Data/Many_hosts/Model_outputs/Days_to_flower/Joint_Phylogeny_Variance.csv")
# effect of species
write.csv(list(posterior.mode(as.mcmc((mcmcfix2$VCV[,c(1)])/rowSums(mcmcfix2$VCV))),
               HPDinterval(as.mcmc((mcmcfix2$VCV[,c(1)])/rowSums(mcmcfix2$VCV)))), 
          file = "./Data/Many_hosts/Model_outputs/Days_to_flower/Species_Variance.csv")
# units
write.csv(list(posterior.mode(as.mcmc((mcmcfix2$VCV[,c(3)])/rowSums(mcmcfix2$VCV))),
               HPDinterval(as.mcmc((mcmcfix2$VCV[,c(3)])/rowSums(mcmcfix2$VCV)))), 
          file = "./Data/Many_hosts/Model_outputs/Days_to_flower/Units_Variance.csv")


# variance components (redundant now?)
write.csv(x = daysvcv <- cbind(Response = "Days to flower",VCVapply(mcmcfix2)),
          file = "./Data/Many_hosts/Model_outputs/Days_to_flower/VCV_days.csv")


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
# for SI, replicates of each host.
fwrite(x = allgrowth[, .(Replicates = .N), by = .(Name)], file = "")

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
survivaldata.2$Time <- as.factor(survivaldata.2$Time)

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
posterior.sd(eha.1)
# posterior mode of the standard deviation of host species conditional on phylogeny
sqrt(posterior.mode(as.mcmc(rowSums(eha.1$VCV[,c(1,2)]))))
sqrt(HPDinterval(as.mcmc(rowSums(eha.1$VCV[,c(1,2)]))))
# posterior modes of the fixed effects
# annual/perennial
posterior.mode(as.mcmc(eha.1$Sol[,"AnnPerAnn"]))
HPDinterval(as.mcmc(eha.1$Sol[,"AnnPerAnn"]))
# and joint effect of functional group
posterior.mode(as.mcmc(c(eha.1$Sol[,"Functional_groupFern"], 
  eha.1$Sol[,"Functional_groupForb"],
  eha.1$Sol[,"Functional_groupLegume"],
  eha.1$Sol[,"Functional_groupWoody"])))
HPDinterval(as.mcmc(c(eha.1$Sol[,"Functional_groupFern"], 
                         eha.1$Sol[,"Functional_groupForb"],
                         eha.1$Sol[,"Functional_groupLegume"],
                         eha.1$Sol[,"Functional_groupWoody"])))


# print table to put into the manuscript:
write.csv(x = specify_decimal(summary(eha.1)$solutions, 4), 
          file = "./Data/Many_hosts/Model_outputs/Survival/EHA_solutions.csv")

# functional group wald test
write.csv(x = aod::wald.test(cov(eha.1$Sol[,5:8, drop=F]), colMeans(eha.1$Sol[,5:8, drop=F]), Terms=1:4)$result$chi2,
          file = "./Data/Many_hosts/Model_outputs/Survival/Functional_group_Wald_Test.csv")
# annual perennial
write.csv(x = aod::wald.test(cov(eha.1$Sol[,3, drop=F]), colMeans(eha.1$Sol[,3, drop=F]), Terms=1)$result$chi2,
          file = "./Data/Many_hosts/Model_outputs/Survival/AnnPer_Wald_Test.csv")
# covariance of functional group. Save
write.csv(x = cov(eha.1$Sol[,5:8, drop=F]), 
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

# save the variance components

write.csv(x = ehavcv <- cbind(Response = "Survival",VCVapply(eha.1)),
          file = "./Data/Many_hosts/Model_outputs/Survival/VCV_eha.csv")

##### Part 5: End of season reproductive nodes #####

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

# some summary stats

allgrowth3[, .(Mean.Nodes = mean(C.R.Nodes.F)), by = "Name"][order(-Mean.Nodes)]
# adjust transplant date
allgrowth3$Transplant.Date <- allgrowth3$Transplant.Date -97

# relevel factors
allgrowth3$AnnPer <- relevel(x = allgrowth3$AnnPer, ref = "Per")
allgrowth3$Functional_group <- relevel(x = allgrowth3$Functional_group, ref = "Grass")

# plot of the functional groups
opar <- par()
par(mfrow = c(1,2), oma = c(1, 0, 0, 0))

barplot_func <- allgrowth3[, .(mean = mean(C.R.Nodes.F), sem = sd(C.R.Nodes.F)/sqrt(.N)), by = .(Functional_group, Name)]
names_keep <- unique(barplot_func[mean > 2,]$Name)

barplot_func1 <- allgrowth3[, .(mean = mean(C.R.Nodes.F), sem = sd(C.R.Nodes.F)/sqrt(.N)), by = .(Functional_group)]
barplot_func2 <- allgrowth3[Name %in% names_keep, .(mean = mean(C.R.Nodes.F), sem = sd(C.R.Nodes.F)/sqrt(.N)), by = .(Functional_group)]

setorder(barplot_func1, mean)
setorder(barplot_func2, mean)

mid1 <- barplot(barplot_func1$mean, plot = FALSE)
mid2 <- barplot(barplot_func2$mean, plot = FALSE)

# save 6 x 10
barplot(barplot_func1$mean, names.arg = c("Woody", "Fern", "Forb", "Grass", "Legume"), ylim = c(0, 62), ylab = "Number of reproductive nodes")
arrows(x0 = mid1, y0 = barplot_func1$mean , x1 = mid1, barplot_func1$mean + barplot_func1$sem, code=3, angle=90, length=0.1)
text(1,59, expression(bold("a")), cex = 2)
barplot(barplot_func2$mean, names.arg = c("Woody", "Fern", "Forb", "Grass", "Legume"), ylim = c(0, 62))
arrows(x0 = mid2, y0 = barplot_func2$mean , x1 = mid2, barplot_func2$mean + barplot_func2$sem, code=3, angle=90, length=0.1)
text(1,59, expression(bold("b")), cex = 2)
mtext("Functional group", line = -2, side = 1, outer = TRUE, cex = 1.5)


par(opar)


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

# save the solutions
write.csv(x = specify_decimal(summary(mcmcfix3)$solutions, 4),
          file = "./Data/Many_hosts/Model_outputs/Reproductive_nodes_end/Reprod_solutions.csv")

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

# variance components
write.csv(x = nodesvcv <- cbind(Response = "End of season reproduction", VCVapply(mcmcfix3)),
          file = "./Data/Many_hosts/Model_outputs/Reproductive_nodes_end/VCV_nodes.csv")

varcomps <- rbind(daysvcv, ehavcv, nodesvcv)[Variable %in% c("Host_Species", "Name"), Variable := "Host Species"][Variable == "animal", Variable := "Phylogeny"]

varcomps2 <- cbind(varcomps[,c(1,2)], VCVglmm::specify_decimal(x = varcomps[,c(3:5)], k = 4))

write.csv(x = varcomps2, file = "./Data/Many_hosts/Model_outputs/varcomps.csv")


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

# some basic stats
# checking annual perennial interactive effect isnt not driven by one host
allgrowth7[, .(mean = mean(Nodes),
               sem = sd(Nodes)/sqrt(.N)), by = .(AnnPer, Name, Time)][Time == 4][order(mean)]
# mean and sem for each host at each time point
# May == 1, June == 2, July == 3, August == 4, September == 5
write.csv(x = allgrowth5[, .(MeanT = mean(Nodes), SEMT = paste("±", sd(Nodes)/sqrt(.N))), by = .(Name, Time)][
  Name %in% c("Lotus corniculatus", "Trifolium pratense", "Cynosurus cristatus")
][order(Name)], file = "./Data/Many_hosts/Model_outputs/Reproduction_over_time/example_time_points.csv")

prior.4<-list(R=list(V=diag(3), nu=0.002), 
              G=list(G1=list(V=diag(3), nu=3, alpha.mu=rep(0,3), alpha.V=diag(3)*1000),
                     G2= list(V=diag(1), nu=3, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))
# run for 13,000 X 20 iterations...
# I should add transplant date here?
mcmcfix4<-MCMCglmm(Nodes~Time*AnnPer+Functional_group + Transplant.Date,
                   random = ~ us(Time):Name + animal,
                   ginverse = list(animal=AinvULT2),
                   rcov=~idh(Time):units,
                   family="poisson", 
                   prior=prior.4, 
                   data=allgrowth7,
                   nitt=13000*50,
                   thin=10*50,
                   burnin=3000*50,
                   verbose = TRUE,
                   pr=TRUE)
summary(mcmcfix4)
# summary output
write.csv(x = specify_decimal(summary(mcmcfix4)$solutions, 4),
          file = "./Data/Many_hosts/Model_outputs/Reproduction_over_time/OT_Solutions.csv")
# functional group
write.csv(x = aod::wald.test(cov(mcmcfix4$Sol[,5:8, drop=F]), colMeans(mcmcfix4$Sol[,5:8, drop=F]), Terms=1:4)$result$chi2,
          file = "./Data/Many_hosts/Model_outputs/Reproduction_over_time/Functional_Group_Wald_Test.csv")
# life history host: time point
write.csv(x = aod::wald.test(cov(mcmcfix4$Sol[,9:10, drop=F]), colMeans(mcmcfix4$Sol[,9:10, drop=F]), Terms=1:2)$result$chi2,
          file = "./Data/Many_hosts/Model_outputs/Reproduction_over_time/AnnPerInt_Wald_Test.csv")
# life history of host
write.csv(x = aod::wald.test(cov(mcmcfix4$Sol[,4, drop=F]), colMeans(mcmcfix4$Sol[,4, drop=F]), Terms=1)$result$chi2,
          file = "./Data/Many_hosts/Model_outputs/Reproduction_over_time/AnnPer_Wald_Test.csv")

# expected number of nodes at time point 1 on perennial 

# not sure this is correct.
exp(3.06301) # = 21 times the expected number of nodes at time point 4 compared to time point 2
exp(0.78719) # 2.2 times the expected number of nodes compared to a perennial grass host at time 2
# interaction at time point four reduces expected number of nodes by 19 times
exp(3.06301+-2.33833)

#  intercept + time4 = expected number of nodes at time4 for perennial grass
exp(summary(mcmcfix4)$solutions[1,1] + #intercept
      summary(mcmcfix4)$solutions[3,1] +
      summary(mcmcfix4)$solutions[9,1]*mean(allgrowth7$Transplant.Date))/ # time 4

exp(summary(mcmcfix4)$solutions[1,1] + #intercept
      summary(mcmcfix4)$solutions[3,1] + # time 4
      summary(mcmcfix4)$solutions[4,1] +
      summary(mcmcfix4)$solutions[11,1] + 
      summary(mcmcfix4)$solutions[9,1]*mean(allgrowth7$Transplant.Date)) # interaction with annual/perrenial
# both of these mean the same thing, that there is a 4.7 times difference in expected number of nodes
# between annual and perennial plants at time point 4, roughly below.
exp(3.06301 +  0.78719 + -2.33833)
# CI 0.14-127.0

# which species are contributing to the significance of the time point 4 annual...
Solapply(mcmcfix4)[Group == "Time4"][order(-Grouped_Value)]


##### Plot 1: Event history analysis #####

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

surv_plot_dat <- survivaldata[to_add, on = "Unique_ID"][Name != "No host"]
surv_plot_dat[, y := (y-1)*-1]

# make time numeric for smooth interpolation
surv_plot_dat$Time <- as.factor(surv_plot_dat$Time)
# fill empty factor levels
completeDT <- function(DT, cols, defs = NULL){
  mDT = do.call(CJ, c(DT[, ..cols], list(unique=TRUE)))
  res = DT[mDT, on=names(mDT)]
  if (length(defs)) 
    res[, names(defs) := Map(replace, .SD, lapply(.SD, is.na), defs), .SDcols=names(defs)]
  res[]
} 
# fill them
surv_plot_dat2 <- completeDT(surv_plot_dat, c("Unique_ID", "Time"))
# invert y
surv_plot_dat2[is.na(y), y := 0 ]
# and carry forward last observation of Family
surv_plot_dat2[, Family := Family[nafill(replace(.I, is.na(Family), NA), "locf")]]
# and name
surv_plot_dat2[, Name := Name[nafill(replace(.I, is.na(Name), NA), "locf")]]

# this plot takes three largest families, two coincide with their own functional group
# May == 1, June == 2, July == 3, August == 4, September == 5

surv_plot_dat3 <- surv_plot_dat2[Family %in% c("Poaceae", "Fabaceae")]
surv_plot_dat3[, Family := ifelse(Family == "Poaceae", yes = "a", no = "b")]
# make time numeric for smooth interpolation
surv_plot_dat3$Time <- as.numeric(surv_plot_dat3$Time)

  (plot_1 <- ggplot(surv_plot_dat3, 
                 aes(x=Time, y=y,group=Name))+
  geom_point(alpha=0.1, position = position_jitter(width = 0.2, height = 0.2))+
  facet_wrap(~Family, scales = "free_x")+
  geom_smooth(aes(x=Time, y=y, col = Family), 
              method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE, size = 1)+
  ylim(c(0,1))+
  theme_bw()+
  ylab(label = "Probability of survival")+
  labs(colour = "Species")+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20, hjust = 0, face = "bold", family = "Helvetica"),
        legend.text = element_text(face = "italic"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
    scale_color_manual(values = alpha(c(cbPalette[8], cbPalette[4]), 0.3)) + 
    scale_x_continuous(breaks = 1:5, labels =  c("May", "June", "July", "August", "September"))+
    new_scale_colour()+
    scale_color_manual(values = alpha(c(cbPalette[8], cbPalette[4]), 1))+
    geom_smooth(aes(x=Time, y=y, col = Family, group = Family), 
                method = "glm", 
                method.args = list(family = "binomial"), 
                se = TRUE, size = 2))
  
  
ggsave(filename = "./Figures/Many_hosts/two_families_survival.pdf", plot = plot_1, 
       device = "pdf", width = 7, height = 6, units = "in")
ggsave(filename = "./Figures/Many_hosts/two_families_survival.jpeg", plot = plot_1, 
       device = "jpeg", width = 7, height = 6, units = "in")


##### Plot 2: Reproductive nodes at end of season with phylogeny #####

bootstrap_constraint_tr <- treeio::read.newick("./Data/Many_hosts/270220_bootstrap.newick", node.label = "support")
# remove the underscores
bootstrap_constraint_tr@phylo$tip.label <- gsub("_", " ", bootstrap_constraint_tr@phylo$tip.label) 
# correct the tip labels
bootstrap_constraint_tr@phylo$tip.label[10] <- "Hyacinthoides non-scripta"
bootstrap_constraint_tr@phylo$tip.label[36] <- "Chenopodium bonus-henricus"
bootstrap_constraint_tr@phylo$tip.label[43] <- "Cystopteris dickieana"

### ROOTING TREE ###

# root tree @ Cystopteris.
bootstrap_constraint_tr2 <- treeio::root(phy = bootstrap_constraint_tr@phylo, outgroup = "Cystopteris dickieana")
# edge cannot be zero, so make it tiny
bootstrap_constraint_tr2$edge.length[1] <- 1e-10

# get support values
support <- setDT(bootstrap_constraint_tr@data)
bootstrap_constraint_tr2$node.label[-1]
# manual editing.
show_boot <- c(rep(NA, 45), support$support, NA)

# export newick here
ape::write.tree(phy = phytools::force.ultrametric(bootstrap_constraint_tr2), file = "./Data/Tree/ultrametric_tree_vis.newick")

# constraint shows some orders
ggtree(phytools::force.ultrametric(constraint))+geom_nodelab()

(p <- ggtree(phytools::force.ultrametric(bootstrap_constraint_tr2))+
  geom_tiplab(size = 4, font = "italic") + theme(axis.line = element_blank()) + 
  xlim_tree(1.5)+
  geom_point2(aes(subset=!is.na(show_boot[-89]), fill = show_boot[-89]), pch=21, size = 3)+
  scale_fill_gradient(low = cbPalette[4], high = cbPalette[8], guide = "legend", name= "Bootstrap Support")+
  theme(legend.position = c(0.1, 0.85))+
  #geom_nodelab()+
  geom_cladelabel(node=54, label="Monocots", align=TRUE, angle = 270, offset =1, hjust = 0.5, offset.text = 0.05, barsize = 1.5)+
  geom_cladelabel(node=71, label="Asterales", align=TRUE, angle = 270, offset =1, hjust = 0.5, offset.text = 0.05, barsize = 1.5)+
  geom_cladelabel(node=83, label="Caryophyllales", align=TRUE, angle = 270, offset =1, hjust = 0.5, offset.text = 0.05, barsize = 1.5)+
  geom_cladelabel(node=64, label="Fabales", align=TRUE, angle = 270, offset =1, hjust = 0.5, offset.text = 0.05, barsize = 1.5))

p
# create the ggtree object
#d <- data.frame(label = bootstrap_constraint_tr2$tip.label)
#(p<-ggtree(phytools::force.ultrametric(constraint.2))+
#    geom_tiplab(size = 4)+ theme(axis.line = element_blank())+ # add higher taxonomic labels here.
#    geom_cladelabel(node=64, label="Legumes", align=TRUE, angle = 270, offset =1, hjust = 0.5, offset.text = 0.05, barsize = 1.5)+
#    geom_cladelabel(node=53, label="Grasses", align=TRUE, angle = 270, offset =1, hjust = 0.5, offset.text = 0.05, barsize = 1.5)
#   )
# data to go with the tree
phylodat<-allgrowth3[,c("Name","C.R.Nodes.F", "Functional_group", "Family")]
phylodat$group <- phylodat$Name

# group the data by host species
phylodat1.1 <- phylodat[, .(Mean.Cum = mean((na.omit(C.R.Nodes.F)+1)),
                            SD.Cum = sd((na.omit(C.R.Nodes.F)+1))/sqrt(.N),
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
    theme(legend.position = c(0.1, 0.85))+
    theme(axis.ticks.x = element_blank()))
# overlay the red observational points

colour_vector <- c(rep(cbPalette[7], 8), cbPalette[1], rep(cbPalette[7], 3), cbPalette[1],
                   rep(cbPalette[7], 5), rep(cbPalette[6], 6), rep(cbPalette[7], 5),
                   rep(cbPalette[4], 8), rep(cbPalette[7], 4), cbPalette[1], rep(cbPalette[5], 3))

(plot_2 <- facet_plot(finalplot, panel = "Cumulative Reproductive Nodes",
           data = phylodat1.1,
           mapping = aes(x=Mean.Cum, group = group),
           geom = geom_point, fill = rev(colour_vector), size = 4, pch=21)+
  theme_bw()+
  theme(strip.text.x = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(size = 15),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        legend.position = c(0.1, 0.85))+
  # add some hilights, the rest can be added in photoshop.
  geom_hilight(node=53, fill=cbPalette[4], alpha=.2, extendto = 1.5)+
  geom_hilight(node=64, fill=cbPalette[6], alpha=.2, extendto = 1.5)+
  geom_hilight(node=55, fill=cbPalette[7], alpha=.2, extendto = 1.5)+ # monocots but not grasses.
  geom_hilight(node=73, fill=cbPalette[7], alpha=.2, extendto = 1.5)+ # above mimulus
  geom_hilight(node=82, fill=cbPalette[7], alpha=.2, extendto = 1.5)+ # galium
  geom_hilight(node=83, fill=cbPalette[7], alpha=.2, extendto = 1.5)+  
  geom_hilight(node=30, fill=cbPalette[7], alpha=.2, extendto = 1.5)+ # mimulus
  geom_hilight(node=34, fill=cbPalette[1], alpha=.2, extendto = 1.5)+ #erica
  geom_hilight(node=29, fill=cbPalette[1], alpha=.2, extendto = 1.5)+ # thymus
  geom_hilight(node=31, fill=cbPalette[7], alpha=.2, extendto = 1.5)+ #plantago
  geom_hilight(node=44, fill=cbPalette[5], alpha=.2, extendto = 1.5) #plantago  
  )
# geom_hline(yintercept = c(8.5, 16.5, 21.5, 27.5), lty =2, alpha = 0.3)
# 54 is galanthus to dactylorhiza
# 60 between vicia and festuca
# pinus 87?
# 78 for erica?
# some edits in photoshop...

ggsave(filename = "./Figures/Many_hosts/Cum_reprod_nodes_phylo.pdf", plot = plot_2, 
       device = "pdf", width = 11, height = 12, units = "in", useDingbats=FALSE)
ggsave(filename = "./Figures/Many_hosts/Cum_reprod_nodes_phylo.jpeg", plot = plot_2, 
       device = "jpeg", width = 9, height = 7, units = "in")

##### Plot 3: Effect of annual/perennial over time (host life history) #####
# reproductive nodes over time (for functional group and annual or perennial!)
plot_3 <- allgrowth5[, .(Nodes = mean(Nodes), uCI = mean(Nodes)+sd(Nodes)/sqrt(.N), lCI = mean(Nodes)-sd(Nodes)/sqrt(.N)), by = c("AnnPer", "Time")] %>%
  ggplot(aes(x = Time, y= Nodes))+ 
  facet_wrap(~AnnPer)+
  geom_point()+
  geom_errorbar(aes(ymax = uCI, ymin =lCI), width = 0.05)+
  geom_line(aes(group = AnnPer))+
  theme_bw()+
  ylab(label = "")
  theme(strip.background = element_blank(), strip.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

ggsave(filename = "./Figures/Many_hosts/ann_per_time.pdf", plot = plot_3, 
       device = "pdf", width = 8, height = 7, units = "in")

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
  theme(strip.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20))

ggsave(filename = "./Figures/Many_hosts/average_host.pdf", plot = plot_4, 
       device = "pdf", width = 13.5, height = 7.5, units = "in")
ggsave(filename = "./Figures/Many_hosts/average_host.jpeg", plot = plot_4, 
       device = "jpeg", width = 13.5, height = 7.5, units = "in")

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
    theme_bw()+
    theme(panel.grid = element_blank())+
    geom_density_ridges()+geom_vline(xintercept = 0, col = "red", lty = 2, size = 1)+
    scale_x_continuous(limits = c(-0.2,1))

ggsave(filename = "./Figures/Many_hosts/phylogenetic_variance.pdf", plot = plot_5, 
       device = "pdf", width = 6, height = 5, units = "in")

