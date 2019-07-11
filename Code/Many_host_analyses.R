#### Fitness of a Euphrasia species on many hosts ####
# Created: 7.8.18 by Max Brown

## Libraries needed ##

library(plyr) 
library(dplyr) 
library(data.table)
library(tibble) 
library(ggplot2) 
library(lattice) 
library(tidyr) 
library(ape)
library(phangorn) 
library(ggtree)
library(ggstance) 
library(MCMCglmm)
library(lme4)
library(broom.mixed)
library(Hmisc)
library(ggrepel)

## functions used ##

# use ggplot interface to plot model of MCMCglmm fixed effects
MCMCfixplot <- function(mod){
  mod$Sol <- mod$Sol[,c(1:mod$Fixed$nfl)]
  rsum<-cbind(B = posterior.mode(mod$Sol[,1:dim(mod$Sol)[2]]), 
              CI = HPDinterval(mod$Sol[,1:dim(mod$Sol)[2]]))
  # put them in a data frame
  rsum<-data.frame(rsum)
  # create a vector of host names
  rsum$variable<- c(rownames(rsum))
  
  # find p values
  summ <- summary(mod)
  pvals <- summ$solutions[,5]
  # add p values to data frame
  rsum$pvals <- pvals
  
  # create a loop to assign asterisks to certain significance levels
  rsum$ptext <- rep(NA, dim(rsum)[1])
  for(i in 1:dim(rsum)[1]){
    if(rsum$pvals[i] <= 0.05 & rsum$pvals[i] > 0.01){
      rsum$ptext[i] <- "*"
    } else if(rsum$pvals[i] <= 0.01 & rsum$pvals[i] > 0.001){
      rsum$ptext[i] <- "**"
    } else if(rsum$pvals[i] <= 0.001 & rsum$pvals[i] > 0){
      rsum$ptext[i] <- "***"
    } else rsum$ptext[i] <- NA
  }
  
  # plot the graph with the confidence intervals, and add asterisks
  ggplot(rsum,aes(x = variable , y = B))+
    geom_pointrange(aes(ymin=lower,
                        ymax=upper))+
    geom_text(aes(label = ptext), nudge_x = 0.2) +
    coord_flip()+
    theme_bw()+
    labs(y = "Posterior modes", x = "Variables")
}
# use ggplot interface to plot MCMCglmm VCV objects
MCMCvarplot <- function(mod){
  rsum<-cbind(B = posterior.mode(mod$VCV[,1:dim(mod$VCV)[2]]), 
              CI = HPDinterval(mod$VCV[,1:dim(mod$VCV)[2]]))
  # put them in a data frame
  rsum<-data.frame(rsum)
  # create a vector of host names
  rsum$variable<- c(rownames(rsum))
  
  data <- as.data.frame(mod$VCV)
  # plot the graph with the confidence intervals, and add asterisks
  ggplot(rsum,aes(x = variable , y = B, fill = variable))+
    geom_pointrange(aes(ymin=lower,
                        ymax=upper))+
    coord_flip()+
    theme_bw()+
    labs(y = "Posterior modes", x = "Variables")
}
# This function is used to calculate the x and y axis for peak
# density on mcmc objects (or any other denisty). Taken from:
# http://ianmadd.github.io/pages/PeakDensityDistribution.html
densMode <- function(x){
  td <- density(x)
  maxDens <- which.max(td$y)
  maxx <- which.max(td$x)
  list(x=td$x[maxDens], y=td$y[maxDens], max.x=td$x[maxx])
}
# use ggplot to plot the random effect posterior distributions
# for group=TRUE check names in the data output.
MCMCranef <- function(mod, group = FALSE){
  if(group == TRUE) {
    ref <- data.frame(mod$Sol[,-c(1:mod$Fixed$nfl)])
    
    rdf <- reshape2::melt(ref)
    
    rdf <- rdf %>%
      mutate(group= gsub("\\..*", "", variable))
    
    xdat <- ddply(rdf,"group", transform, val_mean = signif(densMode(value)$x,3), med.x = signif(densMode(value)$x,3), med.y=signif(densMode(value)$y,3))
    xdat <-xdat[!duplicated(xdat$variable),]
    plot<-ggplot(rdf, aes(value, colour=group))+geom_density()+geom_label(data = xdat, aes(x=med.x,y=med.y,label=val_mean))
    
    print(xdat);print(plot)
  } else {
    ref <- data.frame(mod$Sol[,-c(1:mod$Fixed$nfl)])
    
    rdf <- reshape2::melt(ref)
    
    xdat <- ddply(rdf,"variable", transform, val_mean = signif(densMode(value)$x,3), med.x = signif(densMode(value)$x,3), med.y=signif(densMode(value)$y,3))
    xdat <-xdat[!duplicated(xdat$variable),]
    plot<-ggplot(rdf, aes(value, colour=variable))+geom_density()+geom_label(data = xdat, aes(x=med.x,y=med.y,label=val_mean))+theme(legend.position = "none")
    
    print(xdat);print(plot)
  }
}
# A function modified to the above to calculate the joint 
# posterior of phylogenetic and non phylogenetic effects
# This visualises the random effect point estimates irrespective of whether its phylogenetic or not
MCMCranef_joined <- function(mod){
  ref <- data.frame(mod$Sol[,-c(1:mod$Fixed$nfl)])
  
  n <- data.frame(matrix(ncol = (dim(ref)[2]/2), nrow = 1000))
  for(i in 1:(dim(ref)[2]/2)){
    n[,i]<-rowSums(ref[,c(i, i+(dim(ref)[2]/2))])
  }
  colnames(n) <- unlist(mod$ginverse$animal@Dimnames[1])[order(unlist(mod$ginverse$animal@Dimnames[1]))]
  
  rdf <- reshape2::melt(n)
  xdat <- ddply(rdf,"variable", transform, val_mean = signif(densMode(value)$x,3), med.x = signif(densMode(value)$x,3), med.y=signif(densMode(value)$y,3))
  xdat <-xdat[!duplicated(xdat$variable),]
  xdat$variable <- xdat$variable[c(T, NA)]
  plot <- ggplot(rdf, aes(value, colour=variable))+geom_density()+geom_label(data = xdat, aes(x=med.x,y=med.y,label=variable))+theme(legend.position = "none")
  print(xdat);print(plot)
}
post.phy.ranefs <- function(mod){
  ref <- data.frame(mod$Sol[,-c(1:mod$Fixed$nfl)])
  
  n <- data.frame(matrix(ncol = (dim(ref)[2]/2), nrow = 1000))
  for(i in 1:(dim(ref)[2]/2)){
    n[,i]<-rowSums(ref[,c(i, i+(dim(ref)[2]/2))])
  }
  colnames(n) <- unlist(mod$ginverse$animal@Dimnames[1])[order(unlist(mod$ginverse$animal@Dimnames[1]))]
  n <- as.mcmc(n)
  rsum<-cbind(B = posterior.mode(n[,1:dim(n)[2]]), 
              CI = HPDinterval(n[,1:dim(n)[2]]))
  rsum<-data.frame(rsum)
  rsum$Names <- gsub(pattern = ".", replacement = " ", fixed = T, x = rownames(rsum))
  ggplot(rsum,aes(x = reorder(Names, B) , y = B))+
    geom_pointrange(aes(ymin=lower,
                        ymax=upper))+
    coord_flip()+
    theme_bw()+
    labs(y = "Posterior modes", x = "Variables")
}
# visualise the variance covariance matrix values
VCVdensity <- function(mod){
  ref <- data.frame(mod$VCV)
  rdf <- reshape2::melt(ref)
  ggplot(rdf, aes(x=value, fill=variable))+geom_density(alpha=0.3)
}
# use lattice plot interface to see traces for $Sol and $VCV objects
traces <- function(mod){
  par.settings = list(strip.background=list(col="lightgrey"))
  trace.1<-xyplot(mod$Sol,
                  mean.list = mod$Sol,
                  panel = function(x, y, mean.list) {
                    panel.lines(x, y, col = "black")
                    panel.abline(h = mean(y),
                                 col.line = "red")
                  },
                  par.settings=par.settings)
  trace.2<-xyplot(mod$VCV,
                  mean.list = mod$VCV,
                  panel = function(x, y, mean.list) {
                    panel.lines(x, y)
                    panel.abline(h = mean(y),
                                 col.line = "red")
                  },
                  par.settings=par.settings)
  Newlist <- list(trace.1, trace.2)
  return(Newlist)
}
# Use this function to find variance explained by each component of the model. By Greg Albery.
MCMCRep <- function(Model,scale="original"){
  require(MCMCglmm)
  if(scale=="original"){
    if(unique(Model$family)=="poisson"){
      Beta0<-sapply(1:dim(Model$Sol)[1],function(z) mean(as.matrix(Model$X)%*%as.matrix(Model$Sol[z,1:ncol(Model$X)]))) #Intercept
      mat<-matrix(NA,nrow=dim(Model$VCV)[2],ncol=4) #Setting up a matrix for the VCV
      for(j in 1:dim(Model$VCV)[2]){
        Va<-Model$VCV[,j] #Associated variance
        Ve<-rowSums(Model$VCV) #Inclusive Error variance
        Expected<-exp(Beta0+(0.5*(Ve))) #Expected values
        Repeatability1<-(Expected*(exp(Va)-1))/(Expected*(exp(Ve)-1)+1) #Repeatability
        mat[j,]<-c(colnames(Model$VCV)[j],posterior.mode(Repeatability1),HPDinterval(Repeatability1)[,1],HPDinterval(Repeatability1)[,2])
      }
    }
    if(unique(Model$family)=="gaussian"){
      mat<-matrix(NA,nrow=dim(Model$VCV)[2],ncol=4)
      for(j in 1:dim(Model$VCV)[2]){
        Va<-Model$VCV[,j] #Associated variance
        Ve<-rowSums(Model$VCV)
        Repeatability1<-(Va/Ve)
        mat[j,]<-c(colnames(Model$VCV)[j],round(posterior.mode(Repeatability1),digits=2),round(HPDinterval(Repeatability1)[,1],digits=2),round(HPDinterval(Repeatability1)[,2],digits=2))
      }
    }
  }
  if(scale=="link"){
    if(unique(Model$family)=="poisson"){
      Beta0<-sapply(1:dim(Model$Sol)[1],function(z) mean(as.matrix(Model$X)%*%as.matrix(Model$Sol[z,1:ncol(Model$X)]))) #Intercept
      mat<-matrix(NA,nrow=dim(Model$VCV)[2],ncol=4)
      for(j in 1:dim(Model$VCV)[2]){
        Va<-Model$VCV[,j] #Associated variance
        Ve<-rowSums(Model$VCV)
        Repeatability1<-Va/(Ve+log(1/exp(Beta0)+1))
        mat[j,]<-c(colnames(Model$VCV)[j],round(posterior.mode(Repeatability1),digits=2),round(HPDinterval(Repeatability1)[,1],digits=2),round(HPDinterval(Repeatability1)[,2],digits=2))
      }
    }
    if(unique(Model$family)=="gaussian"){
      mat<-matrix(NA,nrow=dim(Model$VCV)[2],ncol=4)
      for(j in 1:dim(Model$VCV)[2]){
        Va<-Model$VCV[,j] #Associated variance
        Ve<-rowSums(Model$VCV)
        Repeatability1<-(Va/Ve)
        mat[j,]<-c(colnames(Model$VCV)[j],round(posterior.mode(Repeatability1),digits=2),round(HPDinterval(Repeatability1)[,1],digits=2),round(HPDinterval(Repeatability1)[,2],digits=2))      }
    }
  }
  colnames(mat)<-c("Component","Mode","lHPD","uHPD")
  data.frame(mat)
}
# Variance explained for normal distribution, or those on the link scale
MCMCRepnorm <- function(mod, y = "variable"){
  var.a      <- mod$VCV[,y]
  var.e      <- rowSums(mod$VCV)
  postR.link <- var.a/(var.e)
  #EY 		   <- exp(beta0+((var.e+var.a)/2))
  #postR.org  <- EY*(exp(var.a)-1)/(EY*(exp(var.e+var.a)-1)+1)
  # point estimate on link and original scale
  R.link     <- posterior.mode( postR.link )
  #R.org      <- posterior.mode( postR.org )
  # CI's
  CI.link    <- coda::HPDinterval(postR.link)[1,]
  #CI.org     <- coda::HPDinterval(postR.org)[1,]
  
  res 	   <- list(R.link=R.link, CI.link=CI.link)
  #R.org = R.org, CI.org=CI.org)
  res.2 <- unlist(res)*100
  names(res.2) <- c("Point Estimate", "CI Lower", "CI Upper")
  #"Original Scale Point Estimate", "CI Original Scale Lower", "CI Original Scale Upper")
  res.2 <- data.frame(res.2)
  colnames(res.2)[1] <- "%"
  return(res.2)
}
# ultrametric trees for Ainv output
ultm<-function(tree, MPL = FALSE, inverseA = TRUE){
  if(MPL == TRUE){
    MPL <- chronoMPL(tree)
  } else MPL <- tree
  
  ult <- phytools::force.ultrametric(MPL)
  # remove zeroes in the distance matrix
  for(i in 1:length(ult$edge.length)){
    if(ult$edge.length[i] < 1e-16){
      ult$edge.length[i] <- 1e-10
    }
  }
  if(inverseA == TRUE){
    inverseult <- MCMCglmm::inverseA(ult, nodes = "ALL", scale = TRUE)
    Ainv <- as(as.matrix(inverseult$Ainv), "dgCMatrix")
    return(Ainv)
  } 
  return(ult)
}
# manipulating logits and probs
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
probit2prob <- function(probit){
  pnorm(probit)
}


## notes ##

# Format is as follows. Simple stats for each dataset. The model then follows.

# The model is summarised, using summary() and the table is used for the supplementary material.
# Interesting effect sizes and directions can be eyeballed from the model.

# Variance explained by a) host species & phylogeny b) phylogeny are estimated.

# Wald tests of fixed effects for significance of functional group and annual/perennial.

# The phylogeny is constrained using the reference tree


##### tree 1 #####
# constraint tree output
constraint <- read.tree(text = "(Agrostis_capillaris:0.0103868960,Lagurus_ovatus:0.0137334981,((((Cynosurus_cristatus:0.0119501076,Festuca_rubra:0.0143927581)56:0.0014297595,Holcus_lanatus:0.0116621244)61:0.0018419113,Phleum_pratense:0.0106048521)100:0.0080082346,((Zea_mays:0.0376747406,((((Allium_ursinum:0.0696170091,Galanthus_nivalis:0.0548052238)amaryllidaceae/100:0.0011547212,Hyacinthoides_non_scripta:0.0509074406)100:0.0322710961,Dactylorhiza_purpurella:0.0989415971)asparagales/100:0.0128940050,((((((Sorbus_aucuparia:0.0314958348,Fragaria_vesca:0.0582563355)rosaceae/100:0.0349788544,((((Ononis_spinosa:0.0220743766,(Lathyrus_japonicus:0.0294240865,Trifolium_pratense:0.0339786483)100:0.0021655535)100:0.0183607112,Lotus_corniculatus:0.0850694608)100:0.0000019136,Ulex_europaeus:0.0718966448)100:0.0000027395,Vicia_cracca:0.0479406736)fabaceae/100:0.0785149642)100:0.0148394884,(Arabidopsis_thaliana:0.1250612607,Helianthemum_nummularium:0.1469119452)malvales_to_brassicales/100:0.0195592609)100:0.0150299043,((((((Meum_athamanticum:0.0221268430,Anthriscus_sylvestris:0.0221983048)apiaceae/100:0.0502387307,Centranthus_ruber:0.1162127616)100:0.0022623680,((Centaurea_nigra:0.0172220154,(Leucanthemum_vulgare:0.0331192970,Tragopogon_pratensis:0.0182447955)100:0.0000023339)100:0.0048658069,Senecio_vulgaris:0.0209331617)asteraceae/100:0.0535848492)100:0.0118251564,(((Thymus_polytrichus:0.0579324847,Mimulus_guttatus:0.0360518925)100:0.0082342727,Plantago_lanceolata:0.0899883709)100:0.0336066248,(Galium_aparine:0.0100601356,Galium_verum:0.0102328988)galium/100:0.1193345682)100:0.0185788636)100:0.0088343295,Erica_tetralix:0.1160061555)ericales_to_asterales/100:0.0079753034,(((Chenopodium_album:0.0332635734,Chenopodium_bonus_henricus:0.0105604686)chenopodium/100:0.0238607464,(Silene_dioica:0.0014311631,Silene_latifolia:0.0000024274)silene/100:0.0808030168)100:0.0364786972,Rumex_acetosella:0.1210989875)caryophyllales/100:0.0287067285)100:0.0088664052)100:0.0236664505,Papaver_rhoeas:0.1015330721)eudicots/100:0.0318216294,(Pinus_sylvestris:0.2367123693,(Equisetum_arvense:0.3878874112,(Cystopteris_fragilis:0.1233674528,Pteridium_aquilinum:0.1016667919)100:1.0445751830)100:0.2401868160)seedplants/100:0.1314408258)poales_to_asterales/100:0.0454413537)100:0.1494656671)100:0.0242458362,Hordeum_vulgare:0.0237302454)93:0.0060960686)poaceae/100:0.0096602096);")


constraint$tip.label <- gsub("_", " ", constraint$tip.label) 
constraint$tip.label[10] <- "Hyacinthoides non-scripta"
constraint$tip.label[36] <- "Chenopodium bonus-henricus"
constraint$tip.label[43] <- "Cystopteris dickieana"

# constraint.1 remember
constraint.1<-drop.tip(constraint, tip = sptodel)
constraint.1$node.label <- NULL
# root tree
constraint.1 <- root(phy = constraint.1, outgroup = "Cystopteris dickieana", resolve.root = T)
constraint.1$edge.length[1] <- 1e-10

varcons1<-mean(diag(vcv(constraint.1)))

INphylo2 <- inverseA(constraint.1, nodes = "ALL", scale = FALSE) # only scale if tree1 <- chronoMPL(tree1)
# we want the Ainv object

Ainv2 <- INphylo2$Ainv
# convert to matrix

Ainv2<- as.matrix(Ainv2)
# convert to a sparse matrix

Ainv2 <- as(Ainv2, "dgCMatrix")


# comparison of chronoMPL, ultrametric and the difference it makes
# the second is the more realistic and not different from 
# the mean(diag(vcv(tree)))*VCV`animal` method
AinvULT <- ultm(constraint.1, MPL = F, inverseA = T)



##### tree 2 #####
# second tree

# load in tree
constraint.2 <- read.tree(text = "(Agrostis_capillaris:0.0103868960,Lagurus_ovatus:0.0137334981,((((Cynosurus_cristatus:0.0119501076,Festuca_rubra:0.0143927581)56:0.0014297595,Holcus_lanatus:0.0116621244)61:0.0018419113,Phleum_pratense:0.0106048521)100:0.0080082346,((Zea_mays:0.0376747406,((((Allium_ursinum:0.0696170091,Galanthus_nivalis:0.0548052238)amaryllidaceae/100:0.0011547212,Hyacinthoides_non_scripta:0.0509074406)100:0.0322710961,Dactylorhiza_purpurella:0.0989415971)asparagales/100:0.0128940050,((((((Sorbus_aucuparia:0.0314958348,Fragaria_vesca:0.0582563355)rosaceae/100:0.0349788544,((((Ononis_spinosa:0.0220743766,(Lathyrus_japonicus:0.0294240865,Trifolium_pratense:0.0339786483)100:0.0021655535)100:0.0183607112,Lotus_corniculatus:0.0850694608)100:0.0000019136,Ulex_europaeus:0.0718966448)100:0.0000027395,Vicia_cracca:0.0479406736)fabaceae/100:0.0785149642)100:0.0148394884,(Arabidopsis_thaliana:0.1250612607,Helianthemum_nummularium:0.1469119452)malvales_to_brassicales/100:0.0195592609)100:0.0150299043,((((((Meum_athamanticum:0.0221268430,Anthriscus_sylvestris:0.0221983048)apiaceae/100:0.0502387307,Centranthus_ruber:0.1162127616)100:0.0022623680,((Centaurea_nigra:0.0172220154,(Leucanthemum_vulgare:0.0331192970,Tragopogon_pratensis:0.0182447955)100:0.0000023339)100:0.0048658069,Senecio_vulgaris:0.0209331617)asteraceae/100:0.0535848492)100:0.0118251564,(((Thymus_polytrichus:0.0579324847,Mimulus_guttatus:0.0360518925)100:0.0082342727,Plantago_lanceolata:0.0899883709)100:0.0336066248,(Galium_aparine:0.0100601356,Galium_verum:0.0102328988)galium/100:0.1193345682)100:0.0185788636)100:0.0088343295,Erica_tetralix:0.1160061555)ericales_to_asterales/100:0.0079753034,(((Chenopodium_album:0.0332635734,Chenopodium_bonus_henricus:0.0105604686)chenopodium/100:0.0238607464,(Silene_dioica:0.0014311631,Silene_latifolia:0.0000024274)silene/100:0.0808030168)100:0.0364786972,Rumex_acetosella:0.1210989875)caryophyllales/100:0.0287067285)100:0.0088664052)100:0.0236664505,Papaver_rhoeas:0.1015330721)eudicots/100:0.0318216294,(Pinus_sylvestris:0.2367123693,(Equisetum_arvense:0.3878874112,(Cystopteris_fragilis:0.1233674528,Pteridium_aquilinum:0.1016667919)100:1.0445751830)100:0.2401868160)seedplants/100:0.1314408258)poales_to_asterales/100:0.0454413537)100:0.1494656671)100:0.0242458362,Hordeum_vulgare:0.0237302454)93:0.0060960686)poaceae/100:0.0096602096);")
constraint.2$tip.label <- gsub("_", " ", constraint.2$tip.label) 
constraint.2$tip.label[10] <- "Hyacinthoides non-scripta"
constraint.2$tip.label[36] <- "Chenopodium bonus-henricus"
constraint.2$tip.label[43] <- "Cystopteris dickieana"
constraint.2$node.label <- NULL


# root tree
constraint.2 <- root(phy = constraint.2, outgroup = "Pteridium aquilinum", resolve.root = T)
constraint.2$edge.length[1] <- 1e-10

varcons2<-mean(diag(vcv(constraint.2)))

INphylo <- inverseA(constraint.2, nodes = "ALL", scale = FALSE)
# we want the Ainv object
Ainv<-INphylo$Ainv
# convert to matrix
Ainv<-as.matrix(Ainv)
# convert to a sparse matrix
Ainv <- as(Ainv, "dgCMatrix")

AinvULT2 <- ultm(constraint.2, MPL = F, inverseA = T)


## sort out all datasets up here ##

##### At first flowering #####

# load data
floweringsofar<- read.csv("/Users/mbrown/Documents/Edinburgh Ph.D./GROWTH EXPT DATA/FloweringsofarV2.csv", na.strings = "-")
setDT(floweringsofar)
# some data preparation first
# create host no host variable 
floweringsofar$Host.No.Host<-as.numeric(floweringsofar$ID_No > 0)
# relevel so perennial is the baseline
floweringsofar$AnnPer <- relevel(floweringsofar$AnnPer, ref = "Per")
# relevel so grass is baseline
floweringsofar$Functional.group <- relevel(floweringsofar$Functional.group, ref = "Grass")
# make host no host a factor
floweringsofar$Host.No.Host<- as.factor(floweringsofar$Host.No.Host)
# relevel so that 1 is the host baseline
floweringsofar$Host.No.Host <- relevel(floweringsofar$Host.No.Host, ref = "1")
# create 'animal' the phylogenetic random effect variable
floweringsofar$animal <- floweringsofar$Host.Species

# summary stats for corolla length

# 
summary(floweringsofar[, .(Mean.Corolla = mean(na.omit(Corolla.Length))), by = "Host.Species"][order(-Mean.Corolla)])

# corolla length
prior.1 <- list(R=list(V=diag(1), nu=0.002), 
                G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                       G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

# fixed effects estimates for both annual/perennial and functional group
mcmcfix1 <- MCMCglmm(Corolla.Length ~ AnnPer + Functional.group, 
                     random = ~Host.Species + animal, 
                     ginverse = list(animal = AinvULT),
                     data = floweringsofar[floweringsofar$Host.Species != "No host",],
                     prior = prior.1,
                     family= "gaussian",
                     nitt = 13000*10,
                     burnin = 3000*10,
                     thin = 10*10,
                     pr=T,
                     verbose = T)


# summary output
summary(mcmcfix1)

write.csv(x = summary(mcmcfix1)$solutions,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Corolla_length/Corolla_solutions.csv")

# joint phylogenetic distribution
# 66.5%
write.csv(
  list(posterior.mode(as.mcmc(rowSums(mcmcfix1$VCV[,c(1,2)])/rowSums(mcmcfix1$VCV))),
       # 18.5% - 80.4%
       HPDinterval(as.mcmc(rowSums(mcmcfix1$VCV[,c(1,2)])/rowSums(mcmcfix1$VCV)))), 
  file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Corolla_length/Joint_Phylogeny_Variance.csv")

# and the phylogenetic component
# 60.2%
write.csv(
  list(posterior.mode(as.mcmc((mcmcfix1$VCV[,c(2)])/rowSums(mcmcfix1$VCV))),
       # 9.41% - 80.6%
       HPDinterval(as.mcmc((mcmcfix1$VCV[,c(2)])/rowSums(mcmcfix1$VCV)))), 
  file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Corolla_length/Phylogeny_Variance.csv")


# host variance
posterior.mode(as.mcmc((mcmcfix1$VCV[,c(1)])/rowSums(mcmcfix1$VCV)))
HPDinterval(as.mcmc((mcmcfix1$VCV[,c(1)])/rowSums(mcmcfix1$VCV)))


# wald test of annual/perennial
write.csv(x = aod::wald.test(cov(mcmcfix1$Sol[,2, drop=F]), colMeans(mcmcfix1$Sol[,2, drop=F]), Terms=1)$result, 
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Corolla_length/AnnPer_Wald_Test.csv")
# wald test of functional group
write.csv(x = aod::wald.test(cov(mcmcfix1$Sol[,3:6, drop=F]), colMeans(mcmcfix1$Sol[,3:6, drop=F]), Terms=1:4)$result, 
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Corolla_length/Functional_Group_Wald_Test.csv")

# days to flower

prior.2 <- list(R=list(V=diag(1), nu=0.002), 
                G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                       G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

# fixed effects estimates for both annual/perennial and functional group
# the random effects specification allows different variances for each functional group
# but no covariances
mcmcfix2<- MCMCglmm(Days.since.germination ~ AnnPer + Functional.group, 
                    random = ~ Host.Species + animal, 
                    ginverse = list(animal=AinvULT),
                    data = floweringsofar[floweringsofar$Host.Species != "No host",],
                    prior = prior.2,
                    family= "poisson",
                    nitt = 13000*10,
                    burnin = 3000*10,
                    thin = 10*10,
                    pr=T,
                    verbose = T)


# summary output
summary(mcmcfix2)

write.csv(x = summary(mcmcfix2)$solutions,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Days_to_flower/Days_solutions.csv")

# wald test of annual/perennial
write.csv(x = aod::wald.test(cov(mcmcfix2$Sol[,2, drop=F]), colMeans(mcmcfix2$Sol[,2, drop=F]), Terms=1)$result,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Days_to_flower/AnnPer_Wald_Test.csv")
# wald test of functional group
write.csv(x = aod::wald.test(cov(mcmcfix2$Sol[,3:6, drop=F]), colMeans(mcmcfix2$Sol[,3:6, drop=F]), Terms=1:4)$result,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Days_to_flower/Functional_Group_Wald_Test.csv")

# joint phylogenetic distribution
# 33.3%
write.csv(list(posterior.mode(as.mcmc(rowSums(mcmcfix2$VCV[,c(1,2)])/rowSums(mcmcfix2$VCV))),
               # 21.5% - 87.2%
               HPDinterval(as.mcmc(rowSums(mcmcfix2$VCV[,c(1,2)])/rowSums(mcmcfix2$VCV)))), 
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Days_to_flower/Joint_Phylogeny_Variance.csv")

# and the phylogenetic component
# 0.451%
write.csv(list(posterior.mode(as.mcmc((mcmcfix2$VCV[,c(2)])/rowSums(mcmcfix2$VCV))),
               # 1.048705e-05% - 82.7%
               HPDinterval(as.mcmc((mcmcfix2$VCV[,c(2)])/rowSums(mcmcfix2$VCV)))),
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Days_to_flower/Phylogeny_Variance.csv")


# host variance
posterior.mode(as.mcmc((mcmcfix2$VCV[,c(1)])/rowSums(mcmcfix2$VCV)))
HPDinterval(as.mcmc((mcmcfix2$VCV[,c(1)])/rowSums(mcmcfix2$VCV)))

ggplot(as.data.frame((mcmcfix2$VCV[,c(2)])/rowSums(mcmcfix2$VCV)), aes(x = var1))+geom_density(fill="black")+
  theme_bw()

##### event history analysis ######

# May == 1, June == 2, July == 3, August == 4, September == 5
survivaldata<- read.csv("/Users/mbrown/OneDrive - University of Edinburgh/Euphrasia Experiment 1 HOSTS/GROWTH EXPT DATA/Survivalanalysis3.csv")
setDT(survivaldata)
# add time as a column in the data, i.e. at which time points was each individual Euphrasia alive
survivaldata <- survivaldata[, Time := .(Time = 1:.N), by="Unique_ID"][Time < 6,]
table(survivaldata$Time, survivaldata$y)

# SPELL DICKIEANA PROPERLY

# Average proability of death on each host
survivaldata[, .(Mean = mean((y-1)*-1 )), by = c("Name")][order(Mean)]
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
survivaldata$Host.No.Host<-as.numeric(survivaldata$ID_No > 0)
tmp <- as.POSIXlt(survivaldata$Latest.Germ.Date, format = "%d/%m/%y")
survivaldata$Latest.Germ.Date<-tmp$yday

tmp <- as.POSIXlt(survivaldata$Transplant.Date, format = "%d/%m/%y")
survivaldata$Transplant.Date<-tmp$yday

survivaldata$Host.No.Host<- as.factor(survivaldata$Host.No.Host)
survivaldata$Host.No.Host <- relevel(survivaldata$Host.No.Host, ref = "1")
survivaldata$Name <- relevel(survivaldata$Name, ref = "No host")

# model 
survivaldata.1<- data.frame(survivaldata.1) #needed!!
survivaldata.1$animal <- survivaldata.1$Name
levels(survivaldata.1$animal)[11] <- "Cystopteris dickieana"
levels(survivaldata.1$animal)[44] <- "Ulex europaeus"

tidy.g.n$Host.No.Host <- relevel(tidy.g.n$Host.No.Host, ref = "Host")
tidy.g.n$Time <- as.integer(tidy.g.n$Time)
tidy.g.n$Time <- relevel(tidy.g.n$Time, ref = "2")
tidy.g.n$AnnPer <- relevel(tidy.g.n$AnnPer, ref = "Per")
tidy.g.n$Functional.group <- relevel(tidy.g.n$Functional.group, ref = "Grass")


levels(tidy.g.n$Time) <- c("1", "2", "3", "4", "5")

survivaldata.2<-dplyr::left_join(survivaldata.1, unique(tidy.g.n[,c(1,3,4,9,10)]), by = c("Unique_ID", "Time"))
survivaldata.2$Time <- factor(survivaldata.2$Time)
survivaldata.2 <- filter(survivaldata.2, as.numeric(Time) <= 4)


prior.eha <-list(R=list(V=diag(1), fix=1), 
                 G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                        G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))
levels(survivaldata.2$Name)[11] <- "Cystopteris dickieana"
levels(survivaldata.2$Name)[44] <- "Ulex europaeus"
survivaldata.2$animal <- survivaldata.2$Name

survivaldata.2$y2 <- (survivaldata.2$y -1)*-1 # probability of survival

eha.1 <- MCMCglmm(y2 ~ 1 + Time+AnnPer+ Transplant.Date + Functional.group + Nodes,
                  random = ~Name + animal,
                  ginverse = list(animal=AinvULT2),
                  data = survivaldata.2[survivaldata.2$animal != "No host",],
                  prior=prior.eha,
                  family = "threshold",
                  nitt = 13000*10,
                  burnin = 3000*10,
                  thin = 10*10,
                  verbose = T,
                  pr=T)

summary(eha.1)
# print table to put into the manuscript:
write.csv(x = summary(eha.1)$solutions, 
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Event_History_Analysis/EHA_solutions.csv")

plogis(3.615433-1.240765) #  when time is 2, there is 0.914 probability of surival
plogis(3.615433-5.573162) # at time 4, there is a 0.124 probability of survival
exp(-3.429357) # this is the odds ratio. Pr success / Pr failure. So it is 0.032 dead for every one alive.
# note that these two are at the intercept, so this is when there are no reproductive nodes, for a grass host, at time point 1,
# and with the earliest transplant date

# transplant date for example
plogis(0.017166) # for every 1 unit increase in intercept, 1.017004 increase in the odds ratio
logit2prob(0.016861) # that means 50% increase in probability of death for each increase in intercept

logit2prob(-3.429357 + 104*0.016861)

Wald.test.auto <- function(mod){
  if (attributes(mod)$class != "MCMCglmm") {
    stop("Object not of class MCMCglmm")
  }
  fixeff <- as.character(mod$Fixed$formula)[3]
  fixeff2 <- unlist(strsplit(fixeff, "\\+"))
  covars <- gsub(fixeff2, pattern = " ", replacement = "")
  
  if(any(covars == "1")){
    sset <- grep(pattern = "[^1]", covars)
    covars <- covars[sset]
  }
  
  res <- list()
  for (i in 1:length(covars)) {
    fixed <- grep(x = rownames(summary(mod)$solutions), pattern = covars[i])
    varcov <- cov(mod$Sol[, fixed, drop = FALSE])
    coefs <- colMeans(mod$Sol[, fixed, drop = FALSE])
    Terms <- 1:length(fixed)
    res[[i]] <- aod::wald.test(Sigma = varcov, b = coefs, 
                               Terms = Terms)$result
    names(res[[i]]) <- covars[i]
  }
  return(data.frame(res))
}
Wald.test.auto(eha.1)

# functional group wald test
write.csv(x = aod::wald.test(cov(eha.1$Sol[,7:10, drop=F]), colMeans(eha.1$Sol[,7:10, drop=F]), Terms=1:4)$result$chi2,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Event_History_Analysis/Functional_group_Wald_Test.csv")
write.csv(x = aod::wald.test(cov(eha.1$Sol[,5, drop=F]), colMeans(eha.1$Sol[,5, drop=F]), Terms=1)$result$chi2,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Event_History_Analysis/AnnPer_Wald_Test.csv")
# covariance of functional group. Save
write.csv(x = cov(eha.1$Sol[,7:10, drop=F]), 
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Event_History_Analysis/Var_Covar_Functional_group.csv")


# joint phylogenetic distribution
# 27.1%
write.csv(x = list(posterior.mode(as.mcmc(rowSums(eha.1$VCV[,c(1,2)])/rowSums(eha.1$VCV))),
                   HPDinterval(as.mcmc(rowSums(eha.1$VCV[,c(1,2)])/rowSums(eha.1$VCV)))),
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Event_History_Analysis/Joint_Phylogeny_Variance.csv")
# 18.1% - 64.3%
HPDinterval(as.mcmc(rowSums(eha.1$VCV[,c(1,2)])/rowSums(eha.1$VCV)))

# and the phylogenetic component
# 0.512%
posterior.mode(as.mcmc((eha.1$VCV[,c(2)])/rowSums(eha.1$VCV)))*100
# 6.800695e-04% - 54.2%
HPDinterval(as.mcmc((eha.1$VCV[,c(2)])/rowSums(eha.1$VCV)))*100


# host variance
posterior.mode(as.mcmc((eha.1$VCV[,c(1)])/rowSums(eha.1$VCV)))
HPDinterval(as.mcmc((eha.1$VCV[,c(1)])/rowSums(eha.1$VCV)))


# visualise the fixed effects
MCMCfixplot(eha.1)
# visualise the posterior distributions of the random effect means
MCMCranef(eha.1)

# as above but grouped by time
MCMCranef(eha.1, group = TRUE)
# Variances
VCVdensity(eha.1)

##### end of season cumulative nodes #####

allgrowth<- read.csv("/Users/mbrown/OneDrive - University of Edinburgh/Euphrasia Experiment 1 HOSTS/GROWTH EXPT DATA/Allgrowthmeasurements1.csv", na.strings = "-")

# take columns 9 to 36 and collapse into one column
tidy.g<- gather(allgrowth, "Variable", "Value", 9:36)
# filter out the nodes
tidy.g.n<- filter(tidy.g, Variable %in% c("R.Nodes.1", "R.Nodes.2",
                                          "R.Nodes.3", "R.Nodes.4",
                                          "R.Nodes.5"))
# first change variable names to characters
tidy.g.n$Variable <- as.character(tidy.g.n$Variable) 
# then assign each variable factor the appropriate number
tidy.g.n$Variable[tidy.g.n$Variable == "R.Nodes.1"] <- "1" 
tidy.g.n$Variable[tidy.g.n$Variable == "R.Nodes.2"] <- "2" 
tidy.g.n$Variable[tidy.g.n$Variable == "R.Nodes.3"] <- "3"
tidy.g.n$Variable[tidy.g.n$Variable == "R.Nodes.4"] <- "4"
tidy.g.n$Variable[tidy.g.n$Variable == "R.Nodes.5"] <- "5"
# now change back to a numerical variable
tidy.g.n$Variable <- as.numeric(tidy.g.n$Variable) 
# I think the data is ready to be analysed 
# rename columns to make more sense in the model 
colnames(tidy.g.n)[colnames(tidy.g.n) == 'Name'] <- 'Host'
colnames(tidy.g.n)[colnames(tidy.g.n) == 'Variable'] <- 'Time'
colnames(tidy.g.n)[colnames(tidy.g.n) == 'Value'] <- 'Nodes'
tidy.g.n$Nodes <- as.numeric(tidy.g.n$Nodes)
levels(tidy.g.n$Family)[11] <- "Dennstaedtiaceae"

# Only time points 2-4 for modelling purposes as there are too many zeroes otherwise
tidy.g.n1<- tidy.g.n[tidy.g.n$Time %in% c("2", "3", "4"),]
# relevel
tidy.g.n1$Host.No.Host <- relevel(tidy.g.n1$Host.No.Host, ref = "Host")
tidy.g.n1$Time <- as.factor(tidy.g.n1$Time)
tidy.g.n1$Time <- relevel(tidy.g.n1$Time, ref = "2")
tidy.g.n1$AnnPer <- relevel(tidy.g.n1$AnnPer, ref = "Per")
tidy.g.n1$Functional.group <- relevel(tidy.g.n1$Functional.group, ref = "Grass")
# Nodes cumulative at final time point

# siphon off time point two as this is the first time point with any reproductive nodes
attime2<-filter(tidy.g.n1, tidy.g.n1$Time == "2")
# sum nodes from time point 2, 3, 4, 5
cumnodes.1<-filter(tidy.g.n, tidy.g.n$Time == "2")[10]+
  filter(tidy.g.n, tidy.g.n$Time == "3")[10]+
  filter(tidy.g.n, tidy.g.n$Time == "4")[10]+
  filter(tidy.g.n, tidy.g.n$Time == "5")[10]
# create the data frame
cumnodes.1<- unlist(cumnodes.1)
cumnodes.1 <- data.frame(cumnodes.1)
colnames(cumnodes.1)[1] <- "Cumulative.Nodes" 
attime2$Cumulative.Nodes <- cumnodes.1$Cumulative.Nodes
attime2$Cumulative.Nodes <- as.numeric(attime2$Cumulative.Nodes)
# merge the two data frames
attime2.1<-merge(attime2, allgrowth[, c(1, 36)], by = "Unique_ID")

# for phylogenetic effects
attime2.1$animal <- attime2.1$Host

tidy.g.n2 <- tidy.g.n1[tidy.g.n1$Unique_ID != "1323",]
attime2.2<- merge(attime2.1, unique(survivaldata.2[,c(2,3)]))
levels(attime2.3$animal)[44] <- "Ulex europaeus"
levels(attime2.3$animal)[10] <- "Cystopteris dickieana"

attime2.3 <- attime2.2[attime2.2$Host != "No host",]

attime2.3$animal <- factor(attime2.3$animal)

# some summary stats

setDT(attime2.3)[, .(Mean.Nodes = mean(Cumulative.Nodes)), by = "Host"][order(-Mean.Nodes)]
# first 10 hosts
first11<-attime2.3[, .(Mean.Nodes = mean(Cumulative.Nodes)), by = "Host"][order(-Mean.Nodes)]$Host[1:11]

# model #
prior.3<-list(R=list(V=diag(1), nu=0.002),
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                     G2=list(V=diag(1), n=1, alpha.mu=rep(0, 1),alpha.V=diag(1)*1000)))

mcmcfix3<-MCMCglmm(Cumulative.Nodes ~ AnnPer + Functional.group + Transplant.Date,
                   random = ~Host + animal,
                   ginverse = list(animal=AinvULT2),
                   family="poisson", 
                   prior=prior.3, 
                   data=attime2.3[attime2.3$Host != "No host",],
                   nitt=13000*20,
                   thin=10*20,
                   burnin=3000*20,
                   verbose = T,
                   pr=T)
summary(mcmcfix3)
# AnnPer
write.csv(x = aod::wald.test(cov(mcmcfix3$Sol[,2, drop=F]), colMeans(mcmcfix3$Sol[,2, drop=F]), Terms=1)$result$chi2,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Reproductive_nodes_end/AnnPer_Wald_Test.csv")
# Functional group
write.csv(x = aod::wald.test(cov(mcmcfix3$Sol[,3:6, drop=F]), colMeans(mcmcfix3$Sol[,3:6, drop=F]), Terms=1:4)$result$chi2,
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Reproductive_nodes_end/Functional_Group_Wald_Test.csv")


# variance explained output
MCMCRepnorm2(mcmcfix3)
# broadly similar on data scale. Make explicit that they were calculated on the latent scale.
MCMCReppois(mcmcfix3, "units")

# joint phylogenetic distribution
# 90.7%
write.csv(list(posterior.mode(as.mcmc(rowSums(mcmcfix3$VCV[,c(1,2)])/rowSums(mcmcfix3$VCV))),
               # 66.5% - 95.6%
               HPDinterval(as.mcmc(rowSums(mcmcfix3$VCV[,c(1,2)])/rowSums(mcmcfix3$VCV)))),
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Reproductive_nodes_end/Joint_Phylogenetic_Variance.csv")

# and the phylogenetic component
# 74.5%
write.csv(list(posterior.mode(as.mcmc((mcmcfix3$VCV[,c(2)])/rowSums(mcmcfix3$VCV))),
               # 14.8% - 96.1%
               HPDinterval(as.mcmc((mcmcfix3$VCV[,c(2)])/rowSums(mcmcfix3$VCV)))),
          file = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2017_8/Manuscript/Models_Figures/Models/Single_Euphrasia_sp/Reproductive_nodes_end/Phylogenetic_Variance.csv")


# host variance
posterior.mode(as.mcmc((mcmcfix3$VCV[,c(1)])/rowSums(mcmcfix3$VCV)))
HPDinterval(as.mcmc((mcmcfix3$VCV[,c(1)])/rowSums(mcmcfix3$VCV)))

##### reproductive nodes over time #####

# leave phylogeny out of this one?
tidy.g.n2$animal <- tidy.g.n2$Host

levels(tidy.g.n2$animal)[10] <- "Cystopteris dickieana"
levels(tidy.g.n2$animal)[44] <- "Ulex europaeus"
levels(tidy.g.n2$Host)[10] <- "Cystopteris dickieana"
levels(tidy.g.n2$Host)[44] <- "Ulex europaeus"


prior.4<-list(R=list(V=diag(3), nu=0.002), 
              G=list(G1=list(V=diag(3), nu=3, alpha.mu=rep(0,3), alpha.V=diag(3)*1000),
                     G2= list(V=diag(1), nu=3, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))
# run for 13,000 X 20 iterations...
mcmcfix4<-MCMCglmm(round(Nodes)~Host.No.Host+Time*AnnPer+Functional.group,
                   random = ~ us(Time):Host + animal,
                   #ginverse = list(animal=AinvULT2),
                   rcov=~idh(Time):units,
                   family="poisson", 
                   prior=prior.4, 
                   data=tidy.g.n2[tidy.g.n2$Host != "No host",],
                   nitt=13000*20,
                   thin=10*20,
                   burnin=3000*20,
                   verbose = T,
                   pr=T)
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

exp(-4.51464+3.36769+1.20091-2.39262)

model.matrix(mcmcfix4)

##### **PLOTS** #####

## 17.6.19

ggplot(survivaldata.1, aes(x=Time, y=y,group=Name)) + 
  geom_point(alpha=0.1, position = position_jitter(width = 0.2, height = 0.2))+
  facet_wrap(~Family, scales = "free_x")+
  geom_smooth(aes(x=Time, y=y, col = Name), 
              method = "glm", 
              method.args = list(family = "binomial"), 
              se = F, size = 1)+
  ylim(c(0,1))+
  theme_bw()+
  ylab(label = "Probability of death")+
  labs(colour = "Species")+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20))

ggplot(survivaldata.2[Functional.group %in% c("Grass", "Fern", "Forb")], 
       aes(x=Time, y=y,group=Name))+
  geom_point(alpha=0.1, position = position_jitter(width = 0.2, height = 0.2))+
  facet_wrap(~Functional.group, scales = "free_x")+
  geom_smooth(aes(x=Time, y=y, col = Name), 
              method = "glm", 
              method.args = list(family = "binomial"), 
              se = F, size = 1)+
  ylim(c(0,1))+
  theme_bw()+
  ylab(label = "Probability of death")+
  labs(colour = "Species")+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20))


# make a binomial glm to add predictions to graph, in green short dashes
eha.4 <- glm(y ~ Name + Time, data = survivaldata.1,
             family = "binomial")
new.data <- data.frame(survivaldata.1)
new.data <- data.frame(Name = as.factor(sort(rep(unique(survivaldata.1$Name), 100))),
                       Time = as.numeric(rep(seq(1, 4, length = 100), 46, length=100)))
new.data$pred.full <- predict(eha.4, newdata = new.data, type = 'response')

# make a random effect binomial glm, add to graph in blue long dashes

eha.5 <- glmer(y ~  Time + (1 | Name), data = survivaldata.2,
               family = "binomial")
summary(eha.5)

new.data3 <- dplyr::full_join(new.data, unique(tidy.g.n2[,c(3,4, 7)]), by = c("Name" = "Host"))
new.data3$pred.full2 <- predict(eha.5, newdata = new.data, type = 'response')

unique(survivaldata.1$Name)

familykey <- data.table(Name = unique(survivaldata.1$Name),
                        Family = c("Fabaceae", "Poaceae", "Rubiaceae", "Asteraceae",
                                   "Amaranthaceae", "Valerianaceae", "Rubiaceae", "Poaceae", NA, "Poaceae",
                                   "Cistaceae", "Fabaceae", "Fabaceae", "Caryophyllaceae", "Lamiaceae",
                                   "Fabaceae", "Poaceae", "Polygonaceae", "Plantaginaceae", "Poaceae", 
                                   "Papaveraceae", "Cystopteridaceae", "Poaceae", "Brassicaceae", "Fabaceae",
                                   "Asteraceae", "Rosaceae", "Equisetaceae", "Asparagaceae", "Poaceae", "Amaranthaceae",
                                   "Amaryllidaceae", "Amaryllidaceae", "Orchidaceae", "Rosaceae", "Phrymaceae", "Apiaceae",
                                   "Pinaceae", "Apiaceae", "Caryophyllaceae", "Asteraceae", "Fabaceae", "Ericaceae", "Asteraceae",
                                   "Poaceae", "Dennstaedtiaceae"))
survivaldata.1 <- setDT(survivaldata.1)[familykey, on = "Name"]
# combine all three on a plot. Presentation. 12 x 8
ggplot(survivaldata.1[Family %in% c("Poaceae", "Fabaceae", "Asteraceae")], 
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


ggplot(survivaldata.1[Family %in% c("Poaceae", "Fabaceae")],
       aes(x = Time, y=y, group = Family))+
  geom_smooth(aes(x=Time, y=y, colour=Family), 
              method = "glm", 
              method.args = list(family = "binomial"), 
              se = T,
              size = 1, 
              alpha = 0.2)+
  theme_bw()+
  ylab(label = "Probability of death")+
  labs(colour = "Species")+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20))




ggplot(survivaldata.1[Family %in% c("Poaceae", "Fabaceae")], 
       aes(x=Time, y=y,group=Name))+
  geom_point(alpha=0.1, position = position_jitter(width = 0.2, height = 0.2))+
  facet_wrap(~Name, scales = "free_x")+
  geom_smooth(aes(x=Time, y=y), 
              method = "glm", 
              method.args = list(family = "binomial"), 
              se = F,
              col = "red", size = 1)+
  ylim(c(0,1))+
  theme_bw()+
  geom_line(data = new.data,aes(y = pred.full, colour=Name),colour = "blue", lty =2)+
  geom_line(data = new.data3,aes(y = pred.full2, colour=Name),colour = "black", lty =3)+
  theme(strip.background  = element_blank(), strip.text.x = element_text(size = 10))

custplot<-filter(survivaldata.1, Name %in% c("Agrostis capillaris",
                                             "Cynosurus cristatus",
                                             "Festuca rubra",
                                             "Holcus lanatus",
                                             "Hordeum vulgare",
                                             "Lagurus ovatus",
                                             "Phleum pratense",
                                             "Zea mays",
                                             "Lathyrus japonicus",
                                             "Lotus corniculatus",
                                             "Ononis spinosa",
                                             "Trifolium pratense",
                                             "Ulex europaeus ",
                                             "Vicia cracca"))
custplot$Name <- factor(custplot$Name)
levels(custplot$Name) <- c("Agrostis capillaris",
                           "Cynosurus cristatus",
                           "Festuca rubra",
                           "Holcus lanatus",
                           "Hordeum vulgare",
                           "Lagurus ovatus",
                           "Phleum pratense",
                           "Zea mays",
                           "Lathyrus japonicus",
                           "Lotus corniculatus",
                           "Ononis spinosa",
                           "Trifolium pratense",
                           "Ulex europaeus ",
                           "Vicia cracca")
custplot<-custplot[sort.list(custplot$Name),]
custplot$Factor <- c(rep("Poaceae", 735), rep("Fabaceae", 520))

ggplot(custplot, 
       aes(x=Time, y=y, group=Name))+
  geom_point(alpha=0.1, position = position_jitter(width = 0.2, height = 0.2))+
  facet_wrap(~Factor, scales = "free_x")+
  geom_smooth(aes(x=Time, y=y, colour=Name), 
              method = "glm", 
              method.args = list(family = "binomial"), 
              se = F,
              size = 1)+
  ylim(c(0,1))+
  theme_bw()+
  
  #geom_line(data = new.data,aes(y = pred.full),colour = "blue", lty =2)+
  #geom_line(data = new.data3[new.data3$Name %in% unique(custplot$Name),],aes(y = pred.full2 ,colour =Name))+
  theme(strip.background  = element_blank(), strip.text.x = element_text(size = 10))

new.data2 <- new.data[new.data$Name %in% unique(custplot$Name),]
new.data2$Name <- factor(new.data2$Name)
levels(new.data2$Name) <- c("Agrostis capillaris",
                            "Cynosurus cristatus",
                            "Festuca rubra",
                            "Holcus lanatus",
                            "Hordeum vulgare",
                            "Lagurus ovatus",
                            "Phleum pratense",
                            "Zea mays",
                            "Lathyrus japonicus",
                            "Lotus corniculatus",
                            "Ononis spinosa",
                            "Trifolium pratense",
                            "Ulex europaeus ",
                            "Vicia cracca")
new.data2<-new.data2[sort.list(new.data2$Name),]
new.data2$Factor <- c(rep("Poaceae", 800), rep("Fabaceae", 600))


ggplot(new.data2,
       aes(y = pred.full ,colour =Name))+
  geom_smooth(aes(x=Time, y=pred.full, colour=Factor), 
              method = "glm", 
              method.args = list(family = "binomial"), 
              se = T,
              size = 1)+
  theme_bw()
facet_wrap(~Factor)


p.eha1$predict<-predict(eha.1, type = 'response', marginal = NULL, posterior = "all")

str(p.eha1)
p.eha1$predict <- as.vector(p.eha1$predict)
custplot<-filter(p.eha1, p.eha1$Name %in% c("Agrostis capillaris",
                                            "Cynosurus cristatus",
                                            "Festuca rubra",
                                            "Holcus lanatus",
                                            "Hordeum vulgare",
                                            "Lagurus ovatus",
                                            "Phleum pratense",
                                            "Zea mays",
                                            "Lathyrus japonicus",
                                            "Lotus corniculatus",
                                            "Ononis spinosa",
                                            "Trifolium pratense",
                                            "Ulex europaeus ",
                                            "Vicia cracca"))
custplot$Name <- factor(custplot$Name)
levels(custplot$Name) <- c("Agrostis capillaris",
                           "Cynosurus cristatus",
                           "Festuca rubra",
                           "Holcus lanatus",
                           "Hordeum vulgare",
                           "Lagurus ovatus",
                           "Phleum pratense",
                           "Zea mays",
                           "Lathyrus japonicus",
                           "Lotus corniculatus",
                           "Ononis spinosa",
                           "Trifolium pratense",
                           "Ulex europaeus ",
                           "Vicia cracca")
custplot<-custplot[sort.list(custplot$Name),]
custplot$Factor <- c(rep("Poaceae", 735), rep("Fabaceae", 519))

ggplot(custplot, 
       aes(x=Time, y=y, group=Name))+
  geom_point(alpha=0.1, position = position_jitter(width = 0.2, height = 0.2))+
  #facet_wrap(~Factor, scales = "free_x")+
  ylim(c(0,1))+
  theme_bw()+
  geom_smooth(aes(y = predict ,colour =Factor), method = "glm", method.args = list(family = "binomial"), se = F, alpha = 0.1)



# tree plot for cumulative reproductive nodes

# create the ggtree object
(p<-ggtree(phytools::force.ultrametric(constraint.2))+
    geom_tiplab(size = 4)+ theme(axis.line = element_blank())+
    geom_cladelabel(node=64, label="Legumes", align=T, angle = 270, offset =1, hjust = 0.5, offset.text = 0.05, barsize = 1.5)+
    geom_cladelabel(node=53, label="Grasses", align=T, angle = 270, offset =1, hjust = 0.5, offset.text = 0.05, barsize = 1.5))
# data to go with the tree
setDT(attime2.1)
phylodat<-attime2.1[,c("Host","Cumulative.Nodes", "Functional.group", "Family")]
phylodat$group <- phylodat$Host

phylodat1.1 <- phylodat[, .(Mean.Cum = mean(na.omit(Cumulative.Nodes)),
                            SD.Cum = sd(na.omit(Cumulative.Nodes))/sqrt(.N),
                            N = .N,
                            Functional.group = (Functional.group),
                            Family = (Family)), by = "Host"]
phylodat1.1 <- unique(phylodat1.1)
levels(phylodat1.1$Host)[10] <- "Cystopteris dickieana"
levels(phylodat1.1$Host)[44] <- "Ulex europaeus"
phylodat1.1$group = phylodat1.1$Host

# make tip labels pretty
p$data$label <- gsub(pattern = "_", replacement = " ", x = p$data$label)
levels(phylodat$Host)[10] <- "Cystopteris dickieana"
levels(phylodat$Host)[44] <- "Ulex europaeus"

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
facet_plot(finalplot, panel = "Cumulative Reproductive Nodes",
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


# reproductive nodes over time (for functional group and annual or perennial!)
setDT(tidy.g.n)[, .(Nodes = mean(Nodes), uCI = mean(Nodes)+sd(Nodes)/sqrt(.N), lCI = mean(Nodes)-sd(Nodes)/sqrt(.N)), by = c("AnnPer", "Time")] %>%
  ggplot(aes(x = Time, y= Nodes))+ 
  facet_wrap(~AnnPer)+
  geom_point()+
  geom_errorbar(aes(ymax = uCI, ymin =lCI), width = 0.05)+
  geom_line()+
  theme_bw()+
  theme(strip.background = element_blank(), strip.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

colnames(familykey)[1] <- "Host"
# add "Unique_ID" in by = if you want the fancy graphs, and group by Unique_ID
# need an 'average host'
avghost<-setDT(tidy.g.n)[, .(Nodes = mean(Nodes), uCI = mean(Nodes)+sd(Nodes)/sqrt(.N), lCI = mean(Nodes)-sd(Nodes)/sqrt(.N)), by = c("Time")]
avghost2 <- data.table(Time = avghost$Time, 
                       Host = as.factor(rep("Average Host", 5)),
                       Nodes = avghost$Nodes,
                       uCI = avghost$uCI,
                       lCI = avghost$lCI,
                       Family = as.factor(rep("Average", 5)))

overtime <- setDT(tidy.g.n)[, .(Nodes = mean(Nodes), uCI = mean(Nodes)+sd(Nodes)/sqrt(.N), lCI = mean(Nodes)-sd(Nodes)/sqrt(.N)), by = c("Time", "Host")][familykey, on = "Host"][Host %in% first11,] 
overtime <- rbind(overtime, avghost2)
overtime$Host <- factor(overtime$Host, levels = c(as.character(first11), "Average Host"))

overtime %>%
  ggplot(aes(x = Time, y= Nodes, group = Host))+ 
  facet_wrap(~Host)+
  geom_errorbar(aes(ymax = uCI, ymin =lCI), width = 0.05)+
  geom_line(aes(colour = Family))+
  geom_point(alpha=0.5)+
  theme_bw()+
  theme(strip.background = element_blank(), strip.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20))

# reproductive output distributions #

ggplot(attime2.1, aes(x = log(Cumulative.Nodes+1)))+
  geom_density(aes(fill = Host), alpha = 0.1)+
  #xlim(c(0, 10))+
  theme(legend.position = "none")

## corolla ##
corollaplot<-filter(floweringsofar, !is.na(floweringsofar$Corolla.Length)) %>%
  group_by(Host.Species, Family, Functional.group) %>%
  select(c(1:5, "Corolla.Length")) %>%
  summarise(Mean = mean(Corolla.Length),
            SE = sd(Corolla.Length)/sqrt(n()),
            Replicates = n())

ggplot(corollaplot, aes(reorder(x = Host.Species, Mean), y = Mean))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1, angle = 60))+
  geom_errorbar(aes(ymin = Mean - SE, ymax= Mean +SE), width=0.5)+
  geom_point(aes(colour = Functional.group), size=3)+
  ylab(label = "Corolla length (mm)")+
  xlab(label = "Species")+
  scale_colour_viridis_d()

# create the ggtree object
r<-ggtree(phytools::force.ultrametric(constraint.1))+
  geom_tiplab(size = 4)+ theme(axis.line = element_blank())
# data to go with the tree
setDT(floweringsofar)
phylodat2<-floweringsofar[Host.Species != "No host", c("Host.Species","Corolla.Length")]
phylodat2$group <- phylodat2$Host.Species

phylodat2.1 <- phylodat2[, .(Mean.Corolla = mean(na.omit(Corolla.Length)),
                             SD.Corolla = sd(na.omit(Corolla.Length))/sqrt(.N),
                             N = .N), by = "Host.Species"]
phylodat2.1$group <- phylodat2.1$Host.Species

# for making the boxplots
(finalplot2<-facet_plot(r + xlim_tree(1.6), panel = "Corolla Length",
                        data = phylodat2.1,
                        mapping = aes(x=Mean.Corolla, 
                                      xmin = Mean.Corolla-SD.Corolla,
                                      xmax = Mean.Corolla+SD.Corolla,
                                      group = group),
                        geom = geom_errorbarh)+
    theme_tree2())
# overlay the red observational points
(finalplot2.1 <- facet_plot(finalplot2, panel = "Corolla Length",
                            data = phylodat2.1,
                            mapping = aes(x=Mean.Corolla, group = group),
                            geom = geom_point,
                            col="red",
                            size = 3)+
    theme(strip.text.x = element_text(size=20, face = "bold", vjust = -1),
          strip.background = element_rect(colour="white", fill="white"),
          axis.line.x = element_line(colour = "black")))

# Days to flower

# create the ggtree object: geom_text2(aes(label=node)) for labels of nodes.
(r<-ggtree(phytools::force.ultrametric(constraint.1))+
    geom_tiplab(size = 4)+ theme(axis.line = element_blank())+
    geom_cladelabel(node=60, label="Legumes", align=T, angle = 270, offset =1, hjust = 0.5, offset.text = 0.05, barsize = 1.5)+
    geom_cladelabel(node=50, label="Grasses", align=T, angle = 270, offset =1, hjust = 0.5, offset.text = 0.05, barsize = 1.5))

# maybe flip later... flip(tree_view = r, node1 =58, node2 = 65)
# data to go with the tree
setDT(floweringsofar)
phylodat3<-floweringsofar[Host.Species != "No host", c("Host.Species","Days.since.germination", "Functional.group","Family")]
phylodat3$group <- phylodat3$Host.Species

phylodat3.1 <- phylodat3[, .(Mean.Days = mean(na.omit(Days.since.germination)),
                             SD.Days = sd(na.omit(Days.since.germination))/sqrt(.N),
                             N = .N,
                             Functional.group = (Functional.group),
                             Family = (Family)), by = "Host.Species"]
phylodat3.1 <- unique(phylodat3.1)
phylodat3.1$group <- phylodat3.1$Host.Species


# for making the boxplots
(finalplot3<-facet_plot(r + xlim_tree(1.6), panel = "Days to Flower",
                        data = phylodat3.1,
                        mapping = aes(x=Mean.Days,
                                      xmin=Mean.Days-SD.Days,
                                      xmax=Mean.Days+SD.Days,
                                      group = group),
                        geom = geom_errorbarh)+
    theme_tree2())
# overlay the red observational points
(finalplot3.1 <- facet_plot(finalplot3, panel = "Days to Flower",
                            data = phylodat3.1,
                            mapping = aes(x=Mean.Days, group=group),
                            geom = geom_point,
                            col="red",
                            size = 3)+
    theme(strip.text.x = element_text(size=20, face = "bold", vjust = -1),
          strip.background = element_rect(colour="white", fill="white"),
          axis.line.x = element_line(colour = "black")))

#combined
finalplot3.2 <- facet_plot(finalplot3.1, panel = "Corolla Length (mm)",
                           data = phylodat2.1,
                           mapping = aes(x=Mean.Corolla, 
                                         xmin = Mean.Corolla-SD.Corolla,
                                         xmax = Mean.Corolla+SD.Corolla,
                                         group = group),
                           geom = geom_errorbarh)

(finalplot4<-facet_plot(finalplot3.2, panel = "Corolla Length (mm)", data = phylodat2.1, geom = geom_point, 
                        aes(x=Mean.Corolla, group = group), col="red",
                        size = 3)+theme_bw()+  theme(strip.text.x = element_text(size=20, face = "bold", vjust =1),
                                                     strip.background = element_rect(colour="white", fill="white"),
                                                     axis.line.x = element_line(colour = "black"))+
    geom_hline(yintercept = c(6.5, 14.5, 19.5, 25.5), lty =2, alpha = 0.3))

orderoftaxa<- as.data.table(ggplot_build(finalplot4)$data[7][1])$group

# survival

# create the ggtree object
s<-ggtree(phytools::force.ultrametric(constraint.2))+
  geom_tiplab(size = 4)+ theme(axis.line = element_blank())
# data to go with the tree
setDT(survivaldata.2)
phylodat4<-survivaldata[Name != "No host", c("Name","y")]


phylodat4$y2 <- (phylodat4$y-1)*-1

phylodat4 <- phylodat4[, .(Mean.survival = mean(y2),
                           SD.survival = sd(y2)/sqrt(.N)), by = "Name"]
phylodat4$group <- phylodat4$Name

(finalplot4<-facet_plot(s + xlim_tree(1.6), panel = "Survival",
                        data = phylodat4,
                        mapping = aes(x =Mean.survival,
                                      xmin=Mean.survival-SD.survival, 
                                      xmax=Mean.survival+SD.survival,
                                      group = group),
                        geom = geom_errorbarh)
  
  facet_plot(finalplot4, panel = "Survival",
             data = phylodat4,
             mapping = aes(x=Mean.survival, group = group),
             geom = geom_point,
             size =3,
             col = "red")+
    theme_tree2()+
    theme(strip.text.x = element_text(size=20, face = "bold", vjust = -1),
          strip.background = element_rect(colour="white", fill="white"),
          axis.line.x = element_line(colour = "black"))
  
  # ridge plot for the phylogenetic constribution to the four main traits. 
  c1 <- data.table(Density =(mcmcfix1$VCV[,c(2)])/rowSums(mcmcfix1$VCV),
                   Trait = as.factor(rep("Corolla Length (Gaussian)", 1000)))
  d1 <- data.table(Density =(mcmcfix2$VCV[,c(2)])/rowSums(mcmcfix2$VCV),
                   Trait = as.factor(rep("Days to Flower (Poisson)", 1000)))
  s1 <- data.table(Density =(eha.1$VCV[,c(2)])/rowSums(eha.1$VCV),
                   Trait = as.factor(rep("Survival (Threshold, logit)", 1000)))
  r1 <- data.table(Density =(mcmcfix3$VCV[,c(2)])/rowSums(mcmcfix3$VCV),
                   Trait = as.factor(rep("Total Reproductive Output (Poisson)", 1000)))
  final <- rbind(c1, d1, s1, r1)
  final$Density <- as.numeric(final$Density)
  
  final$Trait <- factor(x = final$Trait, levels = c("Days to Flower (Poisson)",
                                                    "Corolla Length (Gaussian)",
                                                    "Total Reproductive Output (Poisson)",
                                                    "Survival (Threshold, logit)"))
  
  ggplot(final, aes(x = Density, y = Trait))+
    theme_bw()+geom_density_ridges()+geom_vline(xintercept = 0, col = "red", lty = 2, size = 1)+
    scale_x_continuous(limits = c(-0.2,1))
  
  
##### Correlations #####
  
  # correlations needed: days to flower:corolla length, days to flower:survival,
  # days to flower:nodes, corolla length:survival, corolla length: nodes, nodes:survival
  setDT(survivaldata.2); setDT(attime2.3); setDT(floweringsofar)
  # make a data table of survival, nodes, corolla length and days to flower
  a <- survivaldata.2[, .(Mean.survival = mean(y2)), by = "animal"][attime2.3[, .(Cumulative.nodes = mean(Cumulative.Nodes)), by = "animal"], on = "animal"]
  corrr <- floweringsofar[, .(Mean.flowering = mean(Days.since.germination), Mean.corolla = mean(Corolla.Length)), by = "animal"][a , on = "animal"]
  
  #corrr$NAs <- apply(corrr, 1, function(x) sum(as.numeric(is.na(x))))
  #corrr <- corrr[NAs == 0,]
  # correlations
  cor(corrr[, -c("animal")], method = "kendall", use = "pairwise.complete.obs")
  
  # only positive correlations are cumulative nodes and corolla length and
  corrr2 <- unique(survivaldata.2[,c("Family", "animal")])[corrr, on = "animal"]
  
  corrr2$subFamily <- ifelse(test = corrr2$Family %in% c("Poaceae", "Fabaceae", "Asteraceae"),
                             yes = corrr2$Family, no = "Other families")
  
  ggplot(corrr2, aes(y = Cumulative.nodes, x = Mean.survival))+geom_point(aes(colour = subFamily), size = 3)+
    theme_bw() + ggrepel::geom_text_repel(aes(label = animal))
  
  # make a data table of survival, nodes, corolla length and days to flower for Unique_ID
  a <- survivaldata.2[, .(Mean.survival = mean(y)), by = "Unique_ID"][attime2.3[, .(Cumulative.nodes = mean(Cumulative.Nodes)), by = "Unique_ID"], on = "Unique_ID"]
  corrr <- floweringsofar[, .(Mean.flowering = mean(Days.since.germination), Mean.corolla = mean(Corolla.Length)), by = "Unique_ID"][a , on = "Unique_ID"]
  
  #corrr$NAs <- apply(corrr, 1, function(x) sum(as.numeric(is.na(x))))
  #corrr <- corrr[NAs == 0,]
  # correlations
  cor(corrr[, -c("Unique_ID")], method = "kendall", use = "pairwise.complete.obs")
  
  # only positive correlations are cumulative nodes and corolla length and
  corrr2 <- unique(survivaldata.2[,c("Family", "Unique_ID")])[corrr, on = "Unique_ID"]
  corrr2$subFamily <- ifelse(test = corrr2$Family %in% c("Poaceae", "Fabaceae", "Asteraceae"),
                             yes = corrr2$Family, no = "Other families")
  
  ggplot(corrr2, aes(x = Mean.survival, y = Cumulative.nodes))+geom_point(aes(colour = subFamily), size = 3)+
    theme_bw()
  
  
  
  
  
##### summary stats #####
  
  allgrowth[, .(Mean_nodes = mean(C.R.Nodes.F),
                SE_Nodes = sd(C.R.Nodes.F)/sqrt(.N)), by = "Functional.group"][order(-Mean_nodes)]
  