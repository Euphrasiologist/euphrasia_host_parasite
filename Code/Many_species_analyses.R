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

# "not" in
`%!in%` <- function(x,y)!('%in%'(x,y))

# little function to convert poisson slopes and intercepts.
ablines <- function(intercept, slope, x1, x2){
  x <- seq(x1,x2,0.01)
  res <- exp(intercept + slope*x)
  res
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

# replicates

sp_reps <- rnodes3[, .(.N), by = .(Host_code, Euphrasia_sp2)][order(Euphrasia_sp2)]
sp_reps <- data.frame(lapply(sp_reps, gsub, pattern = "PLA", replacement = "Plantago lanceolata", fixed = TRUE))
sp_reps<- data.frame(lapply(sp_reps, gsub, pattern = "ACU", replacement = "Agrostis curtisii", fixed = TRUE))
sp_reps<- data.frame(lapply(sp_reps, gsub, pattern = "HLA", replacement = "Holcus lanatus", fixed = TRUE))
sp_reps<- data.frame(lapply(sp_reps, gsub, pattern = "LPE", replacement = "Lolium perenne", fixed = TRUE))
sp_reps<- data.frame(lapply(sp_reps, gsub, pattern = "FOV", replacement = "Festuca ovina", fixed = TRUE))
sp_reps<- data.frame(lapply(sp_reps, gsub, pattern = "OVU", replacement = "Origanum vulgare", fixed = TRUE))
sp_reps<- data.frame(lapply(sp_reps, gsub, pattern = "DFL", replacement = "Deschampsia flexuosa", fixed = TRUE))
sp_reps<- data.frame(lapply(sp_reps, gsub, pattern = "VCH", replacement = "Veronica chamaedrys", fixed = TRUE))
sp_reps<- data.frame(lapply(sp_reps, gsub, pattern = "PMA", replacement = "Plantago maritima", fixed = TRUE))
sp_reps<- data.frame(lapply(sp_reps, gsub, pattern = "UGA", replacement = "Ulex gallii", fixed = TRUE))
sp_reps<- data.frame(lapply(sp_reps, gsub, pattern = "LCO", replacement = "Lotus corniculatus", fixed = TRUE))
sp_reps<- data.frame(lapply(sp_reps, gsub, pattern = "CVU", replacement = "Calluna vulgaris", fixed = TRUE))
sp_reps<- data.frame(lapply(sp_reps, gsub, pattern = "HPU", replacement = "Hypericum pulchrum", fixed = TRUE))

fwrite(x =sp_reps, file = "./Data/Many_species/Replicates.csv")

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



##### Part 2: Add data on reproductive nodes at the end of the season, plus more tidying #####

# read in nodes data
rnodes <- fread("./Data/Many_species/REPRODUCTIVENODES.csv")
# merge with data generated at first flowering. Merged so that UniqueID's that survived to end of season are preserved 
# (not all that first flowered...)
rnodes2 <- expt2dat.2[rnodes, on = "Unique_ID"]
rnodes2$Euphrasia_sp2 <- factor(rnodes2$Euphrasia_sp2)

# create the interaction factors
# between species of Euphrasia and host
rnodes2$Host_given_Euphrasia <- interaction(rnodes2$Euphrasia_sp2, rnodes2$Host_code)

# remove any NA's
rnodes3 <- rnodes2[complete.cases(Euphrasia_sp2)]
# clean up population for no umambiguity
rnodes3 <- rnodes3[, Population := Population2][, -"Population2"]
# change to factors
for_factor <- c("Euphrasia_sp2", "Host_code", "Population", "Host_given_Euphrasia")
for(col in for_factor){
  set(rnodes3, j=col, value=as.factor(rnodes3[[col]]))
}

# tidy the dates
rnodes3[, Trans_date := as.Date(rnodes3$Trans_date, format = "%d/%m/%Y")]
rnodes3[, Germ_date := as.Date(rnodes3$Germ_date, format = "%d/%m/%Y")]
# normalised transplant date
rnodes3[, Norm_Trans_date := as.POSIXlt(rnodes3$Trans_date)$yday-108]
# relevel population 

rnodes3$Population <- relevel(rnodes3$Population, ref = "A1766")

##### Part 3: Model of reproductive nodes as a function of Euphrasia species and look for interaction #####

# prior for interaction model.
prior.manysp <- list(R=list(V=diag(1), nu=0.002), 
                     G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                            #G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                            G3=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))


manysp.2<-MCMCglmm(Reproductive_nodes ~ Euphrasia_sp2 + Population + Norm_Trans_date,
                   random = ~ Host_code + Host_given_Euphrasia,
                   family = "poisson",
                   prior=prior.manysp,
                   data = rnodes3,
                   nitt = 13000*10,
                   burnin = 3000*10,
                   thin = 10*10,
                   pr=TRUE,
                   verbose = TRUE)
summary(manysp.2)

# correlation in host effects
posterior.mode(manysp.2$VCV[,"Host_code"]/(manysp.2$VCV[,"Host_code"]+manysp.2$VCV[,"Host_given_Euphrasia"]))
HPDinterval(manysp.2$VCV[,"Host_code"]/(manysp.2$VCV[,"Host_code"]+manysp.2$VCV[,"Host_given_Euphrasia"]))
# and the inverse for completeness
posterior.mode(manysp.2$VCV[,"Host_given_Euphrasia"]/(manysp.2$VCV[,"Host_code"]+manysp.2$VCV[,"Host_given_Euphrasia"]))
HPDinterval(manysp.2$VCV[,"Host_given_Euphrasia"]/(manysp.2$VCV[,"Host_code"]+manysp.2$VCV[,"Host_given_Euphrasia"]))

# wald tests of significance
# Euphrasia species & population?
write.csv(x = aod::wald.test(cov(manysp.2$Sol[,2:6, drop=F]), colMeans(manysp.2$Sol[,2:6, drop=F]), Terms=1:5)$result$chi2,
          file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/Euphrasia_pop_sig.csv")


write.csv(x = specify_decimal(summary(manysp.2)$solutions, 4),
          file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/Model_solutions.csv")

write.csv(x = data.table(HOST = MCMCReppois(mod = manysp.2, y = "Host_code"),
                         HOST_GIVEN_EUPHRASIA = MCMCReppois(mod = manysp.2, y = "Host_given_Euphrasia"),
                         UNITS = MCMCReppois(mod = manysp.2, y = "units"), keep.rownames = TRUE),
          file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/Variance_Components.csv")

# visualise the variance components
VCVdensity(manysp.2)+xlim(0, 1)

# so there are interactions here, let's dig further
# so this is great, we have the posterior modes for population:host interactions and host

sol_int <- Solapply(manysp.2)[order(Grouped_Value)][Group %in% c("Host_code", "Host_given_Euphrasia")][, .(Variable = Variable, Grouped_Value = exp(Grouped_Value))]

sol_int[, Grouped_Value := specify_decimal(Grouped_Value, 4)]
write.csv(x = sol_int,
          file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/Posterior_Modes.csv")


# significance of random effects
# NOTE: transplant date not added, as non significant in main model
# need to add Obs
rnodes3$Obs <- as.factor(1:nrow(rnodes3)) 

# full model
manysp.2LR1 <- glmer(Reproductive_nodes ~ Euphrasia_sp2 + Population + (1 | Host_code) + (1 | Host_given_Euphrasia) + (1|Obs),
                     family = "poisson", data = rnodes3)
# is Host_code significant?
manysp.2LR2 <- glmer(Reproductive_nodes ~ Euphrasia_sp2 + Population + (1 | Host_given_Euphrasia) + (1 | Obs),
                     family = "poisson", data = rnodes3)
# is Host_given_Euphrasia significant?
manysp.2LR3 <- glmer(Reproductive_nodes ~ Euphrasia_sp2 + Population + (1 | Host_code) + (1 | Obs),
                     family = "poisson", data = rnodes3)
# is Population significant?
#manysp.2LR4 <- glmer(Reproductive_nodes ~ Euphrasia_sp2 + Population + (1 | Host_code) +  (1 | Host_given_Euphrasia) + (1 | Obs),
 #                    family = "poisson", data = rnodes3)

write.csv(
  x = data.table(
    # Host_code is significant
    HOST_CODE = anova(manysp.2LR1, manysp.2LR2),
    # Host_given_code_Euphrasia is significant
    HOST_GIVEN_EUPHRASIA = anova(manysp.2LR1, manysp.2LR3),
    # Population is highly significant
    #HOST_GIVEN_POPULATION = anova(manysp.2LR1, manysp.2LR4), 
    keep.rownames = TRUE
  ), file = "./Data/Many_species/Model_outputs/Host_parasite_interaction/LRTs_of_models.csv"
)





##### Plot 2: Posterior Modes of the interaction model #####

# visualise the posterior modes better.

pms <- Solapply(manysp.2, coda::HPDinterval)[order(`Posterior Mode`)][Group %in% c("Host_code", "Host_given_Euphrasia")]
pms2 <- pms[Group == "Host_given_Euphrasia"]

lis <- strsplit(x = pms2$Variable, split = ".", fixed = TRUE)
pms2[, Euphrasia_species := sapply(lis, "[[", 2)]
pms2[, Host := sapply(lis, "[[", 3)]

pms3 <- pms2[,.(Host, Euphrasia_species, lowerHPD, upperHPD, `Posterior Mode`)]

ggplot(pms3, aes(x = reorder(Euphrasia_species, `Posterior Mode`), y = `Posterior Mode`))+
  facet_wrap(~Host)+
  geom_hline(yintercept = 0, colour = "red", lty = 2, size = 2)+
  geom_errorbar(aes(ymin = lowerHPD, ymax = upperHPD, width = 0.4))+
  geom_point(size = 4)+
  theme_minimal()+
  theme(strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 20),
        axis.text.x = element_text(angle = 60, hjust = 1))

  #scale_x_discrete(limits = c("ACU", "DFL", "FOV", "HLA", "LPE"))



plot_2.1 <- pms %>%
  ggplot(aes(x = reorder(Variable, `Posterior Mode`), y = `Posterior Mode`))+
  geom_point()+
  geom_errorbar(aes(ymin = lowerHPD, ymax = upperHPD))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(filename = "./Figures/Many_species/posterior_interaction_modes_mod_1.pdf", plot = plot_2.1, 
       device = "pdf", width = 21, height = 6, units = "in")

m1 <- data.table(Density =(manysp.2$VCV[,c(1)])/rowSums(manysp.2$VCV),
                 Trait = as.factor(rep("Host contribution", 1000)))
m2 <- data.table(Density =(manysp.2$VCV[,c(2)])/rowSums(manysp.2$VCV),
                 Trait = as.factor(rep("Host:Euphrasia interaction", 1000)))
m3 <- data.table(Density =(manysp.2$VCV[,c(3)])/rowSums(manysp.2$VCV),
                 Trait = as.factor(rep("Residual", 1000)))
final <- rbind(m1,m2,m3)
final$Density <- as.numeric(final$Density)

final$Trait <- factor(x = final$Trait, levels = rev(c("Host:Euphrasia interaction",
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

ggsave(filename = "./Figures/Many_species/posterior_interaction_dist.pdf", plot = plot_2.2, 
       device = "pdf", width = 10, height = 6, units = "in") 

##### Plot 3: Raw data for the manuscript, means and SE's #####
dat <- rnodes3[, .(mean = mean(log(Reproductive_nodes), na.rm = TRUE),
                   sem = sd(log(Reproductive_nodes), na.rm = TRUE)/sqrt(.N),
                   N = .N), by = c("Euphrasia_sp2", "Host_code")]
meandat <- rnodes3[,.(meanH = mean(log(Reproductive_nodes), na.rm = TRUE)), by = .(Host_code)][order(meanH)]

dat$Euphrasia_sp2 <- factor(dat$Euphrasia_sp2, levels = c("Euphrasia anglica", 
                                                          "Euphrasia vigursii",
                                                          "Euphrasia micrantha",
                                                          "Euphrasia tetraquetra"))
dat <- dat[meandat, on = .(Host_code)]
(plot_3.1<- dat %>% 
  
  ggplot(aes(x = reorder(Host_code, mean) , y = mean))+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem, group=Euphrasia_sp2), position = position_dodge(width = 0.9), width=0.6)+
  geom_point(aes(fill = Euphrasia_sp2), position = position_dodge(width = 0.9), size=4, pch=21, stroke=1)+
  facet_wrap(~Euphrasia_sp2, nrow = 1)+
    geom_text(data = data.frame(mean = rep(5,4),
                                Host_code = rep("HPU",4),
                                Euphrasia_sp2 = c("Euphrasia anglica",
                                                  "Euphrasia micrantha",
                                                  "Euphrasia tetraquetra",
                                                  "Euphrasia vigursii")), 
              label = letters[1:4], 
              nudge_x = 0.5, 
              nudge_y = -0.15, 
              size=8)+
  theme_bw()+ theme(strip.text.x = element_blank(),
                    strip.background = element_rect(colour="white", fill="white"),
                    axis.line.x = element_line(colour = "black"),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
                    axis.text.y = element_text(size = 20),
                    axis.title.x.bottom = element_text(size = 20),
                    axis.title.y.left = element_text(size = 20),
                    panel.border = element_rect(colour = "black", fill=NA, size=2),
                    panel.spacing.x = unit(0.5, "lines"),
                    legend.title = element_text(size = 14.7),
                    legend.text = element_text(face = "italic", size = 11), panel.spacing = unit(0, "lines"))+
  xlab(label = "Host Species")+
  #ylab(label = expression(paste(log[e], "(", italic("Euphrasia"), " nodes)")))+
  ylab(label = expression(paste(italic("Euphrasia"), " reproductive nodes")))+
  scale_y_continuous(breaks = c(0,0.6931472, 1.609438,2.302585,3.218876,4.60517), # 1,2, 5, 10, 50, 100
                       labels = round(c(exp(0),exp(0.6931472), exp(1.609438), exp(2.302585), exp(3.218876), exp(4.60517)), 0)) +
  scale_fill_manual(name = "\n      Host Species", 
                    limits = c("Euphrasia anglica", "Euphrasia micrantha", "Euphrasia tetraquetra", "Euphrasia vigursii"),
                    labels = c("HPU = Hypericum pulchrum",
                               "CVU = Calluna vulgaris",
                               "OVU = Origanum vulgare",
                               "HLA = Holcus lanatus\nUGA = Ulex gallii\nPMA = Plantago maritima\nFOV = Festuca ovina\nPLA = Plantago lanceolata\nDFL = Deschampsia flexuosa\nVCH = Veronica chamaedrys\nACU = Agrostis curtisii\nLPE = Lolium perenne\nLCO = Lotus corniculatus"),
                    values = c(cbPalette[2], cbPalette[3], cbPalette[4], cbPalette[5]), 
                    guide = guide_legend(override.aes = list(color = "white", size = 0.001, fill = "white")))
  )

ggsave(filename = "./Figures/Many_species/population_cum_nodes.pdf", plot = plot_3.1, 
       device = "pdf", width = 10, height = 6, units = "in", useDingbats=FALSE)
ggsave(filename = "./Figures/Many_species/population_cum_nodes.jpeg", plot = plot_3.1, 
       device = "jpeg", width = 10, height = 6, units = "in")

dat$Host_code <- factor(dat$Host_code, levels = unique(dat$Host_code[order(dat$meanH)]))

levels(dat$Euphrasia_sp2) <- c("Euphrasia anglica", "Euphrasia micrantha", "Euphrasia tetraquetra", "Euphrasia vigursii")

(plot_3.1.1<- ggplot(dat, aes(x = meanH , y = mean))+
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem, group=Euphrasia_sp2), position = position_dodge(width = 0.9), width=0.1)+
    ggrepel::geom_text_repel(aes(label = Host_code), box.padding = 1.5, max.overlaps = 30)+
    geom_point(aes(fill=Euphrasia_sp2), position = position_dodge(width = 0.9), size=4, pch=21, stroke=1)+
    facet_wrap(~Euphrasia_sp2, nrow = 1)+
    geom_abline(slope=1, intercept=0, col="red", lty=2)+
    geom_text(data = data.frame(mean = rep(5,4),meanH = rep(0,4), Euphrasia_sp2 = c("Euphrasia anglica", 
                                                                                    "Euphrasia micrantha",
                                                                                    "Euphrasia tetraquetra",
                                                                                    "Euphrasia vigursii")), label = letters[5:8], nudge_x = -0.1, 
              nudge_y = 0, size=8)+
    theme_bw()+ theme(strip.text.x = element_blank(),
                      strip.background = element_rect(colour="white", fill="white"),
                      axis.line.x = element_line(colour = "black"),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.text.x = element_text(size = 15),
                      axis.text.y = element_text(size = 20),
                      axis.title.x.bottom = element_text(size = 20),
                      axis.title.y.left = element_text(size = 20),
                      panel.border = element_rect(colour = "black", fill=NA, size=2),
                      panel.spacing.x = unit(0.5, "lines"),
                      legend.title = element_text(size = 20),
                      legend.text = element_text(face = "italic", size = 15), panel.spacing = unit(0, "lines"))+
    ylab(label = expression(paste( italic("Euphrasia"), " reproductive nodes")))+
    xlab(label = expression(paste("Mean ", italic("Euphrasia"), " reproductive nodes per host")))+
    scale_y_continuous(breaks = c(0,0.6931472, 1.609438,2.302585,3.218876,4.60517), # 1,2, 5, 10, 50, 100
                       labels = round(c(exp(0),exp(0.6931472), exp(1.609438), exp(2.302585), exp(3.218876), exp(4.60517)), 0)) +
    scale_x_continuous(breaks = c(0, 0.6931472, 1.609438, 2.70805), 
                       labels = round(c(exp(0), exp(0.6931472), exp(1.609438), exp(2.70805)), 0)) +
    scale_fill_manual(name = expression(paste(italic("Euphrasia"), " species")), 
                      limits = c("Euphrasia anglica", "Euphrasia micrantha", "Euphrasia tetraquetra", "Euphrasia vigursii"),
                      values = c(cbPalette[2], cbPalette[3], cbPalette[4], cbPalette[5])))

(plot_3.1_tes <- cowplot::plot_grid(plot_3.1, plot_3.1.1, nrow = 2, scale = 0.9))

ggsave(filename = "./Figures/Many_species/population_cum_nodes_mean_tes_bt.pdf", plot = plot_3.1_tes, 
       device = "pdf", width = 17, height = 9, units = "in", useDingbats=FALSE)
ggsave(filename = "./Figures/Many_species/population_cum_nodes_mean_tes_bt.jpeg", plot = plot_3.1_tes, 
       device = "jpeg", width = 17, height = 9, units = "in")

#ggsave(filename = "./Figures/Many_species/population_cum_nodes_mean.pdf", plot = plot_3.1.1, 
 #      device = "pdf", width = 10, height = 8, units = "in", useDingbats=FALSE)
#ggsave(filename = "./Figures/Many_species/population_cum_nodes_mean.jpeg", plot = plot_3.1.1, 
 #      device = "jpeg", width = 10, height = 8, units = "in")

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

