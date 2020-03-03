## NOT TO BE UPLOADED TO GITHUB ##

# these are just some other aspects of the species data which included host-parasite co-occurrence

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
