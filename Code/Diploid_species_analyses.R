## Do we see differences in morphology between diploid species? ##
# Created: 9.8.19 by Max Brown #

## Libraries needed ##

library(data.table)
library(ggplot2)
library(MCMCglmm)
library(lme4)
library(VCVglmm)
library(dplyr)

##### Part 1: Data import and tidying #####

# from Many_species_analyses.R
diploids <- fread("./Data/Many_species/expt2dat.2.csv")
# get out the diploids!
diploids1 <- diploids[grepl(pattern = "B", x = Unique_ID),]
diploids2 <- diploids[Population2 %in% c("A1766", "V1761"),][Host_code == "HLA" | Host_code == "PLA",]
diploids3 <- rbind(diploids1, diploids2)
# remove single vigursii...
diploids3 <- diploids3[Population2 != "V1761"]

##### Part 2: Principal Component Analysis #####

# what are the names of the populations?
diploids3[, .(N = .N), by = Population2]

# julian days 
diploids3$Julian <- as.POSIXlt(as.Date(diploids3$`Flower date`), format = "%d/%m%Y")$yday
diploids3$Germ_dateJ <- as.POSIXlt(as.Date(diploids3$Germ_date), format = "%d/%m%Y")$yday
diploids3$Trans_dateJ <- as.POSIXlt(as.Date(diploids3$Trans_date), format = "%d/%m%Y")$yday

diploids3$DaysSinceGerm <- diploids3$Julian - diploids3$Germ_dateJ

# do a PCA for 7 taxonomic traits
# set type of columns
names_prcomp <- c("Node", "Height", "Corolla size", "Teeth", "Branches", "Cauline ratio", "DaysSinceGerm")
for (col in names_prcomp) set(diploids3, j=col, value=as.numeric(diploids3[[col]]))

# locate NA's
which(is.na(diploids3[Euphrasia_sp2 != "Euphrasia rivularis", c("Node", "Height", "Corolla size", "Teeth", "Branches", "Cauline ratio", "DaysSinceGerm")]), arr.ind=TRUE)
# remove NA's
diploids4 <- diploids3[Euphrasia_sp2 != "Euphrasia rivularis", c("Host_code","Population2","Euphrasia_sp2","Node", "Height", "Corolla size", "Teeth", "Branches", "Cauline ratio", "DaysSinceGerm")][-c(243,80),]

pca_diploids <- prcomp(diploids4[,-c("Host_code","Population2","Euphrasia_sp2")], center = TRUE, scale. = TRUE)

# add in the factors we may want to group by
pca_diploids2 <- setDT(as.data.frame(pca_diploids$x))
pca_diploids2$Species <- diploids4$Euphrasia_sp2
pca_diploids2$Population <- diploids4$Population2
pca_diploids2$Hosts <- diploids4$Host_code

# make the plot

ggplot(pca_diploids2, aes(x = PC1, y = PC2))+
  geom_point(aes(colour = Species))+
  stat_ellipse(aes(colour = Species))+
  theme_bw()

##### Part 3: LDA  #####

# 279 observations

train <- sample(1:279, 75)

table(diploids4$Euphrasia_sp2[train])

z <- MASS::lda(Euphrasia_sp2 ~ ., diploids4[,-c("Population2")], prior = c(1,1,1)/3, subset = train)

predz <- predict(z, diploids4[-train, -"Population2"])

ggplot(as.data.frame(predz$x), aes(x = LD1, y = LD2))+
  geom_point(aes(colour = predz$class))

MASS::ldahist(data = predz$x, g=diploids4[,-c("Population2")]$Euphrasia_sp2)

(tapply(X = predz$x[,1], INDEX = predz$class, FUN = mean)[1]+tapply(X = predz$x[,1], INDEX = predz$class, FUN = mean)[2])/2
(tapply(X = predz$x[,1], INDEX = predz$class, FUN = mean)[2]+tapply(X = predz$x[,1], INDEX = predz$class, FUN = mean)[3])/2

calcAllocationRuleAccuracy <- function(ldavalue, groupvariable, cutoffpoints)
{
  # find out how many values the group variable can take
  groupvariable2 <- as.factor(groupvariable[[1]])
  levels <- levels(groupvariable2)
  numlevels <- length(levels)
  # calculate the number of true positives and false negatives for each group
  numlevels <- length(levels)
  for (i in 1:numlevels)
  {
    leveli <- levels[i]
    levelidata <- ldavalue[groupvariable==leveli]
    # see how many of the samples from this group are classified in each group
    for (j in 1:numlevels)
    {
      levelj <- levels[j]
      if (j == 1)
      {
        cutoff1 <- cutoffpoints[1]a
        cutoff2 <- "NA"
        results <- summary(levelidata <= cutoff1)
      }
      else if (j == numlevels)
      {
        cutoff1 <- cutoffpoints[(numlevels-1)]
        cutoff2 <- "NA"
        results <- summary(levelidata > cutoff1)
      }
      else
      {
        cutoff1 <- cutoffpoints[(j-1)]
        cutoff2 <- cutoffpoints[(j)]
        results <- summary(levelidata > cutoff1 & levelidata <= cutoff2)
      }
      trues <- results["TRUE"]
      trues <- trues[[1]]
      print(paste("Number of samples of group",leveli,"classified as group",levelj," : ",
                  trues,"(cutoffs:",cutoff1,",",cutoff2,")"))
    }
  }
}

calcAllocationRuleAccuracy(ldavalue = predz$x[,1], predz$class, cutoffpoints = c(0.264, 0.733))

# so actually LDA can't predict so well.

##### Part 3: #####

ggplot(diploids3[Euphrasia_sp2 != "Euphrasia rivularis"], aes(x = Euphrasia_sp2, y = Node))+
  geom_jitter()+
  geom_boxplot(alpha = 0)

ggplot(diploids3, aes(x = DaysSinceGerm, y = `Corolla size`))+
  geom_jitter(aes(colour = Euphrasia_sp2,
                  shape = Host_code))

