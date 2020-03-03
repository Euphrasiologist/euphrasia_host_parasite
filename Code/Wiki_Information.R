# A few details to add to the Euphrasia Wiki page #
# Max Brown 12.9.19 #

##### Many host questions #####

##### Part 1: Per-day flowering over the season #####

# load data
floweringsofar<- fread("./Data/Many_hosts/FloweringsofarV2.csv")
# Add underscores to spaces
colnames(floweringsofar) <- gsub(pattern = " ", replacement = "_", x = colnames(floweringsofar))

# some data preparation first
# make columns factors
names_factors <- colnames(floweringsofar)
for (col in names_factors) set(floweringsofar, j=col, value=as.factor(floweringsofar[[col]]))
# days to flower and corolla length change to numeric
names_numeric <- c("Days_since_germination", "Corolla_Length")
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

# When is peak Euphrasia flowering time?
floweringsofar$Flower_date <- as.Date(floweringsofar$Flower_date, format = "%d/%m/%Y")
dates <- floweringsofar[, .(Number_flowered = .N), by = "Flower_date"][order(Flower_date)]
dates$Julian <- as.POSIXlt(dates$Flower_date, format = "%d%b%y")$yday

flowered2017 <- ggplot(dates, aes(x = Julian, y = Number_flowered))+
  geom_point(pch = 1)+
  geom_smooth(method = "loess", se = FALSE, col = "black")+
  theme_bw(base_line_size = 0)+
  xlab(label = "Julian Date of Flowering")+
  ylab(label = "Number of individuals that flowered")

ggsave(filename = "./Figures/Wiki/2017germinants.jpeg", plot = flowered2017, 
       device = "jpeg", width = 6, height = 5, units = "in")

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# mode 172
posterior.mode(dates$Julian[!is.na(dates$Julian)])

##### Many species questions #####

##### Part 1: Per-dat flowering over the season #####

# use expt2dat.2 from Many_species_analysis.R

expt2dat.2$`Flower date` <- gsub(x = expt2dat.2$`Flower date`, pattern = "2019", replacement = "2018", fixed = TRUE)
expt2dat.2$`Flower date` <- as.Date(expt2dat.2$`Flower date`, format = "%d/%m/%Y")

dates2 <- expt2dat.2[, .(Flower_date = `Flower date`,
                           Number_flowered = .N), by = "Flower date"][order(Flower_date)]
dates2$Julian <- as.POSIXlt(dates2$Flower_date, format = "%d%b%y")$yday

flowered2018 <- ggplot(dates2, aes(x = Julian, y = Number_flowered))+
  geom_point(pch = 1)+
  geom_smooth(method = "loess", se = FALSE, col = "black")+
  theme_bw(base_line_size = 0)

posterior.mode(dates2$Julian[!is.na(dates2$Julian)])

ggsave(filename = "./Figures/Wiki/2018germinants.jpeg", plot = flowered2018, 
       device = "jpeg", width = 6, height = 5, units = "in")


dates3 <- expt2dat.2[, .(Flower_date = `Flower date`,
                         Number_flowered = .N), by = c("Flower date", "Euphrasia_sp2")][order(Flower_date)]
dates3$Julian <- as.POSIXlt(dates3$Flower_date, format = "%d%b%y")$yday

flowered2018split <- ggplot(dates3[!is.na(dates3$Euphrasia_sp2)], aes(x = Julian, y = Number_flowered))+
  geom_point(pch = 1)+
  stat_smooth(method = "loess", se = FALSE, aes(colour = Euphrasia_sp2), )+
  theme_bw(base_line_size = 0)

ggsave(filename = "./Figures/Wiki/2018germinantssplit.jpeg", plot = flowered2018split, 
       device = "jpeg", width = 6, height = 5, units = "in")

tapply(X = dates3$Julian, INDEX = dates3$Euphrasia_sp2, FUN = function(x) posterior.mode(x))
