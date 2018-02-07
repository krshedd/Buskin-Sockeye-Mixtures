#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Buskin Sockeye Baseline Update
## Kyle Shedd Mon Feb 05 12:05:25 2018
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Introduction ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The goal of this script is to update the Kodiak sockeye baseline for use on
# analysis of Buskin sockeye subsistence, 2014-2017. Grab the Kodiak portion of
# the KMA baseline used for commercial harvest, 2014-2016. Use the same four
# reporting groups as before for this project: 1) Buskin, 2) Lake Louise,
# 3) Saltery, 4) Other Kodiak


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

## Get objects
KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("OLD", "Likelihood Profiles", "OriginalLocusControl.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 103 Loci
## Locus Control
LocusControl <- dget(file = "Objects/LocusControl103.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get Populations
Kodiak57Pops

invisible(sapply(Kodiak57Pops, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCPopsloci103/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
length(objects(pattern = "\\.gcl"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Kodiak76Collections <- unlist(lapply(Kodiak57Pops, function(silly) {unlist(strsplit(x = silly, split = "\\."))} ))

require(xlsx)
KMA762Collections_SamplesSizes <- read.xlsx(file = "Output/KMA762Collections_SampleSizes.xlsx", sheetIndex = 1, header = TRUE)
str(KMA762Collections_SamplesSizes)

rownames(KMA762Collections_SamplesSizes) <- KMA762Collections_SamplesSizes$NA.
Kodiak76Collections_SamplesSizes <- KMA762Collections_SamplesSizes[Kodiak76Collections, -1]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Baseline with 89 Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/Buskin Subsistence Harvest 2010-2017/Baseline/")

BuskinGroups4 <- c("Buskin Lake", "Lake Louise", "Saltery", "Other Kodiak")
BuskinGroupvec4 <- as.numeric(readClipboard())
# c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 2, 2, 4, 4, 4, 4, 3, 3, 3)
sapply(BuskinGroups4, function(RG) {Kodiak57Pops[BuskinGroupvec4 == which(BuskinGroups4 == RG)]} )  # double check

dput(x = Kodiak57Pops, file = "Objects/Kodiak57Pops.txt")
dput(x = Kodiak76Collections, file = "Objects/Kodiak76Collections.txt")
dput(x = Kodiak76Collections_SamplesSizes, file = "Objects/Kodiak76Collections_SamplesSizes.txt")
dput(x = loci89, file = "Objects/loci89.txt")
dput(x = BuskinGroups4, file = "Objects/BuskinGroups4.txt")
dput(x = BuskinGroupvec4, file = "Objects/BuskinGroupvec4.txt")



Kodiak57Pop89LociBaseline <- CreateBaseline.GCL(sillyvec = Kodiak57Pops, loci = loci89, dir = "BAYES/Baseline",
                                                basename = "Kodiak57Pops89Markers", type = "BAYES")
dput(x = Kodiak57Pop89LociBaseline, file = "Objects/Kodiak57Pop89LociBaseline.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Save in Mixture Directory ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/Buskin Subsistence Harvest 2010-2017/Mixtures/")
dir.create("BAYES/Baseline")
file.copy(from = "V:/Analysis/4_Westward/Sockeye/Buskin Subsistence Harvest 2010-2017/Baseline/BAYES/Baseline/Kodiak57Pops89Markers.bse", to = "BAYES/Baseline/")

dput(x = Kodiak57Pops, file = "Objects/Kodiak57Pops.txt")
dput(x = Kodiak76Collections, file = "Objects/Kodiak76Collections.txt")
dput(x = Kodiak76Collections_SamplesSizes, file = "Objects/Kodiak76Collections_SamplesSizes.txt")
dput(x = loci89, file = "Objects/loci89.txt")
dput(x = BuskinGroups4, file = "Objects/BuskinGroups4.txt")
dput(x = BuskinGroupvec4, file = "Objects/BuskinGroupvec4.txt")
dput(x = Kodiak57Pop89LociBaseline, file = "Objects/Kodiak57Pop89LociBaseline.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Leave One Out ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output <- Anderson_etal.GCL(popvec = Kodiak57Pops, loci = loci89, groups = BuskinGroupvec4, group_names = BuskinGroups4)
str(output)

# Jim's new function does not provide Pop to Pop
# Figure this out with "tidy" tools
genefreq.mod <- output$genefreq
genefreq.mod$ind <- rownames(genefreq.mod)

require(tidyverse)

# Convert to tall format
genofreq.mod <- genefreq.mod %>%
  gather(to_pop, likelihood, -from_pop, -from_group, -ind)
str(genofreq.mod)

# Normalize likelihoods within individuals, then average within from and to pops, then make wide
genofreq.df <- genofreq.mod %>% 
  select(-from_group) %>% 
  group_by(ind) %>% 
  mutate(n_likelihood = likelihood / sum(likelihood)) %>% 
  group_by(from_pop, to_pop) %>% 
  summarise(u_likelihood = mean(n_likelihood)) %>% 
  spread(from_pop, u_likelihood)
str(genofreq.df)

# Convert wide df to matrix for lattice
genofreq.mat <- data.matrix(genofreq.df[, -1])
rownames(genofreq.mat) <- genofreq.df$to_pop
str(genofreq.mat)

require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(genofreq.mat[Kodiak57Pops, Kodiak57Pops]), col.regions = new.colors, xlab = "From_Pop", 
          ylab = "To_Pop", main = "Mean Likelihood", at = seq(0, 1, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })
