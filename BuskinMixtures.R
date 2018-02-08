#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Buskin Sockeye Mixtures 2014-2017
## Kyle Shedd Mon Feb 05 13:31:20 2018
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Introduction ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The goal of this script is to analyze Buskin sockeye subsistence, 2014-2017, 
# using a Kodiak Archipelago baseline from the Kodiak portion of
# the KMA baseline used for commercial harvest, 2014-2016. Use the same four
# reporting groups as before for this project: 1) Buskin, 2) Lake Louise,
# 3) Saltery, 4) Other Kodiak


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/Buskin Subsistence Harvest 2010-2017/Mixtures/")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pull genotypes from LOKI 2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
username = "krshedd"

CreateLocusControl.GCL(markersuite = "Sockeye2011_96SNPs", username = username, password = password)

## Save original LocusControl
loci96 <- LocusControl$locusnames
mito.loci <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl.txt")
dput(x = loci96, file = "Objects/loci96.txt")
dput(x = mito.loci, file = "Objects/mito.loci.txt")

#~~~~~~~~~~~~~~~~~~
## Pull all data for each silly code and create .gcl objects for each
Buskin2014_2017 <- paste0("SBUSKS", 14:17)
dput(x = Buskin2014_2017, file = "Objects/Buskin2014_2017.txt")

LOKI2R.GCL(sillyvec = Buskin2014_2017, username = username, password = password)
rm(username, password)

## Save unaltered .gcl's as back-up:
# dir.create("Raw genotypes/")
# dir.create("Raw genotypes/OriginalCollections")
invisible(sapply(Buskin2014_2017, function(silly) {dput(x = get(paste0(silly, ".gcl")), file = paste0("Raw genotypes/OriginalCollections/" , silly, ".txt"))} )); beep(8)

## Original sample sizes by SILLY
collection.size.original <- sapply(Buskin2014_2017, function(silly) get(paste0(silly, ".gcl"))$n)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/Buskin Subsistence Harvest 2010-2017/Mixtures/")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

## Get objects
LocusControl <- dget(file = "Objects/OriginalLocusControl.txt")
KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects
invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)

## Get un-altered mixtures
invisible(sapply(Buskin2014_2017, function(silly) {assign(x = paste0(silly, ".gcl"), value = dget(file = paste0(getwd(), "/Raw genotypes/OriginalCollections/", silly, ".txt")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run all fish from each year as it's own mixture per Op Plan
# Samples were collected as representatively as possible
# Subsampling was random within year


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Buskin2014_2017

Buskin2014_2017_SampleSizes <- matrix(data = NA, nrow = length(Buskin2014_2017), ncol = 5, 
                                      dimnames = list(Buskin2014_2017, c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_Buskin2014_2017_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = Buskin2014_2017, loci = loci96)
min(Original_Buskin2014_2017_SampleSizebyLocus)  ## 267/285
apply(Original_Buskin2014_2017_SampleSizebyLocus, 1, min) / apply(Original_Buskin2014_2017_SampleSizebyLocus, 1, max) # 2016 was 0.83, all others > 0.95


#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_Buskin2014_2017_ColSize <- sapply(paste(Buskin2014_2017, ".gcl", sep = ''), function(x) get(x)$n)
Buskin2014_2017_SampleSizes[, "Genotyped"] <- Original_Buskin2014_2017_ColSize


### Alternate
## Indentify alternate species individuals
Buskin2014_2017_Alternate <- FindAlternateSpecies.GCL(sillyvec = Buskin2014_2017, species = "sockeye")

## Remove Alternate species individuals
RemoveAlternateSpecies.GCL(AlternateSpeciesReport = Buskin2014_2017_Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5)

## Get number of individuals per silly after removing alternate species individuals
ColSize_Buskin2014_2017_PostAlternate <- sapply(paste(Buskin2014_2017, ".gcl", sep = ''), function(x) get(x)$n)
Buskin2014_2017_SampleSizes[, "Alternate"] <- Original_Buskin2014_2017_ColSize-ColSize_Buskin2014_2017_PostAlternate


### Missing
## Remove individuals with >20% missing data
Buskin2014_2017_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = Buskin2014_2017, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_Buskin2014_2017_PostMissLoci <- sapply(paste(Buskin2014_2017, ".gcl", sep = ''), function(x) get(x)$n)
Buskin2014_2017_SampleSizes[, "Missing"] <- ColSize_Buskin2014_2017_PostAlternate-ColSize_Buskin2014_2017_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
Buskin2014_2017_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = Buskin2014_2017, loci = loci96, quantile = NULL, minproportion = 0.95)
Buskin2014_2017_DuplicateCheckReportSummary <- sapply(Buskin2014_2017, function(x) Buskin2014_2017_DuplicateCheck95MinProportion[[x]]$report)

## Remove duplicate individuals
Buskin2014_2017_RemovedDups <- RemoveDups.GCL(Buskin2014_2017_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_Buskin2014_2017_PostDuplicate <- sapply(paste(Buskin2014_2017, ".gcl", sep = ''), function(x) get(x)$n)
Buskin2014_2017_SampleSizes[, "Duplicate"] <- ColSize_Buskin2014_2017_PostMissLoci-ColSize_Buskin2014_2017_PostDuplicate


### Final
Buskin2014_2017_SampleSizes[, "Final"] <- ColSize_Buskin2014_2017_PostDuplicate
Buskin2014_2017_SampleSizes

write.csv(Buskin2014_2017_SampleSizes, file = "Output/Buskin2014_2017_SampleSizes.csv")
dput(x = Buskin2014_2017_SampleSizes, file = "Objects/Buskin2014_2017_SampleSizes.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Combine Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Combine loci
combined.loci.89 <- sapply(grep(pattern = "\\.", x = loci89, value = TRUE), function(locus) {unlist(strsplit(x = locus, split = "\\."))}, simplify = FALSE)
sapply(combined.loci.89, function(loci2combine) {CombineLoci.GCL(sillyvec = Buskin2014_2017, markerset = loci2combine, update = TRUE)} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save PostQC/Combined loci .gcl's as back-up:
dir.create("Raw genotypes/OriginalCollections_PostQC_CombinedLoci")
invisible(sapply(Buskin2014_2017, function(silly) {dput(x = get(paste0(silly, ".gcl")), file = paste0("Raw genotypes/OriginalCollections_PostQC_CombinedLoci/" , silly, ".txt"))} )); beep(8)

## Dput LocusControl
dput(x = LocusControl, file = "Objects/LocusControl98.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get/Create New 4RG MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Flat prior for all
Kodiak57Pops4FlatPrior <- Prior.GCL(groupvec = BuskinGroupvec4, groupweights = rep(1 / 4, 4), minval = 0.01)
Kodiak57PopsInits <- MultiChainInits.GCL(npops = 57, nchains = 5, prop = 0.9)
WASSIPSockeyeSeeds <- dget("V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt")

dput(x = Kodiak57Pops4FlatPrior, file = "Objects/Kodiak57Pops4FlatPrior.txt")
dput(x = Kodiak57PopsInits, file = "Objects/Kodiak57PopsInits.txt")
dput(x = WASSIPSockeyeSeeds, file = "Objects/WASSIPSockeyeSeeds.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create BAYES Files ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dumping Mixture files
Kodiak57Pop89LociMixtureFormat <- CreateMixture.GCL(sillys = Buskin2014_2017[1], loci = loci89, IDs = NULL, mixname = Buskin2014_2017[1],
                                           dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(Kodiak57Pop89LociMixtureFormat, file = "Objects/Kodiak57Pop89LociMixtureFormat.txt")

sapply(Buskin2014_2017, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci89, IDs = NULL, mixname = Mix, dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Buskin2014_2017, function(Mix) {
  CreateControlFile.GCL(sillyvec = Kodiak57Pops, loci = loci89, mixname = Mix, basename = "Kodiak57Pops89Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = BuskinGroupvec4, priorvec = Kodiak57Pops4FlatPrior, initmat = Kodiak57PopsInits, dir = "BAYES/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = Kodiak57Pop89LociMixtureFormat, basefortran = Kodiak57Pop89LociBaseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Buskin2014_2017, function(Mix) {dir.create(paste0("BAYES/Output/", Mix))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/Buskin Subsistence Harvest 2010-2017/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

## Get objects
LocusControl <- dget(file = "Objects/LocusControl98.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize BAYES Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Buskin2014_2017_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(BuskinGroups4), groupnames = BuskinGroups4 ,
  maindir = "BAYES/Output", 
  mixvec = Buskin2014_2017, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)

# dput
dir.create("Estimates objects")
dput(x = Buskin2014_2017_Estimates, file = "Estimates objects/Buskin2014_2017_Estimates.txt")
Buskin2014_2017_Estimates <- dget(file = "Estimates objects/Buskin2014_2017_Estimates.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Baseline
require(xlsx)
require(tidyverse)
kma_baseline <- read.xlsx(file = "V:/Documents/4_Westward/Sockeye/KMA 2014-2016/KMA 2015 Baseline FMS/FMS 16-XX KMA Sockeye 2015 Baseline Tables.xlsx", sheetName = "Table 3. - Collections", startRow = 3, header = TRUE, stringsAsFactors = FALSE)
str(kma_baseline)

pop2colection <- sapply(Kodiak57Pops, function(pop) length(unlist(strsplit(x = pop, split = "\\."))))

buskin_baseline <- kma_baseline %>% 
  filter(ADF.G.code %in% Kodiak76Collections) %>% 
  mutate(Date = format(as.Date(as.numeric(Date), origin = "1899-12-30"), format = "%m/%d/%Y")) %>% 
  rename(Reporting_Group = NA.) %>% 
  mutate(Reporting_Group = BuskinGroups4[rep(BuskinGroupvec4, pop2colection)]) %>% 
  mutate(Collection = Collection - min(Collection) + 1) %>% 
  mutate(Population = Population - min(Population) + 1)
  
str(buskin_baseline)

# Any weird dates?
sort(as.Date(buskin_baseline$Date, format = "%m/%d/%Y"))  # 1905-06-15???, this was for SLUP93 which only has year (1993)

# Replace wayward date for SLUP93
buskin_baseline$Date[buskin_baseline$Date == format(min(as.Date(buskin_baseline$Date, format = "%m/%d/%Y")), "%m/%d/%Y")] <- "1993"

# Write file
write.xlsx(x = buskin_baseline, file = "Tables/Buskin 2014-2017.xlsx", sheetName = "Baseline Collections", col.names = TRUE, row.names = FALSE, showNA = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## QC
# Read in relevant QC data
silly.summary <- read.xlsx(file = "V:/Lab/Genotyping/SNP Projects/Sockeye/Project S177 Buskin Subsistence 2014_2017/QC/Project S177 QC Summary.xlsx", sheetName = "Summary by Silly", stringsAsFactors = FALSE)
str(silly.summary)

conflict.summary <- read.xlsx(file = "V:/Lab/Genotyping/SNP Projects/Sockeye/Project S177 Buskin Subsistence 2014_2017/QC/Project S177 QC Summary.xlsx", sheetName = "Conflicts by Silly", stringsAsFactors = FALSE)
str(conflict.summary)

# Join
QC.df <- left_join(x = silly.summary, y = conflict.summary, by = "NA.") %>% 
  rename(Silly = NA.) %>% 
  mutate(Year = 2014:2017) %>% 
  mutate(HomoHet = Total.Het.Homo + Total.Homo.Het) %>% 
  mutate(HomoHetRate = HomoHet / Total.QC.Genotypes) %>% 
  mutate(ErrorRate = Discrepancy.Rate / 2)

QC_table <- QC.df[, c("Year", "Genotyped", "Total.QC.Fish", "Failure.Rate", "Total.QC.Genotypes", "HomoHet", "HomoHetRate", "Total.Homo.Homo", "Homo.Homo.Rate", "Discrepancy.Rate", "ErrorRate")]
colnames(QC_table) <- c("Year", "Original.n", "QC.n", "Failure.r", "QC.Genotypes", "Homo-Het.n", "Homo-Het.r", "Homo-Homo.n", "Homo-Homo.n", "Discrepancy.r", "Error.r")

DataQC_table <- QC.df[, c("Year", "Genotyped", "Alternate", "Missing", "Duplicate", "Final")]


# Write file
write.xlsx(x = QC_table, file = "Tables/Buskin 2014-2017.xlsx", sheetName = "QC", col.names = TRUE, row.names = FALSE, append = TRUE)
write.xlsx(x = DataQC_table, file = "Tables/Buskin 2014-2017.xlsx", sheetName = "Sample Sizes", col.names = TRUE, row.names = FALSE, append = TRUE)
