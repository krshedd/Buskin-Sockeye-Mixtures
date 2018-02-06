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
