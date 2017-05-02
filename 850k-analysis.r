library(tidyverse)
library(magrittr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

ann.Epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
ann.Epic %>%
  head()

projectDir <- "."

dataDir <- paste(projectDir, "data", sep = "/")

targets <- read.metharray.sheet(dataDir, pattern = "2017-01-20_samplesheet.csv")

rgSet <- read.metharray.exp(targets = targets, verbose = TRUE)

sampleNames(rgSet) <- targets$ID

attr(rgSet, 'annotation')[1] <- 'IlluminaHumanMethylationEPIC'
attr(rgSet, 'annotation')[2] <- 'ilm10b2.hg19'

detP <- detectionP(rgSet)
  detP <- detP[which(detP < 0.05)]

mSetSq <- preprocessQuantile(rgSet)

mSetRaw <- preprocessRaw(rgSet)

plotMDS(getM(mSetSq), top=1000, gene.selection = "common")

# data located at getM(mSetSq); raw data at getM(mSetRaw); ?getM for data extraction




