# Author: Sunni Patton
# Last edited: 10/21/26
# Title: Creating metadata file
# Overview: Built a metadata file based on original file name structure to use in phyloseq object

## Set seed ====
set.seed(123)

## Load libraries ====
library(ggplot2)
library(dplyr)
library(rstatix)
library(broom)

## Set sample names ====
path <- here::here("Data - raw") # Primers have been removed from these reads

### Set path for forward and reverse reads
reads <- sort(list.files(path, pattern = ".fastq", full.names = TRUE))

### Extract sample names from files
sampleNames <- sapply(strsplit(basename(reads), "_bc"), `[`,1) 

## Extract sample name data ====
time <- as.character(sapply(strsplit(sampleNames, "_T"), `[`,2))
time <- as.character(sapply(strsplit(time, "_"), `[`, 1)) 

treatment <- as.character(sapply(strsplit(sampleNames, "G\\d+\\_T\\d+\\_"), `[`, 2)) 
treatment <- as.character(sapply(strsplit(treatment, "_"), `[`, 1)) 

## Create metadata dataframe ====
sample_df <- data.frame(Samples = sampleNames, Time = time, Treatment = treatment)

## Make Treatment and pretreatment 
sample_df$Treatment[sample_df$Treatment == "B2" | 
                      sample_df$Treatment == "B3" | 
                      sample_df$Treatment == "B4" | 
                      sample_df$Treatment == "B5"] <- "No Treatment"

sample_df$Treatment[sample_df$Treatment == "M1"] <-"Mixture Low"
sample_df$Treatment[sample_df$Treatment == "M2"] <-"Mixture High"
sample_df$Treatment[sample_df$Treatment == "A1"] <-"Ampicillin Low"
sample_df$Treatment[sample_df$Treatment == "A2"] <-"Ampicillin High"
sample_df$Treatment[sample_df$Treatment == "S1"] <-"Streptomycin Low"
sample_df$Treatment[sample_df$Treatment == "S2"] <-"Streptomycin High"
sample_df$Treatment[sample_df$Treatment == "C1"] <-"Ciprofloxacin Low"
sample_df$Treatment[sample_df$Treatment == "C2"] <-"Ciprofloxacin High"

sample_df$Pretreatment <- "No Treatment"

sample_df$Treatment[sample_df$Samples == "NC1"] <-"Negative control"
sample_df$Treatment[sample_df$Samples == "NC2"] <-"Negative control"
sample_df$Treatment[sample_df$Samples == "NC3"] <-"Negative control"

sample_df$Time[sample_df$Samples == "NC1"] <-"Negative control"
sample_df$Time[sample_df$Samples == "NC2"] <-"Negative control"
sample_df$Time[sample_df$Samples == "NC3"] <-"Negative control"

sample_df$Pretreatment[sample_df$Samples == "NC1"] <-"Negative control"
sample_df$Pretreatment[sample_df$Samples == "NC2"] <-"Negative control"
sample_df$Pretreatment[sample_df$Samples == "NC3"] <-"Negative control"

## Save metadata file 
write.csv(sample_df, here::here("Output/metadata_full.rds"))

