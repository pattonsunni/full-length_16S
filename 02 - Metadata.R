# Author: Sunni Patton
# Last edited: 7/21/25
# Title: Creating metadata file
# Overview: Built a metadata file based on original file name structure to use in phyloseq object

## Set seed ====
set.seed(123)

## Load libraries ====
library(ggplot2)
library(dplyr)
library(rstatix)
library(broom)

## Extract sample name data ====
### Extract information from samples_short
genotype <- as.character(sapply(strsplit(samples_short, "_T"), `[`,1))

time <- as.character(sapply(strsplit(samples_short, "_T"), `[`,2))
time <- as.character(sapply(strsplit(time, "_"), `[`, 1)) 

treatment <- as.character(sapply(strsplit(samples_short, "G\\d+\\_T\\d+\\_"), `[`, 2)) 
treatment <- as.character(sapply(strsplit(treatment, "_"), `[`, 1)) 

## Create metadata dataframe ====
sample_df <- data.frame(Samples = samples_short, Time = time, Treatment = treatment)

## Make Treatment long 
sample_df$Treatment[sample_df$Treatment == "B2" | 
                         sample_df$Treatment == "B3" | 
                         sample_df$Treatment == "B4" | 
                         sample_df$Treatment == "B5"] <- "Control"

sample_df$Treatment[sample_df$Treatment == "M1"] <-"Mixture Low"
sample_df$Treatment[sample_df$Treatment == "M2"] <-"Mixture High"
sample_df$Treatment[sample_df$Treatment == "A1"] <-"Ampicillin Low"
sample_df$Treatment[sample_df$Treatment == "A2"] <-"Ampicillin High"
sample_df$Treatment[sample_df$Treatment == "S1"] <-"Streptomycin Low"
sample_df$Treatment[sample_df$Treatment == "S2"] <-"Streptomycin High"
sample_df$Treatment[sample_df$Treatment == "C1"] <-"Ciprofloxacin Low"
sample_df$Treatment[sample_df$Treatment == "C2"] <-"Ciprofloxacin High"

## Add Antibiotic column 
sample_df$Antibiotic <- "Negative Control" # Creating a new column in our dataframe and using negative control as a placeholder

### Add correct values
sample_df$Antibiotic[sample_df$Treatment == "Control"] <- "Control"

sample_df$Antibiotic[sample_df$Treatment == "Mixture Low" |
                          sample_df$Treatment == "Mixture High"] <- "Mixture"

sample_df$Antibiotic[sample_df$Treatment == "Ampicillin Low" |
                          sample_df$Treatment == "Ampicillin High"] <- "Ampicillin"

sample_df$Antibiotic[sample_df$Treatment == "Streptomycin Low" |
                          sample_df$Treatment == "Streptomycin High"] <- "Streptomycin"

sample_df$Antibiotic[sample_df$Treatment == "Ciprofloxacin Low" |
                          sample_df$Treatment == "Ciprofloxacin High"] <- "Ciprofloxacin"

## Add Pretreatment column(s)
### Adding these to indicate that all T0 groups are actually pretreatment 

### Add Pretreatment column
sample_df$Pretreatment <- sample_df$Treatment
### Add information to Pretreatment column
sample_df$Pretreatment[sample_df$Time == "0"] <- "Control"

## Save metadata file 
saveRDS(sample_df, here::here("Output/02 - Metadata - Output/sample_df.rds"))
