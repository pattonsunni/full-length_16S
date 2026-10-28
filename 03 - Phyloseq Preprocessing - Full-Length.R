# Author: Sunni Patton
# Last edited: 10/21/26
# Title: Creating initial phyloseq object

## Set seed ====
set.seed(123)

## Load libraries ====
library(phyloseq)
library(decontam)

## Load metadata file ====
read.csv(here::here("Output/metadata_full.rds")) -> metadata

## Load sequence table and taxonomy table ====
readRDS(here::here("seq_table_nomito_full.rds")) -> seq_table
readRDS(here::here("taxonomy_nomito_full.rds")) -> taxonomy

## Make sample numbers into sample names to work with phyloseq ====
x <- metadata$Samples
rownames(metadata) <- x

## Create phyloseq object ====
ps.All <- phyloseq(otu_table(seq_table, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxonomy))
saveRDS(ps.All, here::here("Output/ps.All_full.rds"))

## Inspect phyloseq object ====
ps.All 
# Taxa distribution 
summary(taxa_sums(ps.All@otu_table)) 
# Sample read distribution
summary(sample_sums(ps.All@otu_table)) 

## Decontam ====
ps.All@sam_data$control <- "FALSE"
ps.All@sam_data$control[ps.All@sam_data$Pretreatment == "Negative control"] <- "TRUE"

contam <- isContaminant(ps.All,
                        method = "prevalence",
                        neg = as.numeric(ps.All@sam_data$control), # don't usually do as.numeric() but was giving an error when I didn't have it 
                        threshold = 0.1) 
table(contam$contaminant) # nothing as contaminant -- maybe because NC reads were so low

## Remove NCs ====
subset_samples(ps.All, Pretreatment != "Negative control") -> ps.noneg
ps.noneg <- prune_taxa(taxa_sums(ps.noneg@otu_table) > 0, ps.noneg) # 333 taxa; 2 taxa only found in NCs
sum(ps.noneg@otu_table) # 1128880 reads; 1237 total reads associated with the 2 ASVs only found in NCs
saveRDS(ps.noneg, here::here("ps.noneg_full.rds"))
