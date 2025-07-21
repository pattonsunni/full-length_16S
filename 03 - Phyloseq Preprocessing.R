# Author: Sunni Patton
# Last edited: 7/21/25
# Title: Creating initial phyloseq object

## Set seed ====
set.seed(123)

## Load libraries ====
library(phyloseq)

## Load metadata file ====
readRDS(here::here("Output/02 - Metadata - Output/sample_df.rds")) -> sample_df

## Make sample numbers into sample names to work with phyloseq ====
x <- sample_df$Samples
rownames(sample_df) <- x

## Load sequence table and taxonomy table ====
readRDS(here::here("Output/01 - Read Preprocessing - Output/seq_table_nochim.rds")) -> seq_table_nochim
readRDS(here::here("Output/01 - Read Preprocessing - Output/tax_gg.rds")) -> tax_gg

rownames(seq_table_nochim) <- x


## Create phyloseq object ====
ps.All <- phyloseq(otu_table(seq_table_nochim, taxa_are_rows=FALSE), sample_data(sample_df), tax_table(tax_gg))
saveRDS(ps.All, here::here("Output/03 - Phyloseq Preprocessing - Output/ps.All.rds"))

## Inspect phyloseq object ====
ps.All # only 256 taxa somehow
# Taxa distribution 
summary(taxa_sums(ps.All@otu_table)) 
# Sample read distribution
summary(sample_sums(ps.All@otu_table)) 
