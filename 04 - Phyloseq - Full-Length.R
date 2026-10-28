# Author: Sunni Patton
# Last edited: 10/23/26
# Title: Sample processing in phyloseq
# Overview: Adding phylogenetic tree to phyloseq object and randomly subsampling to an even depth

## Set seed ====
set.seed(123)

## Load libraries ====
library(speedyseq)
library(tidyr)
library(phytools)
library(seqateurs)
library(vegan)
library(microViz)
library(ips)
library(ggplot2)

## Load phyloseq object (contaminants and NCs removed) ====
readRDS(here::here("ps.noneg_full.rds")) -> ps.noneg

## Add ASV to taxonomy table ====
ps.noneg <- ps.noneg %>% mutate_tax_table(ASV = paste0("ASV", 1:333))

## Change DNA sequence to ASV number ====
dna <- Biostrings::DNAStringSet(taxa_names(ps.noneg))
names(dna) <- taxa_names(ps.noneg)
ps.noneg <- merge_phyloseq(ps.noneg, dna)
taxa_names(ps.noneg) <- paste0("ASV", seq(ntaxa(ps.noneg)))

saveRDS(ps.noneg, here::here("ps.noneg_full.rds"))

## See how many ASVs annotated to the genus level
subset_taxa(ps.noneg, Genus != "NA") -> ps.noneg.genus #172 out of 333 ASVs annotated to genus level

## Make rarefaction curve ====
as.matrix(as.data.frame(ps.noneg@otu_table)) -> data
rarecurve(data, step = 100, label = FALSE, abline(v = "9961", col = "red", lty = 2, lwd = 1.5)) #does not fully capture diversity
rarecurve(data, step = 100, label = FALSE, abline(v = "64606", col = "red", lty = 2, lwd = 1.5)) #much better but would need to remove the 2 samples that fall below

# already checked RA to see if the general profile of these two samples are the same, and they are. Just don't want to have these skew anything down the line

## Remove low read samples ====
ps.noneg_sub <- subset_samples(ps.noneg, Samples != "G7_T48_B3_a" & Samples != "G7_T48_C1_3_a")
ps.noneg_sub <- prune_taxa(taxa_sums(ps.noneg_sub@otu_table) > 0, ps.noneg_sub) #322 taxa (11 removed from removing these samples)

subset_taxa(ps.noneg_sub, Genus != "NA") -> ps.noneg_sub.genus # 165 out of 322 annotated to genus level after removing low read samples

## Remove low abundance taxa ====
### Look at abundance and frequency of each taxon
### Doing this makes it easy to see in how many samples each taxon appears (Prevalence), and at what abundance (TotalCounts)
pst <- fast_melt(ps.noneg_sub) # Takes a phyloseq object and returns a melted dataframe with sample ID, variable name, and abundance value for each observation

prevdt <- pst[, list(Prevalence = sum(count > 0), 
                     TotalCounts = sum(count)),
              by = taxaID]
View(prevdt) ## lots of ASVs that are low and only present in one sample
summary(taxa_sums(ps.noneg_sub@otu_table)) # 1st quartile is 21

### Prune taxa that are only present in 1 sample and with total reads < 1st quartile (21) ====
summary(taxa_sums(ps.noneg_sub@otu_table))
keepTaxa <- prevdt[((Prevalence >1 | TotalCounts >21)), taxaID]
ps.pruned <- prune_taxa(keepTaxa, ps.noneg_sub) # 243 taxa remain after pruning (79 ASVs removed)
saveRDS(ps.pruned, here::here("ps.pruned_full.rds"))

### look at taxa removed
subset_taxa(ps.pruned, Genus != "NA") -> ps.pruned.genus 
#114 taxa annotated to genus level out of 243

## Rarefy even depth ====
rarefy_even_depth(ps.pruned, 
                  sample.size = min(sample_sums(ps.pruned)), 
                  trimOTUs = TRUE) -> ps.rare # all ASVs still remain after rarefying
sample_sums(ps.rare) # rarefied to 64,574

saveRDS(ps.rare, here::here("ps.rare_full.rds"))

## want to look into which ASVs were removed from pruning
## I still feel like the choice to prune is valid because these are all supposed to be essentially the same samples
## so if an ASV is in only one sample and at very low abundance, then it probably isn't very relevant