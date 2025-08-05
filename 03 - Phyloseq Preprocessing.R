# Author: Sunni Patton
# Last edited: 7/21/25
# Title: Creating initial phyloseq object

## Set seed ====
set.seed(123)

## Load libraries ====
library(phyloseq)
library(speedyseq)
library(microViz)

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
ps.All_gg <- phyloseq(otu_table(seq_table_nochim, taxa_are_rows=FALSE), sample_data(sample_df), tax_table(tax_gg))
saveRDS(ps.All_gg, here::here("Output/03 - Phyloseq Preprocessing - Output/ps.All_gg.rds"))

ps.All_silv <- phyloseq(otu_table(seq_table_nochim, taxa_are_rows=FALSE), sample_data(sample_df), tax_table(tax_silv))
saveRDS(ps.All_silv, here::here("Output/03 - Phyloseq Preprocessing - Output/ps.All_silv.rds"))

ps.All_gtdb <- phyloseq(otu_table(seq_table_nochim, taxa_are_rows=FALSE), sample_data(sample_df), tax_table(tax_gtdb))
saveRDS(ps.All_gtdb, here::here("Output/03 - Phyloseq Preprocessing - Output/ps.All_gtdb.rds"))

ps.All_refseq <- phyloseq(otu_table(seq_table_nochim, taxa_are_rows=FALSE), sample_data(sample_df), tax_table(tax_refseq))
ps.All_rdp <- phyloseq(otu_table(seq_table_nochim, taxa_are_rows=FALSE), sample_data(sample_df), tax_table(tax_rdp))


## Inspect phyloseq object ====
ps.All_silv # only 274 taxa somehow
# Taxa distribution 
summary(taxa_sums(ps.All_silv@otu_table)) 
# Sample read distribution
summary(sample_sums(ps.All_silv@otu_table)) 
# Remove any taxa not annotated to kingdom level 
## gtdb has 2 ASVs not annotated to kingdom level; one is hitting to A. palmata and the other isn't hitting to anything
sum(taxa_sums(ps.All_silv)) # 876378
subset_taxa(ps.All_silv, Kingdom != "NA") -> ps.All_silv
sum(taxa_sums(ps.All_silv)) # 872,641

## Add ASV column to taxonomy table ====
ps.All_silv <- ps.All_silv %>% mutate_tax_table(ASV = paste0("ASV", 1:274))

## Save taxonomy table (only DNA sequence and ASV columns) ====
### Save taxonomy table as tibble
as_tibble(ps.All_silv@tax_table) -> taxa_tibble
### Save only relevent columns
taxa_df <- data.frame(taxa_tibble$.otu, taxa_tibble$ASV)
### Rename columns
colnames(taxa_df) <- c("Sequence", "ASV")

readr::write_csv(taxa_df, here::here("Output/03 - Phyloseq Preprocessing - Output/sequenceASV_gtdb.csv"))

## Fix taxonomy and add new taxonomy column ====
tax_fix(
  ps.All_silv,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
) -> ps.All_silv

mutate_tax_table(ps.All_silv, ASVnew = paste0(Genus, " ", ASV)) -> ps.All_silv


## Relative abundance ====
ps.All_silv.trans <- transform_sample_counts(ps.All_silv, function(OTU) OTU/sum(OTU))

## Identify top 10 most abundant taxa
top20 <- names(sort(taxa_sums(ps.All_silv.trans), decreasing = TRUE))[1:20] 

## Prune top 10
ps.rare.top20 <- prune_taxa(top20, ps.All_silv.trans)

## Relative abundance ====


plot.relAbund.all <- plot_bar(ps.rare.top20, x="Samples", fill="ASVnew") 

plot.relAbund.all <- plot.relAbund.all + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("Samples") + labs(title = "Relative Abundance")



## Subset only G50 samples ====
subset_samples(ps.All_gtdb, Samples == "G50_T0_M1_3_b" | Samples == "G50_T0_M1_3_d") -> ps.G50_gtdb
ps.G50_gtdb <- prune_taxa(taxa_sums(ps.G50_gtdb@otu_table) > 0, ps.G50_gtdb)

ps.G50_gtdb.trans <- transform_sample_counts(ps.G50_gtdb, function(OTU) OTU/sum(OTU))

## Identify top 10 most abundant taxa
top10 <- names(sort(taxa_sums(ps.G50_gtdb.trans), decreasing = TRUE))[1:10] 

## Prune top 10
ps.rare.top10 <- prune_taxa(top10, ps.G50_gtdb.trans)

plot.relAbund.all <- plot_bar(ps.rare.top10, x="Samples", fill="ASVnew") 

plot.relAbund.all <- plot.relAbund.all + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("Samples") + labs(title = "Relative Abundance")

