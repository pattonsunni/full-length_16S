## to do:
### look at fastqc/multiqc to determine best truncLen (feel like I shouldn't be losing >50% of reads)

# Author: Sunni Patton
# Last edited: 06/03/2026
# Title: Read preprocessing V4

## Load libraries ====
library(here)
library(dada2)
library(dplyr)

## Set seed for reproducibility ====
set.seed(123)

## Pre-Processing ====

### Set path to raw files

path <- here::here("Data - raw - V4Coral") # Primers have been removed from these reads

### Set path for forward and reverse reads
fnFs <- sort(list.files(path, pattern = "_F.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R.fastq.gz", full.names = TRUE))

### Extract sample names from files
sampleNames <- sapply(strsplit(basename(fnFs), "_F"), `[`,1) 

### Set file destination for files after quality filtering
filtFs <- file.path(path, "filterAndTrim_2026", paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filterAndTrim_2026", paste0(sampleNames, "_R_filt.fastq.gz"))

## Filter and Trim ====
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(245, 230),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # truncLen from G7 paper (could run fastqc and multiqc on just this subset of samples and see if anything changes)
### Save output
saveRDS(out, here::here("out_V4_2026.rds"))

## Assess number of reads lost
sum(out[,1])-sum(out[,2]) # 934,385 reads initially, 352,875 reads left (581,510 reads lost)

## Learn errors and infer sample sequence ====
errF <- learnErrors(filtFs,multithread = TRUE) # default nbases (used all reads)
#86,454,375 total bases in 352,875 reads from 15 samples will be used for learning the error rates.
saveRDS(errF, here::here("errF_V4_2026.rds"))

errR <-learnErrors(filtRs,multithread = TRUE) 
#81,161,250  total bases 352,875 reads from 15 samples samples will be used for learning the error rates.
saveRDS(errR, here::here("errR_V4_2026.rds"))

### Ensure sample naming is consistent
names(filtFs)<-sampleNames
names(filtRs)<-sampleNames

## Infer sample sequence ====
dadaForward <- dada(filtFs, err=errF, multithread=TRUE)
saveRDS(dadaForward, here::here("dadaForward_V4_2026.rds"))

dadaReverse <- dada(filtRs, err=errR, multithread=TRUE)
saveRDS(dadaReverse, here::here("dadaReverse_V4_2026.rds"))

## Create contigs and sequence table ====
contigs <- mergePairs(dadaForward, filtFs, dadaReverse, filtRs)
saveRDS(contigs, here::here("contigs_V4_2026.rds"))

## Make sequence table and visualize contig length and frequency
seq_table <- makeSequenceTable(contigs) 
dim(seq_table) # 108 total contigs (but this includes NCs)

table(nchar(getSequences(seq_table))) # 3 252, 101 253, 4 254

## Remove chimeras ====
### no chimeras detected

## Assign taxonomy ====
taxa <- assignTaxonomy(seq_table, here::here("Training Sets/silva_nr99_v138.1_train_set.fa.gz"), multithread=TRUE)
saveRDS(taxa, here::here("taxa_V4_2026.rds")) 

## Remove off-target sequences ====
### New sequence table (no chloroplast)
is.chloroplast <- taxa[,"Order"] %in% "Chloroplast"
seq_table_nochloro <- seq_table[,!is.chloroplast]
dim(seq_table_nochloro) # 104 contigs; 4 identified as chloroplast
sum(seq_table) - sum(seq_table_nochloro) # 4 contigs associated with 968 reads
saveRDS(seq_table_nochloro, here::here("seq_table_nochloro_V4_2026.rds"))

### New taxonomy table (no chloroplast)
taxonomy_nochloro <- taxa[!is.chloroplast,]
dim(taxonomy_nochloro)
saveRDS(taxonomy_nochloro, here::here("taxonomy_nochloro_V4_2026.rds"))

### no reads associated with mitochondria

## Remove sequences not annotated beyond the Kingdom level
is.NA <- taxonomy_nochloro[,"Phylum"] %in% NA
seq_table_noNA <- seq_table_nochloro[,!is.NA]
dim(seq_table_noNA) # only 1 ASV 

sum(seq_table_nochloro) - sum(seq_table_noNA) # only associated with 14 reads
saveRDS(seq_table_noNA, here::here("seq_table_noNA_V4_2026.rds"))












## Resources used 
# https://benjjneb.github.io/LRASManuscript/LRASms_HMP.html

## Set seed ====
set.seed(123)

## Load libraries ====
library(dada2)

## Create path to files ====
path <- here::here("Data - raw")
samples_long <- sort(list.files(path, pattern = "fastq", full.names = TRUE))
