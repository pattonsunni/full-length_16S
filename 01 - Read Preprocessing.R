# 01 - Read Preprocessing
## Samples have been converted from .bam to fastq.gz (no fastqc pr multiqc)

## Set seed ====
set.seed(123)
## Load libraries ====
library(dada2)
## Create path to files ====
path <- here::here("Data - raw")
samples_long <- sort(list.files(path, pattern = "fastq", full.names = TRUE))

## Remove primers ====
### Set primer sequences
F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"

for (i in samples_long) {
  filepath <- paste0(sapply(strsplit(basename(i), ".fastq"), `[`,1), "_rmprimer", ".fastq.gz")
  removePrimers(i, filepath, primer.fwd=F27, primer.rev=rc(R1492), orient=TRUE, verbose=TRUE) 
  print(paste("Done with sample", basename(i)))
  }

## 


samples_short <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))
samples_short <- sapply(strsplit(basename(samples_short), "_bc"), `[`,1) 
# G7_T0_B2_b