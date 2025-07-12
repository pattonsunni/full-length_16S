# 01 - Read Preprocessing
## Samples have been converted from .bam to fastq.gz (no fastqc pr multiqc)


## Set seed
set.seed(123)
## Load libraries 
library(dada2)
## Create path to files
path <- here::here("Data - raw")
samples_long <- sort(list.files(path, pattern = "fastq", full.names = TRUE))
samples_short <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))
samples_short <- sapply(strsplit(basename(samples_short), "_bc"), `[`,1) 

## Remove primers
### Set primer sequences
F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"
## Remove primers from each sample
for (i in samples_long) {
  noprimers <- tempfile(fileext = ".fastq")
  removePrimers(i, noprimers, primer.fwd=F27, primer.rev=rc(R1492), orient=TRUE, verbose=TRUE) 
}

