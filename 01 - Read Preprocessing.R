# 01 - Read Preprocessing
## Samples have been converted from .bam to fastq.gz (no fastqc pr multiqc)


## Resources used 
# https://benjjneb.github.io/LRASManuscript/LRASms_HMP.html
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
## Generally lost a lot of reads. Still have plenty because of how deep sequencing was, but I doubt it's good that
## we lost so many. I read somewhere that allowing for indels can allow for better read retention. For the purposes
## that I need, I don't think I'm too worried about losing reads. Maybe something to come back to later

## Look at length distribution ====
### Make path to files with primers removed
path <- here::here("Data - raw/primers_removed")
hist(nchar(getSequences(path)), 100) # Mostly around 1500 which we would expect

## Filter sequences ====
path <- here::here("Data - raw/primers_removed/filtered")
names <- sort(list.files(path, pattern = ".gz", full.names = FALSE))
filt <- paste0(sapply(strsplit(basename(names), ".fastq.gz"), `[`,1), "_filtered", ".fastq.gz")

track <- filterAndTrim(path, filt, minQ=3, minLen=1000, maxN=0, rm.phix= FALSE, maxEE = 2, verbose = TRUE)

## Dereplicate reads ====
derep <- derepFastq(filt)

## Learn errors ====
err <- learnErrors(derep, BAND_SIZE=32, multithread=TRUE, errorEstimationFunction=dada2:::PacBioErrfun, nbases = 5e8) # band_size is for banded needleman-wunsh; default is 16. Increase if higher indels are expected
# The max qual score of 93 was not detected. Using standard error fitting.
plotErrors(err) # doesn't look good 

## Denoise ====
denoise <- dada(derep, err=err, BAND_SIZE=32, multithread=TRUE)

## Make sequence table ====
### from here to remove chimeras, wasn't in the tutorial; they did assign taxonomy using the denoise object
seq_table <- makeSequenceTable(denoise) 
table(nchar(getSequences(seq_table))) 
sum(seq_table) #1,403,308

### Keep contigs within desired size range
seq_table<-seq_table[,nchar(colnames(seq_table)) %in% 1455:1481]
table(nchar(getSequences(seq_table)))
dim(seq_table) # only 277 contigs
sum(seq_table) # 890,996 reads (lost around 500,000)

seq_table_nochim <- removeBimeraDenovo(seq_table, method="consensus", multithread=TRUE, verbose=TRUE) #Identified 21 bimeras out of 277 input sequences.
dim(seq_table_nochim) # 256 contigs
sum(seq_table) - sum(seq_table_nochim) # lost 14,618 reads

## Assign taxonomy ====
### greengenes2
tax_gg <- assignTaxonomy(seq_table_nochim, "gg2_2024_09_toGenus_trainset.fa", multithread=TRUE) # campylobacterales still not classified beyond order
### silva
tax_silv <- assignTaxonomy(seq_table_nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE) # campylobacterales still not classified beyond order
### GTDB
tax_gtdb <- assignTaxonomy(seq_table_nochim, "GTDB_bac120_arc122_ssu_r202_Genus.fa", multithread=TRUE) # campylobacterales still not classified beyond order
### RDP
tax <- assignTaxonomy(seq_table_nochim, "rdp_19_toGenus_trainset.fa", multithread=TRUE) # says Nitratifractor in nautilales order 
### RefSeq
tax <- assignTaxonomy(seq_table_nochim, "RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa", multithread=TRUE) # campylobacterales still not classified beyond order




saveRDS(seq_table_nochim, here::here("Output/01 - Read Preprocessing - Output/seq_table_nochim.rds"))
saveRDS(tax_gg, here::here("Output/01 - Read Preprocessing - Output/tax_gg.rds"))
saveRDS(tax_silv, here::here("Output/01 - Read Preprocessing - Output/tax_silv.rds"))






samples_short <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))
samples_short <- sapply(strsplit(basename(samples_short), "_bc"), `[`,1) 
