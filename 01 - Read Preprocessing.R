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
err <- learnErrors(derep,multithread=TRUE, errorEstimationFunction=dada2:::PacBioErrfun, nbases = 5e8) # band_size is for banded needleman-wunsh; default is 16. Increase if higher indels are expected
# The max qual score of 93 was not detected. Using standard error fitting.
plotErrors(err) # doesn't look good 


## Edits to learning errors (from github) ====
loessErrfun_mod4 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # jonalim's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 0.95)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_4 <- learnErrors(
  derep,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod4,
  verbose = TRUE
)


## Denoise ====
denoise <- dada(derep, err=errF_4, BAND_SIZE=32, multithread=TRUE)

## Make sequence table ====
### from here to remove chimeras, wasn't in the tutorial; they did assign taxonomy using the denoise object
seq_table <- makeSequenceTable(denoise) 
table(nchar(getSequences(seq_table))) 
sum(seq_table) #1,398,080

### Keep contigs within desired size range
seq_table<-seq_table[,nchar(colnames(seq_table)) %in% 1455:1481]
table(nchar(getSequences(seq_table)))
dim(seq_table) # only 302 contigs
sum(seq_table) # 891,039 reads (lost around 500,000)

seq_table_nochim <- removeBimeraDenovo(seq_table, method="consensus", multithread=TRUE, verbose=TRUE) #Identified 28 bimeras out of 302 input sequences.
dim(seq_table_nochim) # 274 contigs
sum(seq_table) - sum(seq_table_nochim) # lost 18,398 reads

## Assign taxonomy ====
### greengenes2
tax_gg <- assignTaxonomy(seq_table_nochim, "gg2_2024_09_toGenus_trainset.fa", multithread=TRUE) # campylobacterales still not classified beyond order
### silva
tax_silv <- assignTaxonomy(seq_table_nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE) # campylobacterales still not classified beyond order
### GTDB
tax_gtdb <- assignTaxonomy(seq_table_nochim, "GTDB_bac120_arc122_ssu_r202_Genus.fa", multithread=TRUE) # campylobacterales still not classified beyond order
### RDP
tax_rdp <- assignTaxonomy(seq_table_nochim, "rdp_19_toGenus_trainset.fa", multithread=TRUE) # says Nitratifractor in nautilales order 
### RefSeq
tax_refseq <- assignTaxonomy(seq_table_nochim, "RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa", multithread=TRUE) # campylobacterales still not classified beyond order


saveRDS(seq_table_nochim, here::here("Output/01 - Read Preprocessing - Output/seq_table_nochim.rds"))
saveRDS(tax_gg, here::here("Output/01 - Read Preprocessing - Output/tax_gg.rds"))
saveRDS(tax_silv, here::here("Output/01 - Read Preprocessing - Output/tax_silv.rds"))
saveRDS(tax_gtdb, here::here("Output/01 - Read Preprocessing - Output/tax_gtdb.rds"))
saveRDS(tax_rdp, here::here("Output/01 - Read Preprocessing - Output/tax_rdp.rds"))
saveRDS(tax_refseq, here::here("Output/01 - Read Preprocessing - Output/tax_refseq.rds"))







samples_short <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))
samples_short <- sapply(strsplit(basename(samples_short), "_bc"), `[`,1) 
