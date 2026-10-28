# Author: Sunni Patton
# Last edited: 10/20/2026
# Title: Read preprocessing full-length

## Resources used ====
## Resources used for pipeline
# https://benjjneb.github.io/LRASManuscript/LRASms_HMP.html
## Resources used for error calculations
# https://github.com/benjjneb/dada2/issues/1307
# https://github.com/benjjneb/dada2/issues/1892

## Load libraries ====
library(dada2)

## Set seed ====
set.seed(123)

## Remove primers ====
path <- here::here("Data - raw")
samples_long <- sort(list.files(path, pattern = "fastq", full.names = TRUE))

## Plot quality profile 
plotQualityProfile(path, n = 5e8, aggregate = TRUE)

F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"

for (i in samples_long) {
  filepath <- paste0(sapply(strsplit(basename(i), ".fastq"), `[`,1), "_rmprimer", ".fastq.gz")
  removePrimers(i, filepath, primer.fwd=F27, primer.rev=rc(R1492), orient=TRUE, verbose=TRUE, allow.indels = TRUE) 
  print(paste("Done with sample", basename(i)))
}

path <- here::here("Data - raw/Primers removed")
hist(nchar(getSequences(path)), 100) # Mostly around 1500 which we would expect; some very short (~50), some very long (>3000)

plotQualityProfile(path, n = 5e8, aggregate = TRUE)

path_old <- here::here("Data - raw")
hist(nchar(getSequences(path_old)), 100) # Looks like removing primers got rid of a lot of reads that were also around 1500 bp (maybe slightly shorter)

## Filter sequences ====
path <- here::here("Data - raw/Primers removed")
names <- sort(list.files(path, pattern = ".gz", full.names = FALSE))
filt <- paste0(sapply(strsplit(basename(names), "_rmprimer.fastq.gz"), `[`,1), "_filtered", ".fastq.gz")

track <- filterAndTrim(path, filt, minQ=2, minLen=1200, maxLen = 1600, maxN=0, rm.phix= FALSE, maxEE = 2, verbose = TRUE)

path <- here::here("Data - raw/Primers removed/Filtered")
plotQualityProfile(path, n = 5e8, aggregate = TRUE)

## Dereplicate reads ====
derep <- derepFastq(here::here("Data - raw/Primers removed/Filtered"), verbose = TRUE)
## Learn error ====
# Code from Callahan pipeline does not work well here. Will technically run, but plotErrors looks horrible
# Got random code from a GitHub discussion that is meant for other platforms (not necessarily PacBio), and it seems to look a lot better

suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(dplyr))

# THIS IS EXPERIMENTAL
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

errs <- learnErrors(derep, 
                    nbases = 1e9,
                    errorEstimationFunction=loessErrfun_mod4, 
                    BAND_SIZE=32, 
                    multithread=36,
                    verbose=2) # maybe took 10-15 min
## ended up converging after 7 rounds

plotErrors(errs, nominalQ=TRUE) # look for points falling along the line and decreasing error rates as quality score increases
# looks a lot better

## Denoise ====
dd <- dada(derep, err=errs, BAND_SIZE=32, multithread=TRUE) # seconds

## Make sequence table ====
seq_table <- makeSequenceTable(dd) 
table(nchar(getSequences(seq_table))) 
dim(seq_table) # 648 ASVs
sum(seq_table) #1222273
rowSums(seq_table)

## Remove chimeras ====
seq_table_nochim <- removeBimeraDenovo(seq_table, method="consensus", multithread=TRUE, verbose=TRUE) #Identified 43 bimeras out of 648 input sequences.
dim(seq_table_nochim) # 605 contigs
sum(seq_table) - sum(seq_table_nochim) # lost 37971 reads
rowSums(seq_table_nochim)
## Assign taxonomy ====
taxa <- assignTaxonomy(seq_table_nochim, here::here("Training Sets/silva_nr99_v138.1_train_set.fa.gz"), multithread=TRUE)

## Remove off-target sequences ====
### New sequence table (no eukarya)
is.euk <- taxa[,"Kingdom"] %in% "Eukaryota"
seq_table_noeuk <- seq_table_nochim[,!is.euk]
dim(seq_table_noeuk) # 591 contigs; 14 identified as eukaryote
sum(seq_table_nochim) - sum(seq_table_noeuk) # 14 contigs associated with 1945 reads
table(nchar(getSequences(seq_table_noeuk))) 
rowSums(seq_table_noeuk)
sum(seq_table_noeuk) # 1182357
saveRDS(seq_table_noeuk, here::here("seq_table_noeuk_full.rds"))

### New taxonomy table (no eukarya)
taxonomy_noeuk <- taxa[!is.euk,]
dim(taxonomy_noeuk)
saveRDS(taxonomy_noeuk, here::here("taxonomy_noeuk_full.rds"))

### New sequence table (no chloro)
is.chloro <- taxonomy_noeuk[,"Order"] %in% "Chloroplast"
seq_table_nochloro <- seq_table_noeuk[,!is.chloro]
dim(seq_table_nochloro) # 580 contigs; 11 identified as chloroplast
sum(seq_table_noeuk) - sum(seq_table_nochloro) # 11 contigs associated with 3627 reads
table(nchar(getSequences(seq_table_nochloro))) 
rowSums(seq_table_nochloro)
sum(seq_table_nochloro) # 1178730
saveRDS(seq_table_nochloro, here::here("seq_table_nochloro_full.rds"))

### New taxonomy table (no chloro)
taxonomy_nochloro <- taxonomy_noeuk[!is.chloro,]
dim(taxonomy_nochloro)
saveRDS(taxonomy_nochloro, here::here("taxonomy_nochloro_full.rds"))

### New sequence table (no NA kingdom)
is.NA <- taxonomy_nochloro[,"Kingdom"] %in% NA
seq_table_noNA <- seq_table_nochloro[,!is.NA]
dim(seq_table_noNA) # 571 contigs; 9 identified as NA at Kingdom level
## only one ASV from this group had a blast hit (to chocolate vine plant)
sum(seq_table_nochloro) - sum(seq_table_noNA) # 9 contigs associated with 294 reads
table(nchar(getSequences(seq_table_noNA))) 
rowSums(seq_table_noNA)
sum(seq_table_noNA) # 1178436
saveRDS(seq_table_noNA, here::here("seq_table_noNA_full.rds"))

### New taxonomy table (no NA kingdom)
taxonomy_noNA <- taxonomy_nochloro[!is.NA,]
dim(taxonomy_noNA)
saveRDS(taxonomy_noNA, here::here("taxonomy_noNA_full.rds"))

### New sequence table (no NA phylum)
is.NA_phy <- taxonomy_noNA[,"Phylum"] %in% NA
seq_table_noNA_phy <- seq_table_noNA[,!is.NA_phy]
dim(seq_table_noNA_phy) # 508 contigs; 63 identified as NA at phylum level
sum(seq_table_noNA) - sum(seq_table_noNA_phy) # 63 contigs associated with 3486 reads
table(nchar(getSequences(seq_table_noNA_phy))) 
rowSums(seq_table_noNA_phy)
sum(seq_table_noNA_phy) # 1174950
saveRDS(seq_table_noNA_phy, here::here("seq_table_noNA_phy_full.rds"))

### New taxonomy table (no NA phylum)
taxonomy_noNA_phy <- taxonomy_noNA[!is.NA_phy,]
dim(taxonomy_noNA_phy)
saveRDS(taxonomy_noNA_phy, here::here("taxonomy_noNA_phy_full.rds"))

## There are a lot of ASVs that are showing up as mitochondria, and are associated with ~45,000 reads
## manually going through to see if any are potentially bacteria

tax_test <- taxonomy_noNA_phy
tax_test <- as.data.frame(tax_test)
tax_test$Seq <- paste0(rownames(tax_test))
tax_test$Family[tax_test$Seq == "GGTGAACGCTAGTAGTATGCTTTACACATGCGAGTTGAACGAATACGTTTCGTGGCGAACGGGTGAGGAATATATAGAATTTTACCTTTTAATTCAGAATAAACAAATTCTTTTTTGAATGATTTTTAGATGTTAATACTGAATATCATAAAAGCTATCAACGAAACATTTTATAAATAGTTGTTAAAAGATAAATCTATATAAGATTAGGTAGTTGGTAATTTAACAGATTTACCAAGCCAATGATCTTTAGTTGATCTGAGAGGATGATCAACCACACTGGAACTGAGACAAGGTCCAGGCTAAATTTTTGGCCAGCAGTGAGGAATATTGGACAATGAACGAAAGTTTGATCCAGCAATACTATTTGAGTGATGACGGCTATTTTGGTTGTAAAACTCTTTCGGCAGTAATGATAATGACATTATCTGAAGAAGTTTAGTCCCGGCTAATCTCGTGCCAGCAGCCGCGGTAATACGGGAGGGGCAAGCGTTATCCGGAATAATTGGGCGTAAAGGGTCCGTAGACTGTAATGTAAGTTATTCGTTAATTCTTAAAGCTTAACTTTAAATTCGCGTATAAGACTACTTTACTTGAGTTTTTTTCAGAAAAGCAGAATTTTATAAGAAGGGGTGATACCTGTAGATATATAAAGGAATATCATTAGCGAAAGCAGCTTTTTAGTAAAAACTGACGTTGAGGGACGAAAGTGTAGGGAGCGAACAGGATTAGATACCCTGGTAGTCTACACTGTAAATTCTGAGTGCTGTAATTTAATTTTTTAATTAGATTCTTTAAGCTAACGCCATAAGCACTCCGCCTGAGAAGTACTGTCGCAAGATTGGAACTCAAAGGAATTGACGGGGGCCTACTACTAGTGGTGGAGCATGTGGTTTAATGCGATAATCCGCGCAAAACCTTACCAATTTTTGACATAAAATTTTTCATAAAATCAGATAAGATGTTATAGTGCAAGATTTTACAGGTGTTGCACGACTGTCGTCAGCTCGTGTCGTGAGATGTTTGGTTAATTCCTTTAACGAGCGCAACTCTCATTTTTAATGATCAACGCATTTAATTTATAGTTAAAATGCTTTAAACATACTGTCAATGATAAATTGAAGGAAGTGAGAATTAAGTCAAGTCAATGTGGCCCTTATAAATTGGGCTACACACGTGCTACATTGAAAAATACAAAGAGTAACGAAAAATTGAATTTAAGTTAATCTGTTAAAATTTTTCTTAGTTCAGATTGTAGACTGCAACTCGTCTACATGAAGTCGGAATCACTAGTAATCGCAAATCAGCATGTTGCGGTGAATATGTTACTTGGCCTTGTACACACCGCCCGTCACATGCAGAGAGTTAATTTTGCTTGAAAATTTAATACAAATGTAAATTCATTGAAAAATTAGCAATTTGGATG"] <- "Mitochondria_uncultured_bact"
## uncultured alphaproteobacteria E = 0.0, %ID = 86.21%, query cover = 100%

tax_test$Family[tax_test$Seq == "GGTGAACGCTAGTAGTATGCTTTACACATGCGAGTTGAACGAATACGTTTCGTGGCGAACGGGTGAGGAATATATAGAATTTTACCTTTTAATTCAGAATAAACAAATTCTTTTTTGAATGATTTTTAGATGTTAATACTGAATATCATAAAAGCTATCAACGAAACATTTTATAAATAGTTGTTAAAAGATAAATCTATATAAGATTAGGTAGTTGGTAATTTAACAGATTTACCAAGCCAATGATCTTTAGTTGATCTGAGAGGATGATCAACCACACTGGAACTGAGACAAGGTCCAGGCTAAATTTTTGGCCAGCAGTGAGGAATATTGGACAATGAACGAAAGTTTGATCCAGCAATACTATTTGAGTGATGACGGCTATTTTGGTTGTAAAACTCTTTCGGCAGTAATGATAATGACATTATCTGAAGAAGTTTAGTCCCGGCTAATCTCGTGCCAGCAGCCGCGGTAATACGGGAGGGGCAAGCGTTATCCGGAATAATTGGGCGTAAAGGGTCCGTAGACAGTAATGTAAGTTATTCGTTAATTCTTAAAGCTTAACTTTAAATTCGCGTATAAGACTACTTTACTTGAGTTTTTTTCAGAAAAGCAGAATTTTATAAGAAGGGGTGATACCTGTAGATATATAAAGGAATATCATTAGCGAAAGCAGCTTTTTAGTAAAAACTGACGTTGAGGGACGAAAGTGTAGGGAGCGAACAGGATTAGATACCCTGGTAGTCTACACTGTAAATTCTGAGTGCTGTAATTTAATTTTTTAATTAGATTCTTTAAGCTAACGCCATAAGCACTCCGCCTGAGAAGTACTGTCGCAAGATTGGAACTCAAAGGAATTGACGGGGGCCTACTACTAGTGGTGGAGCATGTGGTTTAATGCGATAATCCGCGCAAAACCTTACCAATTTTTGACATAAAATTTTTCATAAAATCAGATAAGATGTTATAGTGCAAGATTTTACAGGTGTTGCACGACTGTCGTCAGCTCGTGTCGTGAGATGTTTGGTTAATTCCTTTAACGAGCGCAACTCTCATTTTTAATGATCAACGCATTTAATTTATAGTTAAAATGCTTTAAACATACTGTCAATGATAAATTGAAGGAAGTGAGAATTAAGTCAAGTCAATGTGGCCCTTATAAATTGGGCTACACACGTGCTACATTGAAAAATACAAAGAGTAACGAAAAATTGAATTTAAGTTAATCTGTTAAAATTTTTCTTAGTTCAGATTGTAGACTGCAACTCGTCTACATGAAGTCGGAATCACTAGTAATCGCAAATCAGCATGTTGCGGTGAATATGTTACTTGGCCTTGTACACACCGCCCGTCACATGCAGAGAGTTAATTTTGCTTGAAAATTTAATACAAATGTAAATTCATTGAAAAATTAGCAATTTGGATG"] <- "Mitochondria_uncultured_bact"
## uncultured alphaproteobacteria E = 0.0, %ID = 86.15%, query cover = 100%

as.matrix(tax_test) -> taxonomy_noNA_phy

### New sequence table (no NA mito)
is.mito <- taxonomy_noNA_phy[,"Family"] %in% "Mitochondria"
seq_table_nomito <- seq_table_noNA_phy[,!is.mito]
dim(seq_table_nomito) # 337 contigs; 171 identified as mito
sum(seq_table_noNA_phy) - sum(seq_table_nomito) # 63 contigs associated with 44933 reads
table(nchar(getSequences(seq_table_nomito))) 
rowSums(seq_table_nomito)
sum(seq_table_nomito) # 1130017
saveRDS(seq_table_nomito, here::here("seq_table_nomito_full.rds"))

### New taxonomy table (no  mito)
taxonomy_nomito <- taxonomy_noNA_phy[!is.mito,]
dim(taxonomy_nomito)
saveRDS(taxonomy_nomito, here::here("taxonomy_nomito_full.rds"))









## didn't end up doing this bc filter and trim already did 1200-1600
### Keep contigs within desired size range
seq_table<-seq_table[,nchar(colnames(seq_table)) %in% 1420:1481] 
## from snail data it looks like the range wasn't small enough to capture aquarickettsia (1421); would expect it to be 1460 with the primer sequences and 1420 without them


