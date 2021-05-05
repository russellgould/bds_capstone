#!/usr/bin/env Rscript
# CHANGE ME to apropriate working directory
setwd("/home/skillinp/")
library(dada2); packageVersion("dada2")
# CHANGE ME to the directory containing the fastq files after unzipping.
path <- "/home/skillinp/"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
#### BE SURE TO CHECK THE THREADS, USE INTEGER TO SPECIFY HOW MANY, MULTITHREAD = TRUE IS AUTO DETERMINED
# On Windows set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, truncLen=c(240), maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
head(out)
#### CHECK THREADS HERE
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
#### CHECK THREADS HERE
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaFs[[1]]
## removing following couple of lines because they don't apply to single end reads
##mergers <- mergePairs(dadaFs, filtFs, verbose=TRUE)
### Inspect the merger data.frame from the first sample
##head(mergers[[1]])
##seqtab <- makeSequenceTable(mergers)
## Replacing previous code with one meant for forward reads only
seqtab <- makeSequenceTable (dadaFs)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#### CHECK THREADS HERE
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
# The following line had to be changed fairly significantly from the pipeline... FYI
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
# CHANGE ME AND BE SURE TO POINT TO THE ACTUAL FILE AND FILE LOCATION
#### CHECK THREADS HERE
taxa <- assignTaxonomy(seqtab.nochim, "/home/skillinp/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
head(track)
# CHANGE ME to save image wherever you want to save it
save.image("/home/skillinp/Dada2PipeImage.R")
