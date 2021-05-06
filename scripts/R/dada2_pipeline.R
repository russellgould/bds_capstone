#!/usr/bin/env Rscript

library("dada2")
# packageVersion("dada2")

args <- commandArgs(trailingOnly = TRUE)

input_seqs_path <- file.path(args[1])
plots_path <- file.path(input_seqs_path, "plots")
stats_path <- file.path(input_seqs_path, "analysis_stats")
objs_path <- file.path(input_seqs_path, "r_objects")
dir.create(plots_path)
dir.create(stats_path)
dir.create(objs_path)

######################################################
# read samples from files

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <-
  sort(list.files(input_seqs_path, pattern = "_R1_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# inspect read quality of forward strand
pdf(file.path(plots_path, "quality_profile.pdf"))
plotQualityProfile(fnFs[1:2])
dev.off()

# Place filtered files in filtered/ subdirectory
filtFs <-
  file.path(input_seqs_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

######################################################
# Filter and trim
out <- filterAndTrim(
  fnFs,
  filtFs,
  truncLen = 210,
  maxN = 0,
  maxEE = 2,
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)
# head(out)

######################################################
# Learn errors
errF <- learnErrors(filtFs, multithread = TRUE)

pdf(file.path(plots_path, "errors.pdf"))
plotErrors(errF, nominalQ = TRUE)
dev.off()

######################################################
# Apply core inference algorithm
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)

save(dadaFs, file = file.path(objs_path, "dadaFs.Rdata"))
# optionally inspect output object
# dadaFs[[1]]

######################################################
# Construct ASV table
seqtab <- makeSequenceTable(dadaFs)
# dim(seqtab)

# Inspect distribution of sequence lengths
# table(nchar(getSequences(seqtab)))

######################################################
# Check for chimeric sequences
seqtab.nochim <-
  removeBimeraDenovo(seqtab,
    method = "consensus",
    multithread = TRUE,
    verbose = TRUE
  )
save(seqtab.nochim, file = file.path(objs_path, "seqtab_nochim.Rdata"))
pct_chimeric <- sum(seqtab.nochim) / sum(seqtab)
writeLines(sprintf("%.4f", pct_chimeric), file.path(stats_path, "percent_chimeric.txt"))

######################################################
# Track read numbers
getN <- function(x) {
  sum(getUniques(x))
}
track <-
  cbind(
    out,
    sapply(dadaFs, getN),
    rowSums(seqtab.nochim)
  )
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <-
  c(
    "input",
    "filtered",
    "denoisedF",
    "nonchim"
  )
rownames(track) <- sample.names
write.csv(track, file.path(stats_path, "read_numbers.csv"), quote = FALSE)
