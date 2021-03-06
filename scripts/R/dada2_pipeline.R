#!/usr/bin/env Rscript

library("dada2")
library("ggplot2")
# packageVersion("dada2")

args <- commandArgs(trailingOnly = TRUE)

input_seqs_path <- file.path(args[1])
output_path <- file.path(args[2])

dir.create(output_path)
plots_path <- file.path(output_path, "plots")
stats_path <- file.path(output_path, "analysis_stats")
objs_path <- file.path(output_path, "r_objects")
dir.create(plots_path)
dir.create(stats_path)
dir.create(objs_path)

######################################################
# read samples from files

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <-
  sort(list.files(input_seqs_path, pattern = "_R1_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)

# save plots of sample quality for all samples
quality_plots_path <- file.path(plots_path, "qualities")
dir.create(quality_plots_path)
for (i in 1:length(fnFs)) {
  plotQualityProfile(fnFs[i])
  ggsave(file.path(quality_plots_path, paste(sample.names[i], ".pdf", sep = "")))
}

# Place filtered files in filtered_samples/ subdirectory
filtFs <-
  file.path(output_path, "filtered_samples", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

######################################################
# Filter and trim
out <- filterAndTrim(
  fnFs,
  filtFs,
  maxN = 0,
  maxEE = 2,
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)
write.table(out, file.path(stats_path, "filter_numbers.tsv"), quote = FALSE, sep = "\t")

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
write.table(track, file.path(stats_path, "all_steps_numbers.tsv"), quote = FALSE, sep = "\t")
