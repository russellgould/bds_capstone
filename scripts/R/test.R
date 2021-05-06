#!/usr/bin/env Rscript

# this script is just to test whether packages were installed correctly

library("dada2"); packageVersion("dada2")
# library("phyloseq"); packageVersion("phyloseq")
# library("ggplot2"); packageVersion("ggplot2")
# library("plyr"); packageVersion("plyr")

# test that arguments are correct and working
# args <- commandArgs(trailingOnly = TRUE)

args <- c("~/large_files/bds")

input_seqs_path <- file.path(args[1])
plots_path <- file.path(input_seqs_path, "plots")
stats_path <- file.path(input_seqs_path, "analysis_stats")
dir.create(plots_path)
dir.create(stats_path)

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

pct.chim <- 0.58983475983475

writeLines(sprintf("%.4f", pct.chim), file.path(stats_path, "percent_chimeric.txt"))
