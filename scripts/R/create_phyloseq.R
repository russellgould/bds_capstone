#!/usr/bin/env Rscript

library("phyloseq")
library("Biostrings")
library("ggplot2")

theme_set(theme_bw())

args <- commandArgs(trailingOnly = TRUE)

# NEED INPUT
# 1) seqtab.nochim: output ASV table from dada2 pipeline
# 2) silva_nr99_path: taxonomy table
# 3) sample_data: input data frame with sample data (depth, lat/long)

######################################################
# Create phyloseq object
load(args[1])
silva_nr99_path <- args[2]

if (!file.exists(silva_nr99_path)) {
  stop("Check that the silva_nr99_v138.1_train_set.fa.gz path is correct!")
}

samp_data <- as.data.frame(read.csv(args[3])[, c(2, 4, 5, 6, 7, 8)])
row.names(samp_data) <- samp_data[, 1]
samp_data <- samp_data[,2:6]

taxa <-
  assignTaxonomy(seqtab.nochim, silva_nr99_path, multithread = TRUE)
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(samp.data)
)

# rename ASVs to have a short label and save sequences in ps object
dna <- DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# save phyloseq object for posterity
save(ps, file = file.path(objs_path, "ps.Rdata"))
######################################################
