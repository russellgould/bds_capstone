#!/usr/bin/env Rscript

library("phyloseq")
library("Biostrings")
library("dada2")

args <- commandArgs(trailingOnly = TRUE)

# NEED INPUT
# 1) silva_nr99_path: path to taxonomy table
# 2) r_objs_dir: path to directory for project R objects
# 3) sample_data: path to file with sample data in it (depth, lat/long)

######################################################
# Create phyloseq object
silva_nr99_path <- args[1]

r_objs_path <- file.path(args[2])

# load the seqtab.nochim object from dada2 output
load(file.path(r_objs_path, "seqtab_nochim.Rdata"))

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
save(ps, file = file.path(r_objs_path, "ps.Rdata"))
######################################################
