#!/usr/bin/env Rscript

library("phyloseq")
library("Biostrings")
library("ggplot2")

theme_set(theme_bw())

args <- commandArgs(trailingOnly = TRUE)

load(args[1])
silva_nr99_path <- "/nfs3/PHARM/David_Lab/MB599/2020_2021/group1/tax/silva_nr99_v138.1_train_set.fa.gz"

######################################################
# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, silva_nr99_path, multithread = TRUE)

ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa)
)
################# The below two lines were the last two from Maude's DADA2 wizard#################
sample_sum_df <- data.frame(depth = sample_sums(ps))
colnames(sample_sum_df) <- c("depth")
############################### They aren't in the online tutorial###############################
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu / sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")
top20 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU / sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
