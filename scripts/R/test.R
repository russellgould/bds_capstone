#!/usr/bin/env Rscript

library("dada2"); packageVersion("dada2")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")

# just test that the seqs can be found
args = commandArgs(trailingOnly = TRUE)
seqs <- args[1]
list.files(path)
