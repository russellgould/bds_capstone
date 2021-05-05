library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")

theme_set(theme_bw())

PS <- readRDS("~/large_files/bds/phyloseq/subset_ps")

ntaxa(PS)
nsamples(PS)

# normalize columns
PSr <- transform_sample_counts(PS, function(x) x / sum(x))
# remove columns with mean < 1e-5
PSfr <- filter_taxa(PSr, function(x) mean(x) > 1e-5, TRUE)

