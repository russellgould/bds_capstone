#!/usr/bin/env Rscript

library("phyloseq")
library("ggplot2")
library("plyr")

# load("~/large_files/bds/ps.Rdata")

GP <- GlobalPatterns

# Remove OTUs that do not show appear more than 5 times in more than half the samples
wh0 <- genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
GP1 <- prune_taxa(wh0, GP)

# Transform to even sampling depth.
GP1 <- transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

# Keep only the most abundant five phyla.
phylum.sum <- tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla <- names(sort(phylum.sum, TRUE))[1:5]
GP1 <- prune_taxa((tax_table(GP1)[, "Phylum"] %in% top5phyla), GP1)

# Let’s start by plotting just the OTUs, and shading the points by Phylum.
GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 <- plot_ordination(GP1, GP.ord, type="taxa", color="Phylum", title="taxa")
print(p1)
p1 + facet_wrap(~Phylum, 3)

# Next, let’s plot only the samples, and shade the points by “SampleType” while also
#   modifying the shape according to whether they are human-associated. There are a few
#   additional ggplot2 layers added to make the plot even nicer…
p2 <- plot_ordination(GP1, GP.ord, type="samples", color="SampleType", shape="human")
p2 + geom_polygon(aes(fill=SampleType)) + geom_point(size=5) + ggtitle("samples")
