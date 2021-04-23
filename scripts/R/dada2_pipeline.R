library(dada2)
# packageVersion("dada2")

path <-
  "~/large_files/bds" # CHANGE ME to the directory containing the fastq files after unzipping.
# list.files(path)

######################################################
# read samples from files

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <-
  sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# inspect read quality of forward strand
plotQualityProfile(fnFs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <-
  file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
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
head(out)

######################################################
# Learn errors
errF <- learnErrors(filtFs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

######################################################
# Apply core inference algorithm
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)

# optionally inspect output object
# dadaFs[[1]]

######################################################
# Construct ASV table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

######################################################
# Check for chimeric sequences
seqtab.nochim <-
  removeBimeraDenovo(seqtab,
    method = "consensus",
    multithread = TRUE,
    verbose = TRUE
  )
dim(seqtab.nochim)
pct_chimeric <- sum(seqtab.nochim) / sum(seqtab)

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
head(track)

######################################################
# Assign taxonomy
library(DECIPHER)
# Create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

##** CHANGE TO THE PATH OF YOUR TRAINING SET **##
load("~/tax/IDTaxa/SILVA_SSU_r132_March2018.RData")

# use all processors
ids <-
  IdTaxa(
    dna,
    trainingSet,
    strand = "top",
    processors = NULL,
    verbose = FALSE
  )
# ranks of interest
ranks <-
  c(
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)
