#!/usr/bin/env bash

# make sure these don't interfere with conda
unset R_LIBS
unset R_LIBS_USER

# this activates the conda environment
eval "$('/dfs/Megraw_Lab/software_downloads/conda/bin/conda' 'shell.bash' 'hook')"
conda activate bds_cgrb

# variables taken as input from submission script
scripts_dir=$1
silva_dir=$2
r_objs_dir=$3
sample_data=$4

# actually call the R script
"$scripts_dir"/R/create_phyloseq.R "$silva_dir" "$r_objs_dir" "$sample_data"
