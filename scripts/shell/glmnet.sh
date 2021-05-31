#!/usr/bin/env bash

# make sure these don't interfere with conda
unset R_LIBS
unset R_LIBS_USER

# this activates the conda environment
eval "$('/dfs/Megraw_Lab/software_downloads/conda/bin/conda' 'shell.bash' 'hook')"
conda activate bds_cgrb

# variables taken as input from submission script
scripts_dir=$1
ps_obj_dir=$2
output_dir=$3

# actually call the R script
"$scripts_dir"/R/glmnet.R "$ps_obj_dir" "$output_dir"
