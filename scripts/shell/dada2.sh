#!/usr/bin/env bash

# our lab defines these, but they interfere with conda paths
unset R_LIBS
unset R_LIBS_USER

eval "$('/dfs/Megraw_Lab/software_downloads/conda/bin/conda' 'shell.bash' 'hook')"
conda activate bds_cgrb

scripts_dir=$1
seqs_dir=$2

"$scripts_dir"/R/test.R "$seqs_dir"
