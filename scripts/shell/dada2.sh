#!/usr/bin/env bash

eval "$('/dfs/Megraw_Lab/software_downloads/conda/bin/conda' 'shell.bash' 'hook')"
conda activate bds

scripts_dir=$1
seqs_dir=$2

"$scripts_dir"/R/test.R "$seqs_dir"
