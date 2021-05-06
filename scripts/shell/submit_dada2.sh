#!/usr/bin/env bash

# make sure this path is correct!
#   everything is relative to this
group1_dir=/nfs3/PHARM/David_Lab/MB599/2020_2021/group1

# where is the code located?
scripts_dir="$group1_dir"/code/scripts

# which input dataset are we using?
seqs_dir="$group1_dir"/SubsetForTrial

# where do we want the output saved?
#  plots, Rdata, etc. are saved here
output_dir="$group1_dir"/subset_output

# these machines have had problems in the past, so we submit to the cluster
#   and avoid using them
queue="*@!(symbiosis*|galls*)"

# submit jobscript
#   -P is # processors
#   -m is amount of memory
#   -q is the queue we submit to
#   -r is where sge saves log files (stderr, stdout)
echo "bash dada2.sh $scripts_dir $seqs_dir $output_dir" |
    SGE_Array \
        -P 50 \
        -m 250G \
        -q "$queue" \
        -r sge_logs_dada2
