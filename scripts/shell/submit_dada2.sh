#!/usr/bin/env bash

queue="*@!(symbiosis*|galls*)"

scripts_dir="/raid1/home/bpp/gouldru/bds/scripts"
seqs_dir="/nfs3/PHARM/David_Lab/MB599/2020_2021/group1/SubsetForTrial"

# submit jobscript
echo "bash dada2.sh $scripts_dir $seqs_dir" |
    SGE_Array \
        -P 1 \
        -m 5G \
        -q "$queue" \
        -r sge_logs_dada2
