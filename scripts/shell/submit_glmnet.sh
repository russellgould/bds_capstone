#!/usr/bin/env bash

# make sure this path is correct!
#   everything is relative to this
group1_dir=/nfs3/PHARM/David_Lab/MB599/2020_2021/group1

# where is the code located?
scripts_dir=/raid1/home/bpp/gouldru/bds/scripts

# which input dataset are we using?
ps_obj_dir="$group1_dir"/dada2_FilteredFullSet_output/r_objects

# these machines have had problems in the past, so we submit to the cluster
#   and avoid using them
queue="*@!(symbiosis*|galls*)"

# submit jobscript
#   -P is # processors
#   -m is amount of memory
#   -q is the queue we submit to
#   -r is where sge saves log files (stderr, stdout)
echo "bash glmnet.sh $scripts_dir $ps_obj_dir $group1_dir" |
    SGE_Array \
        -P 4 \
        -m 500G \
        -q "$queue" \
        -r sge_logs_glmnet
