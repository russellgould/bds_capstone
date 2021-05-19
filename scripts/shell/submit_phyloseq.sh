#!/usr/bin/env bash

# make sure this path is correct!
#   everything is relative to this
group1_dir=/nfs3/PHARM/David_Lab/MB599/2020_2021/group1

# where is the code located?
scripts_dir=/raid1/home/bpp/gouldru/bds/scripts

# where is the silva_nr99 file?
silva_dir="$group1_dir"/tax/silva_nr99_v138.1_train_set.fa.gz

# where do we want the output saved?
#  plots, Rdata, etc. are saved here
r_objs_dir="$group1_dir"/dada2_FilteredFullSet_output/r_objects

# where is the map file with sample data?
sample_data="$group1_dir"/DeepSeaSurveyTrialData.csv

# these machines have had problems in the past, so we submit to the cluster
#   and avoid using them
queue="*@!(symbiosis*|galls*)"

# submit jobscript
#   -P is # processors
#   -m is amount of memory
#   -q is the queue we submit to
#   -r is where sge saves log files (stderr, stdout)
echo "bash create_phyloseq.sh $scripts_dir $silva_dir $r_objs_dir" "$sample_data" |
    SGE_Array \
        -P 10 \
        -m 50G \
        -q "$queue" \
        -r sge_logs_phyloseq
