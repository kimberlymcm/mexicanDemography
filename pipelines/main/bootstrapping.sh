# This pipeline is used to get the bootstrapping results

# Order to run:

# Get the 500 kb blocks that we will break up the bootstrapping with.
# It reads from:
# /home/kfm/kfm_projects/NA/NA_data/getCallableLength/4fold_intron_outgroup_20150826.sites.only.bed
# which I guess is a bed file limited to the relevant sites.
# [kfm@smsx10srw-srcf-d15-37 bootstrapping]$ cat /home/kfm/kfm_projects/NA/NA_data/getCallableLength/4fold_intron_outgroup_20150826.sites.only.bed | wc -l
# 8889201

mkdir -p ../../data/bootstrapping
Rscript ../../scripts/R/get_blocks.r

mkdir -p ../../data/bootstrapping/input_bootstrap
mkdir -p ../../results/bootstrapping

# Make the 1,000 bootstrapped datadict input files
# Below runs makeInputFiles.r
# It reads from: ../NA_CHB_exomes_20150826.data_dict
# It also calculates the callable genome length of each bootstrap replicate
sh ../server/create_bootstrap_input_files.sh

for i in {1..1000}; do
    sed -i $'1 i\\\nHuman\tChimp\tAllele1\tHUI\tMYA\tNAH\tTAR\tTRQ\tCHB\tAllele2\tHUI\tMYA\tNAH\tTAR\tTRQ\tCHB\tchr\tpos\tnum' \
    ../../data/bootstrapping/input_bootstrap/Bootstrap_${i}.txt
done


# Then, run dadi
# This runs: 4pop_dadi_bootstrapping.py
# Divergence rate: 0.0121
# Tb: 0.09 (even though it looks like I stored results for 0.09 and 0.1, so maybe I manually changed this)
sh ../server/run_dadi.sh

# Get the confidence intervals of the results
Rscript ../../scripts/R/get_confidence_intervals.R
