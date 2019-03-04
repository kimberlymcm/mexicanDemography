#!/bin/bash

#SBATCH --account=cdbustam
#SBATCH --time=4-00:00:00
#SBATCH --job-name="Rotate_dadi_models"
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=kimberly.mcmanus@gmail.com
# The following settings are optimal for *most* software, we want one task 
# to have one or more cores available for that task to fork or use threads.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# There are to ways to specify memory, --mem= and --mem-per-cpu=
# --mem is usually enough since total job memory is easy to specify 
# this way.
#SBATCH --mem=5G

echo "Parallel task num ${num}"
date
echo "This job runs the 4 population rotating models on native mexican data."
echo "I ran on host: \$(hostname -s)"
echo "SLURM Environment is:"
env | grep "SLURM" | sort
echo "My limits are:"
ulimit -a

python2.7 ../../scripts/python/run_model_server_rotate.py
