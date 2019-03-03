#!/bin/bash

# Replace MY_PI_SUNetID_or_Project_ID with the PI/project to be charged.
#SBATCH --account=cdbustam

# Set job time to 1 day.
#SBATCH --time=4-00:00:00

# Set a name for the job, visible in `squeue`
#SBATCH --job-name="Create Bootstrapped Files."
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=kimberly.mcmanus@gmail.com

# The following settings are optimal for *most* software, we want one task 
# to have one or more cores available for that task to fork or use threads.
# One node.
#SBATCH --nodes=1
# One task
#SBATCH --ntasks=1
# One CPU/core per task
#SBATCH --cpus-per-task=1

# There are to ways to specify memory, --mem= and --mem-per-cpu=
# --mem is usually enough since total job memory is easy to specify 
# this way.
# 5GB of RAM
#SBATCH --mem=5G

Rscript ../../scripts/R/make_input_files.R
