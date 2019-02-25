#!/bin/bash

# Run the 4 population dadi model across all bootstrapping files
# Loop through all the bootstrap files
# I don't know if this is the best way to do this with SLURM, but
# it is the only way I saw an example for

for num in {1..1000}; do

    cat << EOF | sbatch
#!/bin/bash
#SBATCH --account=cdbustam
#SBATCH --time=4-00:00:00
#SBATCH --job-name="Run_dadi_models_${num}"
#SBATCH --mail-type=FAIL
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
echo "This job runs the 4 population dadi model on bootstrap native mexican data."
echo "I ran on host: \$(hostname -s)"
echo "SLURM Environment is:"
env | grep "SLURM" | sort
echo "My limits are:"
ulimit -a

python ../../scripts/python/4pop_dadi_bootstrapping.py --inFile ${num} --out ${num}
    
EOF
done