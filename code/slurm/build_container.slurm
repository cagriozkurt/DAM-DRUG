#!/bin/bash
#SBATCH --account=mozkurt
#SBATCH --job-name=build-sif
#SBATCH --output=/arf/scratch/mozkurt/build-%j.out
#SBATCH --error=/arf/scratch/mozkurt/build-%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks=20        # barbun minimum
#SBATCH --mem=32G

echo "Building SIF on $SLURMD_NODENAME: $(date)"
apptainer build --force $HOME/containers/scanpy-env.sif $HOME/containers/scanpy-env
echo "Build exit code: $?"
ls -lh $HOME/containers/
exit
