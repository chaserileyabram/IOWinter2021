#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=/home/abram/AnnualReviewsHA/output/testout.out
#SBATCH --error=/home/abram/AnnualReviewsHA/output/testerr.err
#SBATCH --partition=broadwl
#SBATCH --array=1-10
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=6000


module load julia
julia testjulia.jl



