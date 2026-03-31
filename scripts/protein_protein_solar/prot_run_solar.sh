#!/bin/bash
#SBATCH --job-name=solar_bivariate_all_pairs
#SBATCH -o /users/k25046756/Twins_Project/logs/protein_protein_solar/solar_%A_%a.out
#SBATCH -e /users/k25046756/Twins_Project/logs/protein_protein_solar/solar_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=6G
#SBATCH --array=1-200
#SBATCH --cpus-per-task=2

#Ensure that SOLAR is set up/can be ran on bash already 
export PATH=/users/k25046756/Twins_Project/software/solar-eclipse-9.0.1-static-Linux:$PATH

#Set up file paths
solardir=/users/k25046756/Twins_Project/software/solar-eclipse-9.0.1-static-Linux

scriptdir=/users/k25046756/Twins_Project/scripts/protein_protein_solar/solar_chunks

resultsdir=/users/k25046756/Twins_Project/results/protein_protein_solar

cd $resultsdir

#To prevent overwriting in .out file - need a separate directories for each batch
tempdir=$resultsdir/job_chunk_${SLURM_ARRAY_TASK_ID}
mkdir -p $tempdir
cd $tempdir

BATCH_FILE=$scriptdir/bivariate_chunk_${SLURM_ARRAY_TASK_ID}.txt

solar < $BATCH_FILE

