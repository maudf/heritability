#!/bin/bash

## Nombre de Noeuds
#SBATCH --nodes=1

## Nombre de processeur par noeud
#SBATCH --ntasks-per-node=1

## Nom du job
#SBATCH --job-name=LDBLOCKS

## Nom des fichiers de sorties standard et erreur
#SBATCH --output=logs/%x.%j.%a.out
#SBATCH --error=logs/%x.%j.%a.err

## Quantité de RAM par noeud
#SBATCH --mem-per-cpu=4GB

## Temps limite pour lancer le job
#SBATCH --time=00-10:00:00 # days-hh:mm:ss

## Run for each chromosome
#SBATCH --array=1-22

infile=$1 # Root name of the bed/bim/fam file without chromosome number
outfile=$2 # Root name of the plink output file without chromosome number

# Make sure that the input file are properly sorted
plink --bfile $infile.${SLURM_ARRAY_TASK_ID} --make-bed --out $outfile.${SLURM_ARRAY_TASK_ID}

# Compute LD blocks
plink --bfile $outfile.${SLURM_ARRAY_TASK_ID} --blocks no-pheno-req --blocks-max-kb 2000 --out $outfile.${SLURM_ARRAY_TASK_ID}
