#!/bin/bash

## Nombre de Noeuds
#SBATCH --nodes=1

## Nombre de processeur par noeud
#SBATCH --ntasks-per-node=1

## Nom du job
#SBATCH --job-name=SNPAnnot

## Standard output and error files
#SBATCH --output=logs/%x.%j.%a.out
#SBATCH --error=logs/%x.%j.%a.err

## Quantité de RAM par noeud
#SBATCH --mem-per-cpu=20GB

## Temps limite pour lancer le job
#SBATCH --time=00-23:55:55 # days-hh:mm:ss

## Run for each chromosome
#SBATCH --array=1-22
thres=$1

Rscript code/1_CreateSNPAnnotations_merged.R \
-d main/ \
-a Data/Annotations/SNPs_network.mapping \
-b Data/GENOMES/Plink/ \
-n Results/Networks/ \
-l Data/GENOMES/Baseline/ \
-c ${SLURM_ARRAY_TASK_ID} \
-t $thres \
-o Results/LDAnnot/
