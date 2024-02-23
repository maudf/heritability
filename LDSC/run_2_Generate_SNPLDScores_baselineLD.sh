#!/bin/bash
## Nombre de Noeuds
#SBATCH --nodes=1

## Nombre de processeur par noeud
#SBATCH --ntasks-per-node=1

## Nom du job
#SBATCH --job-name=LDscore

## Standard output and error files
#SBATCH --output=logs/%x.%j.%a.out
#SBATCH --error=logs/%x.%j.%a.err

## Quantité de RAM par noeud
#SBATCH --mem-per-cpu=20GB

## Temps limite pour lancer le job
#SBATCH --time=01-00:00:00 # days-hh:mm:ss

## Run for each chromosome
#SBATCH --array=1-22

baseDIR=$1 # Path to the analysis. E.g.: "main/"
scoresDIR=$2 # Path to the scores to analyse. E.g.: "Results/Networks/"
scoreNAME=$3 # Name of the score to analyse
ldscDIR=$4 # Path to the ldsc software
TISSUE=$5 # Tissue to

python2.7 $ldscDIR/ldsc.py \
    --l2 \
    --bfile  ${baseDIR}/LDAnnot/${scoresDIR}/${TISSUE}_${scoreNAME}/snp_1000G.EUR.hg38.${SLURM_ARRAY_TASK_ID} \
    --ld-wind-cm 1 \
    --annot ${baseDIR}/LDAnnot/${scoresDIR}/${TISSUE}_${scoreNAME}/snp_baselineLD.${SLURM_ARRAY_TASK_ID}.annot.gz \
    --out ${baseDIR}/LDResults/${scoresDIR}/${TISSUE}_${scoreNAME}/snp_baselineLD.${SLURM_ARRAY_TASK_ID}

python2.7 $ldscDIR/ldsc.py \
    --l2 \
    --bfile  ${baseDIR}/LDAnnot/${scoresDIR}/${TISSUE}_${scoreNAME}/snp_1000G.EUR.hg38.${SLURM_ARRAY_TASK_ID} \
    --ld-wind-cm 1 \
    --annot ${baseDIR}/LDAnnot/${scoresDIR}/${TISSUE}_${scoreNAME}/snp_score_annot.${SLURM_ARRAY_TASK_ID}.annot.gz \
    --out ${baseDIR}/LDScores/${scoresDIR}/${TISSUE}_${scoreNAME}/snp_score_annot.${SLURM_ARRAY_TASK_ID}
