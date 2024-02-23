#!/bin/bash
## Nombre de Noeuds
#SBATCH --nodes=1

## Nombre de processeur par noeud
#SBATCH --ntasks-per-node=1

## Nom du job
#SBATCH --job-name=H2

## Standard output and error files
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err

## Quantité de RAM par noeud
#SBATCH --mem-per-cpu=50GB

## Temps limite pour lancer le job
#SBATCH --time=02-12:00:00 # days-hh:mm:ss

baseDIR=$1 # Path to the analysis. E.g.: "main/"
scoresDIR=$2 # Path to the scores to analyse. E.g.: "Results/Networks/"
scoreNAME=$3 # Name of the score to analyse
sumstatDIR=$4 # Path to the sumstats to analyse. E.g.: "Data//"
weightfileDIR=$5 # Path to the files with hm3 weights "Data/Genomes/Weights"
ldscDIR=$6 # Path to the ldsc software
TISSUE=$7
SUMSTAT=$8


python2.7 $ldscDIR/ldsc.py \
      --h2 ${sumstatDIR}/${SUMSTAT}.sumstats.gz \
      --ref-ld-chr ${baseDIR}/LDScores/${scoresDIR}/${TISSUE}_${scoreNAME}/snp_score_annot.,${baseDIR}/LDScores/${scoresDIR}/${TISSUE}_${scoreNAME}/snp_baselineLD. \
      --frqfile-chr ${baseDIR}/LDAnnot/${scoresDIR}/${TISSUE}_${scoreNAME}/snp_1000G.EUR.hg38. \
      --w-ld-chr ${weightfileDIR}/weights.hm3_noMHC. \
      --overlap-annot \
      --print-coefficients \
      --print-delete-vals \
      --out ${baseDIR}/LDResults/${scoresDIR}/${TISSUE}_${scoreNAME}/snp_baseline.score_annot.${SUMSTAT}
