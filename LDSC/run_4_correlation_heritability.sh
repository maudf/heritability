#!/bin/bash
##Partition type
#SBATCH --partition=long

## Nombre de Noeuds
#SBATCH --nodes=1

## Nombre de processeur par noeud
#SBATCH --ntasks-per-node=1

## Standard output and error files
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err

## Nom du job
#SBATCH --job-name=CORRTRAIT

## Quantité de RAM par noeud
#SBATCH --mem-per-cpu=50GB

## Temps limite pour lancer le job
#SBATCH --time=02-12:00:00 # days-hh:mm:ss

baseDIR=$1 # Path to the analysis. E.g.: "~/Documents/Heritability_eqtls_networks/"
ldscDIR=$2 # Path to the ldsc software
SUMSTAT1=$3 # First trait
SUMSTAT2=$4 # Second trait


echo $SUMSTAT
python2.7 $ldscDIR/ldsc.py \
 --rg ${baseDIR}/Data/SUMSTATS/${SUMSTAT1}.sumstats.gz,${baseDIR}/Data/SUMSTATS/${SUMSTAT2}.sumstats.gz \
 --ref-ld-chr ${baseDIR}/Data/GENOMES/Baseline/snp_baselineLD. \
 --frqfile-chr ${baseDIR}/Data/GENOMES/Plink/Seq. \
 --w-ld-chr ${baseDIR}/Data/GENOMES/Weights/weights. \
 --overlap-annot \
 --out ${baseDIR}/Results/LDResults/Correlation/correlation.${SUMSTAT1}.${SUMSTAT2}
