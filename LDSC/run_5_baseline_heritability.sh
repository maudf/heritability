#!/bin/bash
## Nombre de Noeuds
#SBATCH --nodes=1

## Nombre de processeur par noeud
#SBATCH --ntasks-per-node=1

## Nom du job
#SBATCH --job-name=HERIT

## Standard output and error files
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err

## Quantité de RAM par noeud
#SBATCH --mem-per-cpu=50GB

## Temps limite pour lancer le job
#SBATCH --time=00-23:00:00 # days-hh:mm:ss

baseDIR=$1 # Path to the analysis. E.g.: "~/Documents/Heritability_eqtls_networks/"
ldscDIR=$2 # Path to the ldsc software
SUMSTAT=$3 # Trait to analyse


echo $SUMSTAT
mkdir -p ${baseDIR}/Results/LDResults/Heritability/
python2.7 $ldscDIR/ldsc.py \
 --h2 ${baseDIR}/Data/SUMSTATS/${SUMSTAT}.sumstats.gz \
 --ref-ld-chr ${baseDIR}/Data/GENOMES/Baseline/snp_baselineLD. \
 --frqfile-chr ${baseDIR}/Data/GENOMES/Plink/snp_1000G.EUR.hg38. \
 --w-ld-chr ${baseDIR}/Data/GENOMES/Weights/weights.hm3_noMHC. \
 --overlap-annot \
 --out ${baseDIR}/Results/LDResults/Heritability/heritability.${SUMSTAT}
