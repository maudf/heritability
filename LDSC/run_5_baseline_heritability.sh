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
scoresDIR=$2 # Path to the scores to analyse. E.g.: "Results/Networks/Weighted/"
sumstatDIR=$3 # Path to the sumstats to analyse. E.g.: "Data/UKKB"
weightfileDIR=$4 # Path to the files with hm3 weights "Data/1000G.EUR.QC/1000G_Phase3_weights_hm3_no_MHC"
ldscDIR=$5 # Path to the ldsc software
SUMSTAT=$6 # Trait to analyse


echo $SUMSTAT
python2.7 $ldscDIR/ldsc.py \
 --h2 ${baseDIR}/${sumstatDIR}/${SUMSTAT}.sumstats.gz \
 --ref-ld-chr ${baseDIR}/${scoresDIR}/snp_baselineLD. \
 --frqfile-chr ${baseDIR}/${scoresDIR}/snp_1000G.EUR.hg38. \
 --w-ld-chr ${baseDIR}/${weightfileDIR}/weights.hm3_noMHC. \
 --overlap-annot \
 --out ${baseDIR}/${scoresDIR}/heritability.${SUMSTAT}
