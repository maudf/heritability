#!/bin/bash

## Nombre de Noeuds
#SBATCH --nodes=1

## Nombre de processeur par noeud
#SBATCH --ntasks-per-node=1

## Standard output and error files
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err

## Quantité de RAM par noeud
#SBATCH --mem-per-cpu=30G

#SBATCH --time=3-00:00:00 # days-hh:mm:ss

### Set variables

Rscript extract_snps_annot.R -d main/ -i Results/Networks/ -o Data/Annotations/SNPs_network 
