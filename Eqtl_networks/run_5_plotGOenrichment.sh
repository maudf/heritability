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

Rscript enrichment_cluster_topGO.R \
-d main/ \
-t Adipose_Subcutaneous \
-a Data/Annotations/tissues_annotation.Rds \
-m 142 \
-p 0.01 \
-i Results/GOanalysis/
