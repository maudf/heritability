# General information
This series of scripts can be used to replicate the analyses and plots from Stone et al., 2024 ().

# Required setting
Before you start, please organise your work directory as follow, and check that you have installed all the required softwares and R packages.

## Filesystem
You must start with a file system organised as follow:

```
main
├── Data
    ├── EQTL
    │   ├── Adipose_Subcutaneous_eqtls.Rds
    │   └── ...
    ├── GENOMES
    │   ├── Baseline
    │   │   ├── snp_baselineLD.1.annot.gz
    │   │   └── ...
    │   ├── Plink
    │   │   ├── Seq.1.bed
    │   │   ├── Seq.1.bim
    │   │   ├── Seq.1.fam
    │   │   ├── Seq.1.frq
    │   │   └── ...
    │   ├── Weights
    │   │   ├── weights.1.l2.M
    │   │   ├── weights.1.l2.M_5_50
    │   │   ├── weights.1.l2.ldscore.gz
    │   │   └── ...
    └── SUMSTATS
        ├── PASS_Alzheimers_Jansen2019.sumstats.gz
        └── ...
```

The scripts will then create other folders and files with additional annotation information and results.

## Required softwares
plink
plink2
LDSC

## Required R packages
### General packages
require("getopt")
require("scriptName")

### Data management packages
require("data.table")
require("tidyverse")
require("stringr")
require("magrittr")

### Graphic plots
require("scales")
require("ggplot2")
require("gplots")
require("plyr")
require("RColorBrewer")
require("dendextend")
require("wordcloud")
require("wesanderson")

### Network analysis packages
require("netZooR")
require("Matrix")
require("igraph")


### Gene Ontology enrichment analysis packages
require("topGO")
require("org.Hs.eg.db")

# Execute the scripts
## Build eQTL networks and analyse their structure
```bash
$ cd Eqtl_networks
$ sbatch run_1_compute_corescores.sh
$ sbatch run_2_compute_ld_blocks.sh
$ sbatch run_3_extract_annot_snps_in_network.sh 
$ sbatch run_4_analyseGOenrichment.sh
$ sbatch run_5_plotGOenrichment.sh
$ cd ..
```
## Compute LDSC scores for multi-tissue analysis
Here are the instructions for a threshold of 0.75. This example includes 2 traits because it showcases the trait heritability correlation analysis.

```bash
$ mainpath="~/main"
$ ldscpath="/usr/local/bin/ldsc"
$ trait1="PASS_Alzheimers_Jansen2019"
$ trait2="PASS_BreastCancer"
$ tissue="Adipose_subcutaneous"
$ cd LDSC
$ sbatch run_1_CreateSNPAnnotations_merged.R  0.75
$ sbatch run_2_Generate_SNPLDScores_baselineLD.sh ${mainpath}/Results/ 0.75/ all_scores ${ldscpath} ${tissue}
$ sbatch run_3_Run_SNP_LDSC_score.sh ${mainpath}/Results/ 0.75/ all_scores ${mainpath}/Data/SUMSTATS ${mainpath}/Data/GENOMES/Weights/ ${ldscpath} ${tissue} ${trait1} 
$ sbatch run_3_Run_SNP_LDSC_score.sh ${mainpath}/Results/ 0.75/ all_scores ${mainpath}/Data/SUMSTATS ${mainpath}/Data/GENOMES/Weights/ ${ldscpath} ${tissue} ${trait2} 
```

## Compute LDSC genetic correlation and heritability between traits
$ mainpath="~/main"
$ ldscpath="/usr/local/bin/ldsc"
$ traits=("PASS_Alzheimers_Jansen2019" "PASS_BreastCancer" "PASS_Type_2_Diabetes")
$ for trait1 in ${traits[0..11]};
do
  sbatch run_5_baseline_heritability.sh  ${mainpath}/ ${ldscpath} ${trait1}  ;
  for trait2 in ${traits[1..12]}; 
  do
    sbatch run_4_correlation_heritability.sh ${mainpath}/ ${ldscpath} ${trait1}  ${trait2} ;
  done;
done;
sbatch run_5_baseline_heritability.sh  ${mainpath}/ ${ldscpath} ${traits[12]} ;

## Compute LDSC scores for multi-tissue analysis
Here are the instructions for a multiple threshold analysis for one trait and one tissue (Whole_Blood)

```bash
$ mainpath="~/main"
$ ldscpath="/usr/local/bin/ldsc"
$ trait1="PASS_Alzheimers_Jansen2019"
$ tissue="Whole_Blood"
$ scorename="all_scores"
$ cd LDSC
$ for thres in (0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95);
do
  sbatch run_1_CreateSNPAnnotations_merged.R ${thres} ;
  sbatch run_2_Generate_SNPLDScores_baselineLD.sh ${mainpath}/Results/ ${thres}/ ${scorename} ${ldscpath} ${tissue};
  sbatch run_3_Run_SNP_LDSC_score.sh ${mainpath}/Results/ ${thres}/ ${scorename} ${mainpath}/ata/SUMSTATS ${mainpath}/Data/GENOMES/Weights/ ${ldscpath} ${tissue} ${trait1};
  sbatch run_5_baseline_heritability.sh  ${mainpath}/Results/ ${thres}/ ${mainpath}/Data/SUMSTATS/ ${mainpath}/Data/GENOMES/Weights/ ${ldscpath} ${trait1} ;
done
```
## Extract correlation and heritability values
```bash
$ mainpath="~/main"
$ grep "^p1" e$(ls ${mainpath}/Results/LDResults/Correlation/correlation.* | head -1) >${mainpath}/Results/LDResults/Correlation/summary_correlation_traits.txt
$ for f in $(ls ${mainpath}/Results/LDResults/Correlation/correlation.*) ;
do
  grep "^${mainpath}/" ${f} | sed -e 's!/[^ ]*/PASS_!!g' | sed -e 's/.sumstats.gz//g' >>${mainpath}/Results/LDResults/Correlation/summary_correlation_traits.txt ;
done
```
## Extract heritability values
```bash
$ mainpath="~/main"
$ echo "Trait h2 se.h2" >${mainpath}/Results/LDResults/summary_heritability.txt
$ for f in $(ls ${mainpath}/Results/LDResults/Heritability/heritability.*log) ;
do
 b=$(grep "Total Observed scale h2" $f | sed -e 's/Total Observed scale h2: //g' | sed -e 's/[\\(\\)]//g')
 a=$(echo $f | sed -e 's!^.*PASS_!!g' | sed -e 's/.log//g'
 echo "$a $b" >>${mainpath}/Results/LDResults/summary_heritability.txt
done
```
