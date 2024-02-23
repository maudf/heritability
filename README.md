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
    │   │   ├── baselineLD.1.annot.gz
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
``
$ cd Eqtl_networks
$ sbatch run_1_compute_corescores.sh
$ sbatch run_2_compute_ld_blocks.sh
$ sbatch run_3_extract_annot_snps_in_network.sh
$ sbatch run_3_extract_annot_snps_in_network.sh
$ sbatch run_4_analyseGOenrichment.sh
$ sbatch run_5_plotGOenrichment.sh

$cd ../LDSC
$ sbatch run_1_CreateSNPAnnotations_merged.R
$ sbatch run_2_Generate_SNPLDScores_baselineLD.sh
$ sbatch run_3_Run_SNP_LDSC_score.sh
$ sbatch run_4_correlation_heritability.sh
sbatch run_5_baseline_heritability.sh
