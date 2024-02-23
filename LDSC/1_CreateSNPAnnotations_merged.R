###############################################
#### Step 1: Create annot files by chromosome
###############################################

### Load libraries
require("getopt")
require("scriptName")
require(data.table)
require(tidyverse)
require(stringr)
require(magrittr)
require(scales)
cat("Libraries loaded\n")

## Set variables
spec = matrix(c('help','h', 0, "logical", "Print help message",
                'directory','d', 1, "character", "Path to main folder",
                'annot','a', 1, "character", "SNP annotation file",
                'bim','b', 1, "character", "BIM file with SNP positions in BP and CM",
                'network','n', 1, "character",  "Path to Network folder with scores",
                'ldbaseannot','l', 1, "character",  "Path to LDbase annot files",
                'chr','c', 1, "character", "Chromosome name to study",
                'threshold','t', 1, "numeric", "Threshold for scores to determine 0 and 1 categories",
                'output','o', 1, "character", "Path to the output folder"
              ), byrow=TRUE, ncol=5)

opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)}

if ( is.null(opt$directory    ) ) {opt$directory    = '~/Documents/Heritability_eqtl_networks/'     }
if ( is.null(opt$annot    ) ) {opt$annot    = 'Data/Annotations/SNPs_network.mapping'     }
if ( is.null(opt$bim    ) ) {opt$bim    = 'Data/Genome/Plink/Seq.'     }
if ( is.null(opt$chr    ) ) {opt$chr    = '1'     }
if ( is.null(opt$network    ) ) {opt$network    = 'Results/Networks/'     }
if ( is.null(opt$ldbaseannot    ) ) {opt$ldbaseannot    = 'Data/Genome/Baseline/'     }
if ( is.null(opt$threshold    ) ) {opt$threshold    = 0.75  }
if ( is.null(opt$output    ) ) {opt$output    = 'Results/LDAnnot/'     }

cat("Running", current_filename(), "with options:\n",
    "directory =", opt$directory, "\n",
    'annot =',opt$annot, "\n",
    'bim =',opt$bim, "\n",
    'network =',opt$network, "\n",
    'ldbaseannot =',opt$ldbaseannot, "\n",
    'chr =',opt$chr, "\n",
    'threshold =',opt$threshold, "\n",
    "output =",opt$output, "\n",)

###Check if annotation file exists
if(!file.exists(opt$directory, opt$annot))){
  stop(paste0(opt$directory, opt$annot,
              " not found, missing snp annotation file. Exit !\n"))
}

###Check if bim file exists
if(!file.exists(opt$directory, opt$bim, opt$chr, '.bim'))){
  stop(paste0(opt$directory, opt$bim, opt$chr, '.bim',
              " not found, missing BIM file. Exit !\n"))
}

###Check if network folder exists
if(!dir.exists(opt$directory, opt$network))){
  stop(paste0(opt$directory, opt$network,
              " not found, missing network folder. Exit !\n"))
}

if(!dir.exists(opt$directory, opt$ldbaseannot))){
  stop(paste0(opt$directory, opt$ldbaseannot,
              " not found, missing LDbase annotation files folder. Exit !\n"))
}
### Create output folder if do not exists
if(!dir.exists(dirname(paste(opt$directory, opt$output, sep="/")))){
  dir.create(dirname(paste(opt$directory, opt$output, sep="/")),
  recursive = T, showWarnings = F)
}


### Define functions

readSNPPos <- function(snp_pos_file){
  if(str_detect(snp_pos_file, '.gz$')){
    snp_pos_file <- str_c('zcat ', snp_pos_file)}
  snp_pos_file %>% fread %>% as.data.frame %>% set_colnames(c('node',	'chr', 'pos'))}

################


### Get the SNP positions and annotations
cat("Loading SNP annotations at ", opt$directory, opt$annot, "\n", sep="")
bed_map <- readRDS(paste0(opt$directory, opt$annot))
bed_map <- data.table(bed_map)
setkey(bed_map,SNP)
cat("Annotations loaded\n")

cat("Loading SNP positions at ", opt$directory, opt$bim, opt$chr, '.bim', "\n", sep="")
bim_in <- fread(paste0(opt$directory, opt$bim, opt$chr, '.bim'))
names(bim_in) <- c('CHR','SNP','CM','BP','ref','alt')
cat("Positions loaded\n")


### Get the score data to analyse
cat("Getting scores files at ", opt$directory, opt$network, "...\n", sep="")
files<-list.files(path=paste0(opt$directory, opt$network), pattern = paste0( "*_scores.Rds"))
cat("List of files: \n")
print(files)
data_id <- gsub(paste0("_scores.Rds"), "", files)
cat("List of tissues:\n")
print(data_id)
## Loop through all dataset files to prepare .annot for each dataset
for (f in 1:length(data_id)){#grep("Lung", data_id):length(data_id)){
  s <- data_id[f]
  cat("Running script for dataset", s, "...\n")
  ##Loop through all scores files to prepare .annot for each dataset
  scores <- readRDS(paste0(opt$directory, opt$network,  data_id[f], "_", prefix,".Rds")))
  scores[, network := 1]
  scores=scores[,c('SNP','network', 'outdegree', 'Qik')]
  scores$outdegree.annot <- as.numeric(scores$outdegree >= quantile(scores$outdegree, threshold, na.rm=T))
  scores$corescore.annot <- as.numeric(scores$Qik >= quantile(scores$Qik, threshold, na.rm=T))
  scores[, outdegree:=NULL]
  scores[, Qik:=NULL]
  colnames(scores)=c('SNP','network', 'outdegree', 'corescore')
  setkey(scores, SNP)

  scores_mapped <- merge(scores, bed_map, by='SNP')

  ## Pull annotations for each chromosome
  scores_mapped[,chr_ind := substring(scores_mapped$CHR, 4, nchar(scores_mapped$CHR))]
  cat("Generating annotation files for chromosome", opt$chr, "...\n")

  ## Prepare annots
  cat("Loading scores...\n")
  scores_chr <- scores_mapped[CHR == opt$chr,]
  scores_chr$BP <- as.numeric(scores_chr$BP)
  key_col <- c('BP')

  ##Remove duplicates
  bed_chr_dup <- merge(bim_in, scores_chr[,c('network','outdegree.annot','corescore.annot','BP')],
                      by=key_col, all.x=T)
  set(bed_chr_dup, which(is.na(bed_chr_dup[,'network'])),7L,0)
  set(bed_chr_dup, which(is.na(bed_chr_dup[,'outdegree'])),8L,0)
	set(bed_chr_dup, which(is.na(bed_chr_dup[,'corescore'])),9L,0)
  bed_chr_dedup <- bed_chr_dup[,.SD[which.max(corescore)], by=SNP]
  bed_chr_dedup1 <- bed_chr_dedup[,.SD[which.max(corescore)], by=BP]
  bed_chr_dedup2 <- bed_chr_dedup1[,.SD[which.max(outdegree)], by=SNP]
  bed_chr <- bed_chr_dedup2[,.SD[which.max(outdegree)], by=BP]

  ## Get relevant annot files and save
  dir.create(paste0(opt$directory, opt$output, threshold, "/", data_id[f], '_all_scores/'),
            recursive = T, showWarnings = F) # create output file

  ## Prepare frq, bim/bed/fam
  cat("Preparing bed, bim, fam, and frq files...\n")
  write.table(data.frame(bed_chr$SNP, stringsAsFactors=F),
                  file=paste0(opt$directory, opt$output, threshold, "/",
                  data_id[f], '_all_scores/snp_kept_', opt$chr, '.snplist'),
                  row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")

  #Subset plink files by tissue/chrom
  cat("subsetting relevant SNPs...\n")
  system(paste0('plink2 --bfile ', opt$directory, opt$bim, opt$chr, ' \\
           --extract ', opt$directory, opt$output, threshold, "/", data_id[f],
           '_all_scores/snp_kept_', opt$chr, '.snplist --sort-vars \\
           --make-pgen --out ', opt$directory, opt$output, threshold, "/", data_id[f],
           '_all_scores/snp_1000G.EUR.tmp.', opt$chr ))

  system(paste0('plink2 --pfile ', opt$directory, opt$output, threshold, "/", data_id[f],
          '_all_scores/snp_1000G.EUR.tmp.', opt$chr  ' \\
           --make-bed --out ', opt$directory, opt$output, threshold, "/", data_id[f],
          '_all_scores/snp_1000G.EUR.hg38.', opt$chr ))

  #Prepare frequency files
  frqfile=(paste0(opt$directory, opt$bim, opt$chr,'.frq'))
  frq <- fread(frqfile)

  snpfile=paste0(opt$directory, opt$output, opt$threshold, "/", data_id[f],
          '_all_scores/snp_kept_', opt$chr, '.snplist')
  snp_list=scan(snpfile, what=character(0))

  frq=frq[frq$SNP %in% snp_list,]

  bim_sorted <- read.table(paste0(opt$directory, opt$output, opt$threshold, "/", data_id[f],
  '_all_scores/snp_1000G.EUR.hg38.', opt$chr, '.bim'), stringsAsFactors=F)

  write.table(frq[bim_sorted$V2,],
               file=(paste0(opt$directory, opt$output,  opt$threshold, "/", data_id[f],
               '_all_scores/snp_1000G.EUR.hg38.',opt$chr,'.frq')),
               row.names = FALSE, sep="\t", quote = FALSE)

  ### Print annotation files for baseline and complete model
  cat("Preparing SNP annotation files...\n")
  baseld <- fread(paste0(opt$directory, opt$lsbaseannot, 'baselineLD.', opt$chr,'.annot.gz'))

  baseldtoprint=baseld[baseld$SNP %in% bed_chr$SNP,]
  baseldtoprint=data.frame(baseld, stringsAsFactors=F)
  rownames(baseldtoprint)=baseldtoprint$SNP
  write.table(baseldtoprint[bim_sorted$V2,],
              file=gzfile(paste0(opt$directory, opt$output, threshold, "/", data_id[f],
                            '_all_scores/snp_baselineLD.', opt$chr,'.annot.gz')),
                  row.names = FALSE, sep="\t", quote = FALSE)


  bed_chr_tmptoprint=data.frame(bed_chr, stringsAsFactors=F)
  rownames(bed_chr_toprint)=bed_chr_toprint$SNP
  write.table(bed_chr_toprint[bim_sorted$V2,],
                  file=gzfile(paste0(opt$directory, opt$output, threshold, "/", data_id[f],
                                '_all_scores/snp_score_annot.', opt$chr,'.annot.gz')),
                  row.names = FALSE, sep="\t", quote = FALSE)
}
