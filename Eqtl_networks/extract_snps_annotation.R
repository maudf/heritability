### 08-08-2022
### Maud Fagny
### extract_snps_annot.R
### extract SNP annotation from a data file

### Load libraries
require("getopt")
require("scriptName")

### Load arguments
spec = matrix(c('help','h', 0, "logical",
                'directory','d', 1, "character",
                'input','i', 1, "character",
                'output','o', 1, "character"
		              ), byrow=TRUE, ncol=4)

opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)}

if ( is.null(opt$directory    ) ) {opt$directory    = '~/Documents/Heritability_eqtl_networks/'     }
if ( is.null(opt$input    ) ) {opt$input    = 'Results/Networks/'     }
if ( is.null(opt$output    ) ) {opt$output    = 'Data/Annotations/SNPs_network'     }

cat("Running", current_filename(), "with options:\n",
    "directory =", opt$directory, "\n",
    "input =", opt$input, "\n",
    "output =", opt$output, "\n")

###Check if input file exists
if(!file.exists(opt$directory, opt$input))){
  stop(paste0(opt$directory, opt$input,
              " not found, missing input file. Exit !\n"))
}

### Create output folder if do not exists
if(!dir.exists(dirname(paste(opt$directory, opt$output, sep="/")))){
  dir.create(dirname(paste(opt$directory, opt$output, sep="/")),
  recursive = T, showWarnings = F)
}


### Set variables
setwd(base.dir)
files=list.files(paste0(opt$directory, opt$input), pattern = "*_eqtlnetworks.Rds")

### Define function to get snp list from score file. Change to fit your need
get.snp.list <- function(filename){
  data <- readRDS(filename)
  out=unique(data$red.memb$red.names)
  return(out)
}


### Extract list of SNPs
list.snps=NULL
for (f in files){
  print(f)
  snps=get.snp.list(paste0(opt$directory, opt$input, f))
  list.snps=unique(sort(c(list.snps, snps)))
}

### Split SNP information
all.snps=t(apply(data.frame(list.snps), 1, function(x){strsplit(x, "_")[[1]]}))
colnames(all.snps)=c("CHR", "BP", "ref", "alt", "build")
snps.annot=data.frame("SNP"=list.snps, all.snps[, c("CHR", "BP", "ref", "alt")])
snps.annot=snps.annot[order(snps.annot$CHR, snps.annot$BP), ]

### Save SNP annotations in file
saveRDS(snps.annot, file=paste0(opt$directory, opt$output, ".Rds"))
write.table(snps.annot, file=paste0(opt$directory, opt$output, ".mapping"),
            quote=F, row.names=F, sep="\t")
