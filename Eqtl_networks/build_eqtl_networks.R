### 29-04-2022

### Load libraries
require("netZooR")
require("Matrix")
require("getopt")
require("scriptName")
require("data.table")
require("igraph")

###Set Variables
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
if ( is.null(opt$input    ) ) {opt$input    = 'Data/EQTL/'     }
if ( is.null(opt$output    ) ) {opt$output    = 'Results/Networks/'     }

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


###Load data
qtl=readRDS(paste0(opt$directory, opt$input))

### Create condor network object
elist <- data.frame("red"=qtl$snps[qtl$FDR<=0.2],
                    "blue"=qtl$gene[qtl$FDR<=0.2])
condor.object <- createCondorObject(elist)
network <- condorCluster(condor.object, cs.method = "LCS", project = TRUE,
                                low.memory = TRUE, fast.compute=TRUE,
                                deltaQmin = "default", # threshold for convergence of modularity (either 'default', or a number lower than 1)
                                stepred=10000, #number of red members computed together
                                stepblue=200)
### Compute core scores
qscore <- condorQscore(network)
saveRDS(qscore, file = paste0(opt$directory, opt$output, "_eqtlnetworks.Rds"))

### Extract outdegree and corescore

extract.scores.outdegree <- function(data){
  tmp=tapply(data$edges$weight, data$edges$red, sum)
  out=data.frame("SNP"=names(tmp),
                 "outdegree"=tmp, stringsAsFactors = F)
  return(out)
}

extract.scores.cs <- function(data){
  out=data.frame("SNP"=as.character(data$qscores$red.qscore$red.names),
                 "Qik"=data$qscores$red.qscore$Qik, stringsAsFactors = F)
  return(out)
}

outdegree=extract.scores.outdegree(data=qscore)
outdegree=data.table(outdegree)
setkey(outdegree, SNP)
corescores=extract.scores.cs(data=qscore)
corescores=data.table(corescores)
setkey(corescores, SNP)
scores=merge(outdegree, corescores)
saveRDS(scores, file=paste0(opt$directory, opt$output, "_scores.Rds"))
