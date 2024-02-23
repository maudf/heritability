### Maud Fagny
### 2015/12/21
### analyze_cluster_GO.R
### make modularity plots and gene ontology analysis
### ________________________________________________________

### Source code
require(topGO)
require(org.Hs.eg.db)
require(data.table)
require("scriptName")
require("getopt")


### Set variables
spec = matrix(c('help','h', 0, "logical", "Print help message",
                'directory','d', 1, "character", "Path to main folder",
                'tissue','t', 1, "character", "tissue to analyse",
                'annottissue','a', 1, "character", "tissue annotation file",
                'network','n', 1, "character",  "Path to Network folder",
                'genes ','g', 1, "numeric", "Minimum number of genes in enriched category",
                'pval','p', 1, "numeric", "Threshold for pvalue to dtermine significant enrichment",
                'output','o', 1, "character", "Path to the output folder"
              ), byrow=TRUE, ncol=5)

opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)}

if ( is.null(opt$directory    ) ) {opt$directory    = 'eqtl_networks/'     }
if ( is.null(opt$tissue    ) ) {opt$tissue    = 'Adipose_Subcutaneous'     }
if ( is.null(opt$annottissue    ) ) {opt$annottissue    = 'Data/Annotations/tissues_annotation.Rds'     }
if ( is.null(opt$network    ) ) {opt$network    = 'Results/Networks/'     }
if ( is.null(opt$genes    ) ) {opt$genes    = 3     }
if ( is.null(opt$pval    ) ) {opt$threshold    = 0.05  }
if ( is.null(opt$output    ) ) {opt$output    = 'Results/GOanalysis/'     }

cat("Running", current_filename(), "with options:\n",
    "directory =", opt$directory, "\n",
    'tissue =',opt$tissue, "\n",
    'annottissue =',opt$annottissue, "\n",
    'network =',opt$network, "\n",
    'genes =',opt$genes, "\n",
    'pval =',opt$pval, "\n",
    "output =",opt$output, "\n",)

###Check if annotation file exists
if(!file.exists(opt$directory, opt$annottissue))){
  stop(paste0(opt$directory, opt$annottissue,
              " not found, missing tissue annotation file. Exit !\n"))
}

###Check if network folder exists
if(!dir.exists(opt$directory, opt$network))){
  stop(paste0(opt$directory, opt$network,
              " not found, missing network folder. Exit !\n"))
}

### Load data
communities.file <- paste0(opt$directory, opt$network, opt$tissue, "eqtlnetworks.Rds")
tissue.annot.path <- paste0(opt$directory, opt$annottissue)
communities=readRDS(communities.file)
tissue.annot <- readRDS(tissue.annot.path)
go.results.file <- paste0(opt$directory, opt$output, "Rfiles/topGO_results_", tissue, ".Rds")
go.results.txt <- paste0(opt$directory, opt$output, "Tables/topGO_results_", tissue, "_" )
go.results.pdf <- paste0(opt$directory, opt$output, "Figures/topGO_results_", tissue, ".pdf")

### Create results directory
if(!dir.exists(dirname(paste(opt$directory, opt$output, "Rfiles", sep="/")))){
  dir.create(dirname(paste(opt$directory, opt$output, "Rfiles"), sep="/")),
  recursive = T, showWarnings = F)
}

if(!dir.exists(dirname(paste(opt$directory, opt$output, "Tables", sep="/")))){
  dir.create(dirname(paste(opt$directory, opt$output, "Tables"), sep="/")),
  recursive = T, showWarnings = F)
}
if(!dir.exists(dirname(paste(opt$directory, opt$output, "Figures", sep="/")))){
  dir.create(dirname(paste(opt$directory, opt$output, "Figures"), sep="/")),
  recursive = T, showWarnings = F)
}

### Compute GO term enrichment
allResfinal=NULL # Will store GO enrichment analysis results
GOtermsGenes=NULL # Will store gene ID for each enriched GO Term
pdf(go.results.pdf, width=11.5, height=8)
for (com in sort(unique(communities$blue.memb$com))){
  cat("Running GO analysis for com", com, "...\n")
  genes.names=gsub("\\.[0-9][0-9]*", "", communities$blue.memb$blue.names)
  all.genes <-   factor(as.integer(communities$blue.memb$com == com))

  names(all.genes)=genes.names
  if(sum(communities$blue.memb$com == com)>=5){ ## Run topGO only if the list of genes to be analyzed contain more than 5 genes

    allRes=list() # List of GO enrichment analysis results
    for(ont in c("BP", "MF")){
      cat("and ontology", ont, "...\n")
      topGO = new("topGOdata", description=paste0("com_",com), ontology= ont,
                  allGenes = all.genes, nodeSize = 5,
                  annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl",)
      test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
      resultElimFisher <- getSigGroups(topGO, test.stat)
      allRes[[`ont`]] <- GenTable(topGO,
                                  elim = resultElimFisher, orderBy = "elim",
                                  ranksOf = "elim", topNodes = length(resultElimFisher@score),
                                  numChar=3000)
      allRes[[`ont`]] <- allRes[[`ont`]][allRes[[`ont`]]$elim<=pval & allRes[[`ont`]]$Significant>=min.genes.signif,]
      allRes[[`ont`]]$ontology=rep(ont, nrow(allRes[[`ont`]]))
      allRes[[`ont`]]$com=rep(com, nrow(allRes[[`ont`]]))

      ## Plot topGO output
      showSigOfNodes(topGO, score(resultElimFisher), firstSigNodes = 10,
                     useInfo ='all')
      mtext(side=3, line=-1, text = paste("Community", com, "-", ont))

      ## Write gene list for each enriched GO ID
      myterms <- allRes[[`ont`]]$GO.ID
      mygenes  <- genesInTerm(topGO, myterms)
      if(length(myterms)>0){
        for (j in 1:length(myterms)) {
          myterm <- myterms[j]
          mygenesforterm <- mygenes[myterm][[1]]
          myfactor <- mygenesforterm %in%
            names(all.genes[communities$blue.memb$com == com]) # find the genes that are in the list of genes of interest
          mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
          mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
          GOtermsGenes <- rbind(GOtermsGenes,
                                data.frame("Ontology"=ont, "com"=com,
                                          "GOTerm"=myterm,"Genes"=mygenesforterm2))
        }
      }
    }
    allResfinal=rbind(allResfinal, allRes$BP, allRes$MF)

    ## write only significant results (3 genes minimum in the group + )
  }
}
dev.off()
topGO.all=list(enrichment=allResfinal, genelist=GOtermsGenes)

write.table(allResfinal[(allResfinal$ontology=="BP"),],
              file=paste0(go.results.txt, "BP_enrichment.txt"),
              quote=FALSE,row.names=FALSE, sep="\t")
write.table(GOtermsGenes[GOtermsGenes$Ontology=="BP",],
              file=paste0(go.results.txt, "BP_term_genes.txt"),
              quote=FALSE,row.names=FALSE, sep="\t")
write.table(allResfinal[(allResfinal$ontology=="MF"),],
            file=paste0(go.results.txt, "MF_enrichment.txt"),
            quote=FALSE,row.names=FALSE, sep="\t")
write.table(GOtermsGenes[GOtermsGenes$Ontology=="MF",],
            file=paste0(go.results.txt, "MF_term_genes.txt"),
            quote=FALSE,row.names=FALSE, sep="\t")
saveRDS(topGO.all, file=go.results.file)
