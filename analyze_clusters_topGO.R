### Maud Fagny
### 2015/12/21
### analyze_cluster_GO.R
### make modularity plots and gene ontology analysis
### ________________________________________________________

### Source code
require(topGO)
#require(GOstats, lib.loc = "/shared/ifbstor1/home/mfagny/R/x86_64-conda-linux-gnu-library/4.2/")
require(org.Hs.eg.db)
require(data.table)

### Set variables
args=commandArgs(trailingOnly=TRUE)
tissue <- args[1] #"Adipose_Subcutaneous"
rootdir <- args[2] #"Documents/Heritability_eqtl_networks/"
networkrootfilename <- args[3] #"Results/Networks/Weighted/network_qscores_"
resultsfilesdir <- args[4] #"Results/GOanalysis/"
tissue.annot.file <- args[5] #"Data/Annotations/tissues_annotation.rds"
min.genes.signif <- args[6] #3
pval <- args[7] # 0.05

### Load data
communities.file <- paste0(rootdir, networkrootfilename, tissue, ".Rds")
tissue.annot.path <- paste0(rootdir, tissue.annot.file)
communities=readRDS(communities.file)
tissue.annot <- readRDS(tissue.annot.path)
go.results.file <- paste0(rootdir, resultsfilesdir, "Rfiles/topGO_results_", tissue, ".Rds")
go.results.txt <- paste0(rootdir, resultsfilesdir, "Tables/topGO_results_", tissue, "_" )
go.results.pdf <- paste0(rootdir, resultsfilesdir, "Figures/topGO_results_", tissue, ".pdf")

### Create results directory
dir.create(paste0(rootdir, resultsfilesdir, "Rfiles/"), recursive = T, showWarnings = F) 
dir.create(paste0(rootdir, resultsfilesdir, "Tables/"), recursive = T, showWarnings = F) 
dir.create(paste0(rootdir, resultsfilesdir, "Figures/"), recursive = T, showWarnings = F) 

allResfinal=NULL
GOtermsGenes=NULL
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

