### plot_GO_module.R
### Load libraries
require("scriptName")
require("getopt")
require("ggplot2")
require("gplots")
require("plyr")
require("RColorBrewer")
require("dendextend")
require("wordcloud")
require("wesanderson")

### Set variables
spec = matrix(c('help','h', 0, "logical", "Print help message",
                'directory','d', 1, "character", "Path to main folder",
                'tissue','t', 1, "character", "tissue to analyse",
                'annottissue','a', 1, "character", "tissue annotation file",
                'module','m', 1, "integer",  "ID of the module to analyse",
                'pval','p', 1, "numeric", "Threshold for pvalue to dtermine significant enrichment",
                'input','i', 1, "character", "Path to the input folder"
              ), byrow=TRUE, ncol=5)

opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)}

if ( is.null(opt$directory    ) ) {opt$directory    = 'eqtl_networks/'     }
if ( is.null(opt$tissue    ) ) {opt$tissue    = 'Adipose_Subcutaneous'     }
if ( is.null(opt$annottissue    ) ) {opt$annottissue    = 'Data/Annotations/tissues_annotation.Rds'     }
if ( is.null(opt$module    ) ) {opt$module    = 142    }
if ( is.null(opt$pval    ) ) {opt$threshold    = 0.01 }
if ( is.null(opt$input    ) ) {opt$input    = 'Results/GOanalysis/'     }

cat("Running", current_filename(), "with options:\n",
    "directory =", opt$directory, "\n",
    'tissue =',opt$tissue, "\n",
    'annottissue =',opt$annottissue, "\n",
    'module =',opt$module, "\n",
    'pval =',opt$pval, "\n",
    "input =",opt$input, "\n",)

###Check if annotation file exists
if(!file.exists(opt$directory, opt$annottissue))){
  stop(paste0(opt$directory, opt$annottissue,
              " not found, missing tissue annotation file. Exit !\n"))
}

###Check if input folder exists
if(!dir.exists(paste0(opt$directory, opt$input))){
  stop(paste0(opt$directory, opt$input,
              " not found, missing input folder. Exit !\n"))
}

###Check if output folder exists
if(!dir.exists(paste(opt$directory, opt$input, "Figures", sep="/"))){
  dir.create(paste(opt$directory, opt$input, "Figures", sep="/"), recursive = T, showWarnings = F)
}


### Load data
tissue.annot.path <- paste0(opt$directory,  opt$annottissue)
tissue.annot <- readRDS(tissue.annot.path)

### Load data
go.results.rds <- paste0(opt$directory, opt$input, "Rfiles/topGO_results_", opt$tissue, ".Rds")
go.results.bp <- paste0(opt$directory, opt$input, "Tables/topGO_results_", opt$tissue, "_BP_enrichment.txt")
go <- read.table(go.results.bp, header=T, stringsAsFactors = F, quote="", sep="\t")

### Plot wordcloud plot
wordcloud.pdf <- paste0(opt$directory, opt$input, "Figures/Gene_Ontology_wordcloud_", opt$tissue, "_", opt$module, ".pdf")

go.com=go[go$com==module & go$elim<=pval,]
pdf(wordcloud.pdf, h=6, w=6)
words = unlist(strsplit(go.com$Term, " "))
words=words[!(words %in% c("of", "the", "a", "in"))]
freq=table(words)/length(words)
wordcloud(names(freq), freq, scale=c(3,.5),
          colors=wes_palette("Zissou1", 10, type = "continuous"))
dev.off()

### Plot bubble plot

bubble.plot.unique.top10 <- paste0(opt$directory, opt$input, "Figures/Gene_Ontology_top10_", opt$tissue, "_", opt$module, ".pdf")

if(nrow(go.com)>10){go.com=go.com[1:10,]}
logp <- -log10(as.numeric(as.character(go.com$elim)))
Term <- factor(go.com$Term, rev(as.character(go.com$Term)))
main <- paste(tissue.annot[tissue], "- Module", module)
OddsRatio <- as.numeric(as.character(go.com$Significant))/as.numeric(as.character(go.com$Expected))
OddsRatio[OddsRatio==Inf]<- max(OddsRatio[OddsRatio!=Inf])
res=data.frame(logp, Term, OddsRatio)
pdf(bubble.plot.unique.top10, h=4, w=8)
g<-ggplot(res, aes(x=logp, y=Term, size=OddsRatio))
print(g+geom_point( shape=21, colour="black", fill="dodgerblue3")+
        guides(fill=guide_legend(title.theme = element_text(size = 15)))+
        scale_size_continuous(range = c(0,
                                        30)*0.5,
                              breaks = seq( 0, 30, 5))+
        labs(
          x = "-log10(p-value)",
          y = "",
          size = "Significant/\nExpected"
        )+
        theme_bw()+
        xlim(0, max(res$logp)*1.3)+
        #ylim(0, length(res$logp)*1.1)+
        theme(axis.text.y = element_text(colour="black", size = rel(1.5)))+
        theme(axis.text.x = element_text(colour="black", size = rel(1.5)))+
        theme(axis.title = element_text(colour="black", size = rel(1.5)))+
        theme(legend.text = element_text(colour="black", size = rel(1.5)))+
        theme(plot.margin=unit(c(1, 0, 2, 0), unit="pt"),
              legend.key.width=unit(2, unit="pt"),
              legend.key.height=unit(8, unit="pt"),
              legend.key = element_blank()) # test
)
dev.off()
