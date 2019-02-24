#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#----------------------------------R API -------------------------------------------------------------------
library(recount);library(GenomicRanges);library(limma);library(edgeR);library(DESeq2);library(regionReport);library(clusterProfiler)
library(org.Hs.eg.db);library(gplots);library(derfinder);library(rtracklayer);library(GenomicFeatures);library(bumphunter)
library(derfinderPlot);library(devtools);library(xtable);library(Peptides)

#----------------------Data and Programs for Classroom Lecture-------------------------------------
source("https://bioconductor.org/biocLite.R") # additional functions
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
source("https://bioconductor.org/biocLite.R") # additional functions
biocLite("BSgenome.Hsapiens.UCSC.hg19")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tx.Hsapiens.seqs <- extractTranscriptSeqs(Hsapiens, TxDb.Hsapiens.UCSC.hg19.knownGene, use.names=TRUE)

getChromInfoFromUCSC("hg19")
seqlevels(txdb)<-"chr8"
GR.8 <- transcripts(txdb) 
length(GR.8)
GR.8.Gene.List <- transcriptsBy(txdb, by = "gene")
GR.8.Genes<-genes(txdb, columns="gene_id", filter=NULL)
names(GR.8.Gene.List)[1:5]
tx.Hsapiens.8.seqs <- extractTranscriptSeqs(Hsapiens,txdb, use.names=TRUE)
tx.Hsapiens.8.seqs.translations<-suppressWarnings(translate(tx.Hsapiens.8.seqs))

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
seqlevels(txdb)<-"chr18"
GR.18 <- transcripts(txdb) 
length(GR.18)
GR.18.Gene.List <- transcriptsBy(txdb, by = "gene")
names(GR.18.Gene.List)[1:5]
tx.Hsapiens.18.seqs <- extractTranscriptSeqs(Hsapiens,txdb, use.names=TRUE)
tx.Hsapiens.18.seqs.translations<-suppressWarnings(translate(tx.Hsapiens.18.seqs))

bind.potential.8.df<-Protein.Protein.Interaction.Index(tx.Hsapiens.8.seqs.translations,2.48)
bind.potential.8.pct<-sum(bind.potential.8.df[,1])/length(tx.Hsapiens.8.seqs.translations)

#---------------------------------------Tables-----------------------------------
dir.create('renderReport-example', showWarnings = FALSE, recursive = TRUE)
#---------------------------------------Figures----------------------------------
Figure.1<-hist(bind.potential.8.df[,1])
#--------------------------------------References--------------------------------

citation('regionReport')

#--------------------------------------Function Library--------------------------
Protein.Protein.Interaction.Index<-function(x,threshold)
{
  bind.potential<-NULL
  protein.names<-NULL
  bind.potential.df<-data.frame()
  for(i in 1:length(x))
  {
    if(boman(x[[i]]) > threshold)
      bind.potential[i]<-1
    else 
      bind.potential[i]<-0
  }
  bind.potential.df<-cbind(bind.potential)
  colnames(bind.potential.df)<-c("Bind Potential")
  return(bind.potential.df)
}

