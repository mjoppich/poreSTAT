#!/home/proj/biosoft/software/R/R-3.2.2/bin/Rscript
############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-07-02 13:10:08
# 
# descr: differential expression analysis of RNA-seq data
# 
############################################################

# call: Rscript de_rseq.R <exprs.file> <pdat.file> <fdat.file> <de.method> <out.file>
#
# choose <de.method> out of {'limma', 'edgeR', 'DESeq'}

if(length(commandArgs()) != 10) message("usage: Rscript de_rseq.R <exprs.file> <pdat.file> <fdat.file> <de.method> <out.file>")
stopifnot(length(commandArgs()) == 10)

message("Loading EnrichmentBrowser")
print(.libPaths())
suppressWarnings(suppressPackageStartupMessages(library(EnrichmentBrowser)))

exprs.file <- commandArgs()[6]
pdat.file <- commandArgs()[7]
fdat.file <- commandArgs()[8]
de.method <- commandArgs()[9]
out.file <- commandArgs()[10]

message("Reading data ...")
eset <- read.eset(exprs.file, pdat.file, fdat.file)

message("Removing genes with low read count ...")
eset <- eset[rowSums(exprs(eset), na.rm=TRUE) > ncol(eset),]

message("DE analysis ...")
eset <- de.ana(eset, de.method=de.method, padj.method="none")

de.tbl <- fData(eset)[,sapply(c("FC.COL","ADJP.COL"), config.ebrowser)]
de.tbl <- cbind(de.tbl, p.adjust(de.tbl[,2], method="BH"))
de.tbl <- cbind(featureNames(eset), de.tbl)
colnames(de.tbl) <- c("GENE.ID", "log2FC", "RAW.PVAL", "ADJ.PVAL")

write.table(de.tbl, file=out.file, row.names=FALSE, quote=FALSE, sep="\t")
message("DE table written to ", out.file)
