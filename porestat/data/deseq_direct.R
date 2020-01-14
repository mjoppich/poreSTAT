#!/home/proj/biosoft/software/R/R-3.2.2/bin/Rscript
############################################################
#
# author: Ludwig Geistlinger
# date: 2015-07-02 13:10:08
#
# descr: differential expression analysis of RNA-seq data
#
############################################################

# call: Rscript de_rseq.R <exprs.file> <pdat.file> <fdat.file> <out.file>
#


if (length(commandArgs()) != 9) message("usage: Rscript de_rseq.R <exprs.file> <pdat.file> <fdat.file> <out.file>")
stopifnot(length(commandArgs()) == 9)

message("Loading EnrichmentBrowser")
message("R paths")
message(.libPaths())


rUserLib = Sys.getenv("R_LIBS_USER")
message("R user path")
message(rUserLib)
dir.create(rUserLib, recursive=TRUE)


message("Adding R user Path")
.libPaths( c(  rUserLib , .libPaths()) )

message("Complete lib Paths")
message(.libPaths())



requiredPackages = c('DESeq2', "EnrichmentBrowser")
for (rPackage in requiredPackages) {
    if (! require(rPackage, character.only = TRUE))
    {
        message("Not installed: ", rPackage)
        message("Trying to install to ", rUserLib)

        message("Sourcing biocLite")
        source("https://bioconductor.org/biocLite.R")
        message("Running biocLite for ", rPackage, "in", rUserLib)
        biocLite(rPackage, lib=rUserLib,    lib.loc=.libPaths())
    }

}

#message("Loading ", rPackage)
#suppressWarnings(suppressPackageStartupMessages(library(rPackage)))

exprs.file <- commandArgs()[6]
pdat.file <- commandArgs()[7]
fdat.file <- commandArgs()[8]
out.file <- commandArgs()[9]

message("Reading data ...")
eset <- readSE(exprs.file, pdat.file, fdat.file)


message("DE analysis ...")
dds <- DESeqDataSet(eset, design = ~GROUP)
dds <- DESeq(dds)
( res <- results(dds) )
outres = res[! is.na(res$padj),]


finalres = data.frame(PROBEID=rownames(outres), FC=outres$log2FoldChange, PVAL=outres$pvalue, ADJ.PVAL=outres$padj)

write.table(finalres, file=out.file, row.names=F, quote=F, sep="\t")


