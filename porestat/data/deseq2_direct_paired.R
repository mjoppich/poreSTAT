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


if (length(commandArgs()) != 8) message("usage: Rscript de_rseq.R <exprs.file> <pdat.paired.file> <out.file>")
stopifnot(length(commandArgs()) == 8)

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



requiredPackages = c('DESeq2')
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
out.file <- commandArgs()[8]

message("Reading data ...")
cts <- read.csv(exprs.file,sep="\t", header=TRUE)
rownames(cts) = cts$gene
cts$gene = NULL

cts = as.matrix(cts)
print(head(cts))

sampleTable = read.csv(pdat.file, sep="\t", header=TRUE)
rownames(sampleTable) = sampleTable$name
sampleTable$name = NULL

sampleTable$condition <- factor(sampleTable$condition)
sampleTable$sample <- factor(sampleTable$sample)

print(sampleTable)

cts = cts[,rownames(sampleTable)]

message("DE analysis ...")

print("Input matrix")
maxzero = 0.75 * nrow(sampleTable)
print(dim(cts))
cts = cts[!rowSums(cts == 0) >= maxzero,]
print(dim(cts))


dds <- DESeqDataSetFromMatrix(countData =cts, colData = sampleTable, design = ~sample+condition)
dds <- DESeq(dds)
( res <- results(dds) )
outres = res[! is.na(res$padj),]
outres <- outres[order(outres$padj),]


finalres = data.frame(PROBEID=rownames(outres), FC=outres$log2FoldChange, PVAL=outres$pvalue, ADJ.PVAL=outres$padj)

print(dim(finalres[finalres$ADJ.PVAL < 0.05, ]))

write.table(finalres, file=out.file, row.names=F, quote=F, sep="\t")


