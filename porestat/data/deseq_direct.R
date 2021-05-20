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

#exprs.file = "diffregs/neutrophils_excl.umi.all.diffreg/count_expr"
#pdat.file = "diffregs/neutrophils_excl.umi.all.diffreg/count_p_data"
#fdat.file = "diffregs/neutrophils_excl.umi.all.diffreg/count_f_data"
#out.file = "diffregs/neutrophils_excl.umi.all.diffreg/count_out_data_DirectDESeq2"

message("Reading data ...")
eset <- readSE(exprs.file, pdat.file, fdat.file)


message("DE analysis ...")
dds <- DESeqDataSet(eset, design = ~GROUP)
dds <- DESeq(dds)
( res <- results(dds) )
outres = res[! is.na(res$padj),]

normExpr = counts(dds, normalized=T)
cn = colnames(normExpr)
normExpr2 = cbind(rownames(normExpr), normExpr)
colnames(normExpr2) = c("Geneid", cn)

finalres = data.frame(PROBEID=rownames(outres), FC=outres$log2FoldChange, PVAL=outres$pvalue, ADJ.PVAL=outres$padj)

write.table(finalres, file=out.file, row.names=F, quote=F, sep="\t")

out_dir_name = dirname(out.file)
out_base_name = basename(out.file)
write.table(normExpr2, file=paste(out_dir_name, "/", "norm_expr.", out_base_name, sep=""), row.names=F, quote=F, sep="\t")


fresSorted = finalres[order(finalres$ADJ.PVAL), ]
library(ggplot2)

geneList = fresSorted$PROBEID[1:min(10, length(fresSorted$PROBEID))]

if ("ENSG00000011422" %in% finalres$PROBEID)
{
    geneList = c("ENSG00000011422", geneList)
}

for (gene in geneList)
{
    d <- plotCounts(dds, gene=gene, intgroup="GROUP", returnData=T)
    p = ggplot(d, aes(x=GROUP, y=count)) + geom_point(position=position_jitter(w=0.1,h=0))

    statdf = fresSorted[fresSorted$PROBEID == gene,]
    statView = toString(statdf)

    p = p + labs(title=gene, subtitle=statView)
    ggsave(paste(out.file, ".", gene, ".png", sep=""), plot=p)
}