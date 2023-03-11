############################################################
#
# author: Markus Joppich
# date: 2020
#
# descr: differential expression analysis of RNA-seq data using DESeq2
#
############################################################

# call: Rscript de_rseq.R <exprs.file> <pdat.file> <fdat.file> <out.file>
#


if (length(commandArgs()) != 9) message("usage: Rscript de_rseq.R <exprs.file> <pdat.file> <fdat.file> <out.file>")
stopifnot(length(commandArgs()) == 9)

message("Loading Libraries")
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



requiredPackages = c('DESeq2', "svglite", "summarytools")
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

#exprs.file = "diffregs/myh11_2_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/count_df"
#exprs.file = "diffregs/kupffer_cells/kupffer_cells.umiexon.diffreg/count_df"
#pdat.file = "diffregs/myh11_2_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/count_p_data"
#pdat.file = "diffregs/kupffer_cells/kupffer_cells.umiexon.diffreg/count_p_data"
#fdat.file = "diffregs/myh11_1_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/count_f_data"
#out.file = "diffregs/myh11_2_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/count_out_data_delboy"
#out.file = "diffregs/kupffer_cells/kupffer_cells.umiexon.diffreg/count_out_data_delboy"


out_dir_name = dirname(out.file)
out_base_name = basename(out.file)
fig.width = 7
fig.height = 5

message("Reading data ...")

pdat <- read.delim(pdat.file, header=F)
colnames(pdat) = c("id", "group")
pdat$group = paste("Grp", pdat$group, sep="")

counts <- read.delim(exprs.file, row.names = 1)
counts = counts[, pdat$id]

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=pdat, 
                              design=~group, tidy = FALSE)


message("DE analysis ...")

dds <- estimateSizeFactors(dds)

nQuant = 7
calcQuantiles = c(0.25, 0.5, 0.75, 0.8, 0.9, 0.925, 0.95)
ncountQuantiles = quantile(c(as.matrix(counts(dds, normalized=TRUE))), probs=calcQuantiles)

print(paste("Filtering normalized counts by quantile=", calcQuantiles[nQuant], "and value=", ncountQuantiles[nQuant]))

keep = rowSums( counts(dds, normalized=TRUE) >= ncountQuantiles[nQuant] ) >= 3



dds <- dds[keep,]


dds <- DESeq(dds)

print(head(dds$baseMean))


res <- results(dds, alpha=0.05)

print(summary(res))

res <- res[order(res$pvalue),]
res$padj[is.na(res$padj)] = 1

resNoFilt <- results(dds, independentFiltering=FALSE, alpha=0.05, cooksCutoff=FALSE)
resNoFilt <- resNoFilt[order(resNoFilt$pvalue),]

write.table(resNoFilt, file="test.tsv", row.names=T, quote=F, sep="\t")

plotname = paste(out_dir_name, "/", "DirectDESeq2.FilterThreshold.", out_base_name, ".svg", sep="")
svglite::svglite(file = plotname, width = fig.width, height = fig.height)
plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
dev.off()

addmargins(table(filtering=(res$padj < .05),
                 noFiltering=(resNoFilt$padj < .05)))


plotname = paste(out_dir_name, "/", "DirectDESeq2.topgenes.", out_base_name, ".svg", sep="")
svglite::svglite(file = plotname, width = fig.width, height = fig.height)
par(mfrow=c(2,3))

for (geneName in rownames(head(res, n=6)))
{
    plotCounts(dds, gene=geneName, intgroup="group")
}
dev.off()

plotname = paste(out_dir_name, "/", "DirectDESeq2.ma.", out_base_name, ".svg", sep="")
svglite::svglite(file = plotname, width = fig.width, height = fig.height)
plotMA(res, ylim=c(-2,2))
dev.off()

resLFC <- lfcShrink(dds, coef="group_Grp1_vs_Grp0", type="apeglm")
resNorm <- lfcShrink(dds, coef="group_Grp1_vs_Grp0", type="normal")
resAsh <- lfcShrink(dds, coef="group_Grp1_vs_Grp0", type="ashr")

plotname = paste(out_dir_name, "/", "DirectDESeq2.ma_by_methods.", out_base_name, ".svg", sep="")
svglite::svglite(file = plotname, height = 6, width = 15)
pltXLim = 10 ** ceiling(log10(max(counts(dds, normalized=T))))
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,pltXLim); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
dev.off()


plotname = paste(out_dir_name, "/", "DirectDESeq2.disp_estimates.", out_base_name, ".svg", sep="")
svglite::svglite(file = plotname, width = fig.width, height = fig.height)
plotDispEsts(dds)
dev.off()


normExpr = counts(dds, normalized=T)
cn = colnames(normExpr)
normExpr2 = cbind(rownames(normExpr), normExpr)
colnames(normExpr2) = c("Geneid", cn)


outres = res[! is.na(res$padj),]
finalres = data.frame(PROBEID=rownames(outres), FC=outres$log2FoldChange, PVAL=outres$pvalue, ADJ.PVAL=outres$padj)

write.table(finalres, file=out.file, row.names=F, quote=F, sep="\t")

out_dir_name = dirname(out.file)
out_base_name = basename(out.file)
write.table(normExpr2, file=paste(out_dir_name, "/", "norm_expr.", out_base_name, sep=""), row.names=F, quote=F, sep="\t")
