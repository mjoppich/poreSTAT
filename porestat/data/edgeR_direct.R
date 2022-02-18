############################################################
#
# author: Markus Joppich
# date: 2021
#
# descr: differential expression analysis of RNA-seq data using limma-voom
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



requiredPackages = c('edgeR', "svglite")
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

#exprs.file = "diffregs/myh11_1_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/count_noiseq"
#pdat.file = "diffregs/myh11_1_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/count_p_data"
#fdat.file = "diffregs/myh11_1_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/count_f_data"
#out.file = "diffregs/myh11_1_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/count_out_data_DirectLimmaVoom"

out_dir_name = dirname(out.file)
out_base_name = basename(out.file)
fig.width = 7
fig.height = 5

#pdat <- read.delim("diffregs/myh11_1_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/pdat2", header=F)
pdat <- read.delim(pdat.file, header=F)
colnames(pdat) = c("samples", "group")
pdat$group = paste("Grp", pdat$group, sep="")

counts <- read.delim(exprs.file, row.names = 1)
counts = counts[, pdat$samples]
d0 <- DGEList(counts,group=pdat$group)

cutoff <- 2
drop <- apply(cpm(d0), 1, max) >= cutoff
d <- d0[drop,] 
dim(d) # number of genes left

d <- calcNormFactors(d, method="TMM")




plotname = paste(out_dir_name, "/", "DirectEdgeR.mds.", out_base_name, ".svg", sep="")

svglite::svglite(file = plotname, width = fig.width, height = fig.height)
plotMDS(d, col = as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:length(unique(d$samples$group)), pch=20)
dev.off()

d1 <- estimateCommonDisp(d, verbose=T)
d1 <- estimateTagwiseDisp(d1)

#plotname = paste(out_dir_name, "/", "DirectEdgeR_bcv.", out_base_name, ".pdf", sep="")
#pdf(plotname, width=7, height=5)
#plotBCV(d1)
#dev.off()

#et12 <- exactTest(d1)
#de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
#de1tags12 <- rownames(d1)[as.logical(de1)] 
#plotname = paste(out_dir_name, "/", "DirectEdgeR_smear.", out_base_name, ".pdf", sep="")
#pdf(plotname, width=7, height=5)
#plotSmear(et12, de.tags=de1tags12)
#abline(h = c(-2, 2), col = "blue")
#dev.off()


design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="auto")
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotname = paste(out_dir_name, "/", "DirectEdgeR.bcv2.", out_base_name, ".svg", sep="")
svglite::svglite(file = plotname, width = fig.width, height = fig.height)
plotBCV(d2)
dev.off()

fit <- glmFit(d2, design.mat)
lrt12 <- glmLRT(fit, contrast=c(1,-1))

de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.05)
de2tags12 <- rownames(d2)[as.logical(de2)]
plotname = paste(out_dir_name, "/", "DirectEdgeR.smear2.", out_base_name, ".svg", sep="")
svglite::svglite(file = plotname, width = fig.width, height = fig.height)
plotSmear(lrt12, de.tags=de2tags12)
abline(h = c(-2, 2), col = "blue")
dev.off()

deRes = as.data.frame(topTags(lrt12, n=Inf, adjust.method = "BH"))

finalres = data.frame(
    PROBEID=rownames(deRes),
    FC=-deRes$logFC,
    PVAL=deRes$PValue,
    ADJ.PVAL=deRes$FDR)

finalres <- finalres[order(finalres$PVAL),]

    
write.table(finalres, file=out.file, row.names=F, quote=F, sep="\t")


