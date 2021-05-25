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

d <- calcNormFactors(d)



snames <- colnames(counts) # Sample names

pdatvec = pdat$group
names(pdatvec) = pdat$samples
expgrps = pdatvec[snames]

group <- interaction(expgrps)



plotname = paste(out_dir_name, "/", "limma_voom.mds.", out_base_name, ".svg", sep="")

svglite::svglite(file = plotname, width = fig.width, height = fig.height)
plotMDS(d, col = as.numeric(group))
dev.off()

mm <- model.matrix(~0 + group)

plotname = paste(out_dir_name, "/", "limma_voom.mean_var.", out_base_name, ".svg", sep="")
svglite::svglite(file = plotname, width = fig.width, height = fig.height)
y <- voom(d, mm, plot = T)
dev.off()

fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupGrp0 - groupGrp1, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "P", n = Inf,adjust.method="BH")
head(top.table, 20)

deRes = as.data.frame(top.table)

finalres = data.frame(
    PROBEID=rownames(deRes),
    FC=deRes$logFC,
    PVAL=deRes$P.Value,
    ADJ.PVAL=deRes$adj.P.Val)
    
write.table(finalres, file=out.file, row.names=F, quote=F, sep="\t")



