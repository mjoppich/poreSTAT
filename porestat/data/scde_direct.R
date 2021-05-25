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



requiredPackages = c('MAST', "svglite")
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

exprs.file = "diffregs/myh11_1_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/count_noiseq"
pdat.file = "diffregs/myh11_1_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/count_p_data"
fdat.file = "diffregs/myh11_1_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/count_f_data"
out.file = "diffregs/myh11_1_cx3cr1/cx3cr1_wt_ko.umi_exon.diffreg/count_out_data_DirectLimmaVoom"

out_dir_name = dirname(out.file)
out_base_name = basename(out.file)
fig.width = 7
fig.height = 5

pdat <- read.delim(pdat.file, header=F)
colnames(pdat) = c("samples", "group")
pdat$group = paste("Grp", pdat$group, sep="")

counts <- read.delim(exprs.file, row.names = 1)
counts = counts[, pdat$samples]


#https://www.nature.com/articles/nmeth.2967
#devtools::install_version('flexmix', '2.3-13')
#BiocManager::install("scde")
require(scde)
cd <- clean.counts(counts, min.lib.size=100, min.reads = 1, min.detected = 1)
stages <- factor(pdat$group,levels=c("Grp0","Grp1"))
names(stages) <- colnames(counts)

o.ifm <- scde.error.models(counts = cd, groups = stages, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 0)

o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)

ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  stages, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# convert Z-score and corrected Z-score to 2-sided p-values 
ediff$pvalue <- 2*pnorm(-abs(ediff$Z))
ediff$p.adjust <- 2*pnorm(-abs(ediff$cZ))

# sort and look at top results
ediff <- ediff[order(abs(ediff$Z), decreasing  =  TRUE), ]
head(ediff, n=10)

finalres = data.frame(
    PROBEID=rownames(ediff),
    FC=ediff[["mle"]],
    PVAL=ediff[["pvalue"]],
    ADJ.PVAL=ediff[["p.adjust"]])

write.table(finalres, file=out.file, row.names=F, quote=F, sep="\t")

