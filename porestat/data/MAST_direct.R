############################################################
#
# author: Markus Joppich
# date: 2021
#
# descr: differential expression analysis of RNA-seq data using MAST
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

#tutorial
#https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2019-02-04/labs/Differential_gene_expression.html
#BiocManager::install("MAST")
#https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5

fdata <- data.frame(rownames(counts))
sca <- FromMatrix(log2(as.matrix(counts)+1),pdat,fdata)

cdr2 <-colSums(assay(sca)>0)
colData(sca)$cngeneson <- scale(cdr2)

cond<-factor(colData(sca)$group)
colData(sca)$condition<-cond
zlmCond <- zlm(~condition + cngeneson, sca)
summaryCond <- summary(zlmCond, doLRT='conditionGrp1') 

summaryDt <- summaryCond$datatable


fcHurdle <- merge(summaryDt[contrast=='conditionGrp1' & component=='H',.(primerid, `Pr(>Chisq)`)], 
     summaryDt[contrast=='conditionGrp1' & component=='logFC', 
               .(primerid, coef, ci.hi, ci.lo)], by='primerid') 
# change name of Pr(>Chisq) to fdr and sort by fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle <- fcHurdle[order(fdr),]

finalres = data.frame(
    PROBEID=fcHurdle[["primerid"]],
    FC=fcHurdle[["coef"]],
    PVAL=fcHurdle[["Pr(>Chisq)"]],
    ADJ.PVAL=fcHurdle[["fdr"]])
    
finalres <- finalres[order(finalres$PVAL),]

write.table(finalres, file=out.file, row.names=F, quote=F, sep="\t")


quit(status=0, save='no')


#https://academic.oup.com/bioinformatics/article/34/18/3223/4983067
BiocManager::install("DEsingle")

library(DEsingle)
results <- DEsingle(counts = counts, group = factor(pdat$group),parallel=TRUE)

finalres = data.frame(
    PROBEID=rownames(results),
    FC=results[["norm_foldChange"]],
    PVAL=results[["pvalue"]],
    ADJ.PVAL=results[["pvalue.adj.FDR"]])



#https://academic.oup.com/bioinformatics/article/35/24/5155/5514046
devtools::install_github("cz-ye/DECENT")
library("DECENT")

# DECENT without spike-ins
de.table <- decent(data.obs = counts,
                   X = ~as.factor(pdat$group), 
                   use.spikes = F,
                   CE.range = c(0.02, 0.1) # specify the range of the ranked random capture efficiency
                   )


#https://www.sciencedirect.com/science/article/pii/S0888754321000744
devtools::install_github("sam-uofl/SwarnSeq")
library("SwarnSeq")

fGroup = as.numeric(as.factor(pdat$group))
results2 <- SwarnSeq(CountData=counts, RNAspike.use=FALSE, CE.range=c(0.1, 0.4), parallel=T, norm.method="TMM", group=fGroup, CellCluster=fGroup, maxit=500, eps=1e-10,CellAuxil=NULL,muoffset=NULL, phioffset=NULL, weights=NULL)

finalres = data.frame(
    PROBEID=rownames(results2),
    FC=results2[["mle"]],
    PVAL=results2[["pvalue"]],
    ADJ.PVAL=results2[["p.adjust"]])

write.table(finalres, file=out.file, row.names=F, quote=F, sep="\t")

