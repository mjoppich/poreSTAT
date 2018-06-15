if (length(commandArgs()) != 7) message("usage: Rscript noiseq_diffreg.R <exprs.file> <out.file>")
stopifnot(length(commandArgs()) == 7)

message("Loading NOISeq")
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



requiredPackages = c('NOISeq', 'data.table')
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


infile <- commandArgs()[6]
outfile <- commandArgs()[7]


df = read.csv(infile, header=TRUE, sep="\t", check.names=FALSE)
df2 <- df[, -1]
rownames(df2) = df[, 1]


sampleNames = print(colnames(df2))
print(sampleNames)


myfactors = data.frame( Samples=colnames(df2))
mydata = readData(df2, factors=myfactors)

myresults <- noiseq(mydata, factor="Samples",k=NULL,norm="uqua",pnr=0.5,nss=5, v = 0.02,lc=0, replicates="no")

der = degenes(myresults, q=0.2, M=NULL)

dfr = setDT(der, keep.rownames = TRUE)[]

#> colnames(der)
#[1] "rn"
#[2] "X170329_2d_sequencing_run_1.1_p12_pooled_mean"
#[3] "X170330_2d_sequencing_run_chip1_run2_khg.h70_pooled_mean"
#[4] "M"
#[5] "D"
#[6] "prob"
#[7] "ranking"



colnames(dfr)[1] <- "GENE.ID"
colnames(dfr)[4] <- "log2FC"


dfr$RAW.PVAL = 2*pnorm(-abs(scale(dfr$prob)))
adjPVals = p.adjust(dfr$RAW.PVAL, method='BH', n=length(dfr$RAW.PVAL))
dfr$ADJ.PVAL = adjPVals


for (i in 1:length(sampleNames))
{
    colnames(dfr)[i+1] = sampleNames[i]
}

print(colnames(dfr))


write.table(dfr, file=outfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=colnames(dfr))
