if (length(commandArgs()) != 7) message("usage: Rscript noiseq_diffreg.R <exprs.file> <out.file> ... conditions ...")
stopifnot(length(commandArgs()) >= 7)

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

allArgs = commandArgs()
print(allArgs)

infile <- commandArgs()[6]
outfile <- commandArgs()[7]

conditionData = allArgs[8:length(allArgs)]
print(conditionData)

df = read.csv(infile, header=TRUE, sep="\t", check.names=FALSE)
df2 <- df[, -1]
rownames(df2) = df[, 1]


sampleNames = colnames(df2)
print(sampleNames)


myfactors = data.frame( Samples=colnames(df2), Conditions=conditionData )

print(myfactors)

mydata = readData(df2, factors=myfactors)

print(mydata)

if (length(conditionData) == 2)
{
    print("no replicates")
    myresults <- noiseq(mydata, factor="Samples",k=NULL,norm="uqua",pnr=0.2, nss=5, v = 0.02, replicates="no")
} else {
    print("with replicates")
    myresults <- noiseq(mydata, factor="Conditions",k=NULL,norm="uqua",pnr=0.2, replicates="technical", conditions=conditionData)
}



der = degenes(myresults, q=0.0, M=NULL)
dfr = setDT(der, keep.rownames = TRUE)[]

#> colnames(der)
#[1] "rn"
#[2] "s1_pooled_mean"
#[3] "s2_pooled_mean"
#[4] "M"
#[5] "D"
#[6] "prob"
#[7] "ranking"



colnames(dfr)[1] <- "GENE.ID"
colnames(dfr)[4] <- "log2FC"



dfr$RAW.PVAL = 2*pnorm(-abs(scale(dfr$ranking)))
adjPVals = p.adjust(dfr$RAW.PVAL, method='BH', n=length(dfr$RAW.PVAL))
dfr$ADJ.PVAL = adjPVals

uCond = unique(conditionData)

for (i in 1:length(uCond))
{
    colnames(dfr)[i+1] = uCond[i]
}

print(colnames(dfr))

write.table(dfr, file=outfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=colnames(dfr))

DE.plot(myresults, q = 0.8, graphic = "MD", log.scale = TRUE)
