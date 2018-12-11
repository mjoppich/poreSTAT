# Title     : MS-EmpiRE DE
# Objective : perform DE analysis using MS-EmpiRE
# Created by: mjopp
# Created on: 10/12/2018



if (length(commandArgs()) != 10) message("usage: Rscript de_rseq.R <exprs.file> <pdat.file> <fdat.file> <de.method> <out.file>")
stopifnot(length(commandArgs()) == 10)

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



requiredPackages = c('devtools', 'Biobase', 'msEmpiRe')
for (rPackage in requiredPackages) {
    if (! require(rPackage, character.only = TRUE))
    {

        if (rPackage == "msEmpiRe")
        {
            install_github("zimmerlab/MS-EmpiRe", lib=rUserLib)
        } else if (rPackage == "devtools")
        {
            message("Not installed: ", rPackage)
            message("Trying to install to ", rUserLib)

            install.packages(rPackage, lib=rUserLib)
            require(rPackage, character.only = TRUE)
        } else {
            message("Not installed: ", rPackage)
            message("Trying to install to ", rUserLib)

            message("Sourcing biocLite")
            source("https://bioconductor.org/biocLite.R")
            message("Running biocLite for ", rPackage, "in", rUserLib)
            biocLite(rPackage, lib=rUserLib,    lib.loc=.libPaths())
        }

    }

}

#message("Loading ", rPackage)
#suppressWarnings(suppressPackageStartupMessages(library(rPackage)))




exprs.file <- commandArgs()[6]
pdat.file <- commandArgs()[7]
fdat.file <- commandArgs()[8]
de.method <- commandArgs()[9]
out.file <- commandArgs()[10]

p_data <- read.csv(pdat.file, sep="\t", header=FALSE, stringsAsFactors = FALSE)

print(p_data)

p_data$V1 = sapply(p_data$V1, make.names)

print(p_data)

tmpout = tempfile()
write.table(p_data, file=tmpout, row.names=FALSE, quote=FALSE, sep="\t")

print(tmpout)


print("Reading files")
#loading data from installed data sets
data <- msEmpiRe::read.standard(fdat.file, tmpout, prot.id.generator=function(pep) pep)


print("Extracting Conditions1")
#extract the first two conditions
conditions <- msEmpiRe::extract_conditions(data)
print("Extracting Conditions2")
print(conditions)
print("Extracting Conditions3")
conditions <- conditions[, c(1,2)]


print(conditions)
print(head(data))

#removing peptides that are detected in less than 2 samples per condition
data <- msEmpiRe::filter_detection_rate(data, condition=conditions)

#normalize
data <- msEmpiRe::normalize(data)

#analyse
result <- de.ana(data)

print(head(result))

myres <- data.frame(
   GENE.ID = result$prot.id,
   log2FC = result$log2FC,
   RAW.PVAL = result$p.val,
    ADJ.PVAL = result$p.adj
)

print(head(myres))

write.table(myres, file=out.file, row.names=F, quote=F, sep="\t")

