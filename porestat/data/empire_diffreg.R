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


filter_detection_rate <- function(ffdata, condition, rate=2)
{
  x <- exprs(ffdata)

  s <- (x > 0) %*% condition

    checkDF = subset(s, !(colnames(condition)[1] == 0 & colnames(condition)[2] == 0))
    checkDF = subset(checkDF, (colnames(condition)[1] >= rate & colnames(condition)[2] >= rate) | (colnames(condition)[1] == 0 | colnames(condition)[2] == 0))
    f <- apply(checkDF, 1, any)

  print(sprintf("Removed %d peptides because of low detection rate, i.e. detected in less than %d samples)", nrow(x) - sum(f), rate))

  samples_to_keep <- rowSums(condition) > 0


  exp <- exprs(ffdata)[f, samples_to_keep]
  f_dat <- data.frame(fData(ffdata)[f, ])
  colnames(f_dat) <- colnames(fData(ffdata))
  p_dat <- data.frame(condition=pData(ffdata)$condition[samples_to_keep], row.names=rownames(pData(ffdata))[samples_to_keep])
  p_dat <- Biobase::AnnotatedDataFrame(p_dat)
  r <- Biobase::ExpressionSet(exp, p_dat)
  fData(r) <- f_dat
  return(r)
}

filderDetectionRate <- function(data, rate=2, condition=NULL)
{
  x <- exprs(data)

  if(is.null(condition))
  {
    condition <- extract_conditions(data)
  }

  s <- (x > 0) %*% condition

    s = data.frame(s)

    colnames(s) = c("V1", "V2")

    checkDF = s
    f = ((checkDF$V1!=0 | checkDF$V2!=0) & (checkDF$V1 >= rate & checkDF$V2 >= rate) | (checkDF$V1 == 0 & checkDF$V2 >= rate) | (checkDF$V2 == 0 & checkDF$V1 >= rate))

    print("f")
    print(head(f))

    #f <- apply(s >= 0, 1, all)

    write.table(checkDF[f,], file="/tmp/markusout", row.names=T, quote=F, sep="\t")

  print(sprintf("Removed %d peptides because of low detection rate, i.e. detected in less than %d samples)", nrow(x) - sum(f), rate))

  samples_to_keep <- rowSums(condition) > 0

  exp <- exprs(data)[f, samples_to_keep]

    exp[exp == 0] <- 1

    print(head(exp))

  f_dat <- data.frame(fData(data)[f, ])
  colnames(f_dat) <- colnames(fData(data))
  p_dat <- data.frame(condition=pData(data)$condition[samples_to_keep], row.names=rownames(pData(data))[samples_to_keep])
  p_dat <- Biobase::AnnotatedDataFrame(p_dat)
  r <- Biobase::ExpressionSet(exp, p_dat)
  fData(r) <- f_dat

    print(head(exprs(r)))

  return(r)
}

#removing peptides that are detected in less than 2 samples per condition
#data <- msEmpiRe::filter_detection_rate(data, condition=conditions)
data = filderDetectionRate(data, conditions, rate=3)

print(data)


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