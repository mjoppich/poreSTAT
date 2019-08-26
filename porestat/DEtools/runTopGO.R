library(biomaRt)
library(org.Mm.eg.db)
library(topGO)
args = commandArgs(trailingOnly=TRUE)

filename = args[1]
mode=args[2]

print(paste("Running in mode", mode, sep=" "))


indata = read.table(filename, header=TRUE, sep="\t")

if (mode == "all")
{
    

    all_genes = indata[indata$ROB_ADJ.PVAL<=1.0,]$ROB_ADJ.PVAL
    names(all_genes) = indata[indata$ROB_ADJ.PVAL<=1.0,]$id

    print("All Mode - take all")
    diffGenes <- function(allScores) {
        return(allScores < 0.05)
    }

} else if (mode == "up")
{

    signThresh = 0.5

    print("Significance Threshold")
    print(signThresh)

    topDiffDownGenes <- function(allScores) { return(allScores >= signThresh)}


} else if (mode == "down")
{
    allGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL<=0.05 & indata$ROB_log2FC < 0,]$id)
    bgGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL>0.05 | (indata$ROB_ADJ.PVAL<=0.05 & indata$ROB_log2FC >= 0),]$id)


    print(length(allGeneIDs))
    print(length(bgGeneIDs))
}


#
# biocLite("clusterProfiler")

for (GODB in c("BP")) { #, "MF", "CC"

    print(head(all_genes))
    print(length(all_genes))

    topDiffDownGenes <- function(allScores) { return(allScores <= -signThresh)}

    GOdata = new("topGOdata", description="bla", ontology = GODB, allGenes=all_genes, geneSel=diffGenes, annot=annFUN.org, mapping = "org.Mm.eg.db", ID = "Ensembl", nodeSize=5)
    print(GOdata)
    #test.stat.ks = new("classicScore", testStatistic=GOKSTest, name="KS tests")
    #test.stat.count <- new("classicCount", testStatistic=GOFisherTest, name = "Fisher test")

    #downResultKS = runTest(GOdata, algorithm = "classic", statistic = "ks")
    #downResultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    ResultKS.weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
    ResultCount = runTest(GOdata, algorithm = "classic", statistic = "fisher")

    usedGOCount = length(usedGO(GOdata))
    print("used GO count")
    print(usedGOCount)

    #allres = GenTable(GOdata, topNodes = 100, classicKS=downResultKS, classicFisher=downResultCount, weight01KS=downResultKS.weight01, elimKS = downResultKS.elim, orderBy = "weight01KS")
    allres = GenTable(GOdata, topNodes = usedGOCount,  classicFisher=ResultCount, weight01KS=ResultKS.weight01, ranksOf="classicFisher", orderBy = "classicFisher")
    allres["classicFisher_ADJ"] = p.adjust(allres$classicFisher, method = "BH", n = usedGOCount)
    write.table(allres, file=paste(filename,"GeneOntology",GODB, "down.fisher.tsv", sep="."), sep="\t", quote=F, row.names=FALSE)

    allres = GenTable(GOdata, topNodes = usedGOCount,  classicFisher=ResultCount, weight01KS=ResultKS.weight01, ranksOf="weight01KS", orderBy = "weight01KS")
    write.table(allres, file=paste(filename,"GeneOntology",GODB,"down.weightks.tsv", sep="."), sep="\t", quote=F, row.names=FALSE)


}