
suppressMessages(suppressWarnings(require(clusterProfiler)))
suppressMessages(suppressWarnings(require(annotables)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(RDAVIDWebService)))
suppressMessages(suppressWarnings(require(qvalue)))

args = commandArgs(trailingOnly=TRUE)


filename = args[1]
organism = args[2]
mode=args[3]

print(paste("Running in mode", mode, sep=" "))

allGeneIDs = NULL
bgGeneIDs = NULL
indata = read.table(filename, header=TRUE, sep="\t")

minFC = 1.0

allGenes = indata
allGenes = allGenes[allGenes$ROB_ADJ.PVAL != 1.0 & allGenes$ROB_log2FC != 0.0 & !is.na(allGenes$id) ,]
allGenes = allGenes[!is.na(allGenes$id),]
allGeneIDs = as.vector(allGenes$id)


if (nrow(allGenes) == 0)
{
    print("NO GENE ID SELECTED. TERMINATING")
    quit(status=0, save='no')
}


annotTable = NULL;
readableState=T
keyType = "ENTREZID"

if (organism == "org.Mm.eg.db")
{
    annotTable = grcm38
    egid = annotTable %>% dplyr::filter(ensgene %in% allGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()


} else if (organism == "org.Hs.eg.db")
{
    annotTable = grch38
    egid = annotTable %>% dplyr::filter(ensgene %in% allGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()

} else if (organism == "org.Sc.sgd.db")
{

    suppressMessages(suppressWarnings(require(org.Sc.sgd.db)))


    allEntrez = sapply(allGeneIDs, function(x) {get(x, org.Sc.sgdENTREZID)} )
    egid = data.frame("ensgene"=allGeneIDs, "entrez"=allEntrez)

    keyType = "GENENAME"

}

#egid = grcm38 %>% dplyr::filter(ensgene %in% allGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()
#head(egid)

print("Dims of egid")
print(dim(egid))

egid = merge(x=egid, y=allGenes, by.x="ensgene", by.y="id")

tdf = data.frame(etz=egid$entrez, ens=egid$ensgene, lfc=egid$ROB_log2FC, pval=egid$ROB_ADJ.PVAL)
print(tdf)

entrezGenes = egid[!is.na(egid$entrez) & !is.na(egid$ensgene) & !is.null(egid$ensgene) & !is.null(egid$ROB_log2FC)& !is.nan(egid$ROB_log2FC) &!is.na(egid$ROB_log2FC),]

#entrezGenesC = entrezGenes$ROB_ADJ.PVAL
entrezGenesC = as.vector(entrezGenes$ROB_log2FC)


if (keyType == "GENENAME")
{
    names(entrezGenesC) = entrezGenes$ensgene
} else {
    names(entrezGenesC) = entrezGenes$entrez
}

entrezGenesC = sort(entrezGenesC, decreasing = TRUE)

print("entrez genes c")
head(entrezGenesC)

entrezGenesC <- entrezGenesC[!duplicated(names(entrezGenesC))]

print("num entrez genes c")
print(length(entrezGenesC))
print("Dup genes")
print(length(entrezGenesC[duplicated(names(entrezGenesC))]))

#like
#https://www.r-bloggers.com/kegg-enrichment-analysis-with-latest-online-data-using-clusterprofiler/

for (GODB in c("BP", "MF", "CC")) { #


    #by DOSE as fgsea has problems with all genes being in one class...
    kk <- gseGO(entrezGenesC, organism, keyType=keyType, ont=GODB, pvalueCutoff=0.5, pAdjustMethod="BH", by="DOSE")

    if (is.null(kk))
    {
        print("NO DATA RETURNED. TERMINATING")
        quit(status=0, save='no')
    }

    rs = as.data.frame(kk)

    if (nrow(rs) == 0)
    {
        print("NO DATA in DF. TERMINATING")
        quit(status=0, save='no')
    }

    rsc = colnames(rs)
    rsc[1] = "GO.ID"
    colnames(rs) = rsc

    if (match('qvalues',colnames(rs)))
    {
        names(rs)[names(rs)=="qvalues"] <- "qvalue"
    }

    mode="all"
    write.table(rs, file=paste(filename,"GeneOntology", GODB, mode,"gsea.tsv", sep="."), sep="\t", quote=F, row.names=FALSE)

    mode="down"
    rsd = rs[rs$NES < 0,]
    write.table(rsd, file=paste(filename,"GeneOntology", GODB, mode,"gsea.tsv", sep="."), sep="\t", quote=F, row.names=FALSE)

    mode="up"
    rsd = rs[rs$NES > 0,]
    write.table(rsd, file=paste(filename,"GeneOntology", GODB, mode,"gsea.tsv", sep="."), sep="\t", quote=F, row.names=FALSE)
}

print("finished")
quit(status=0, save='no')