
suppressMessages(suppressWarnings(require(clusterProfiler)))
suppressMessages(suppressWarnings(require(annotables)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(RDAVIDWebService)))
suppressMessages(suppressWarnings(require(qvalue)))

suppressMessages(suppressWarnings(require(org.Hs.eg.db)))
suppressMessages(suppressWarnings(require(org.Mm.eg.db)))



args = commandArgs(trailingOnly=TRUE)

filename = args[1]
organism = args[2]
mode=args[3]
minPVal = as.numeric(args[4])
minFC = as.numeric(args[5])

print(paste("Running in mode", mode, sep=" "))

allGeneIDs = NULL
bgGeneIDs = NULL
indata = read.table(filename, header=TRUE, sep="\t")

minFC = 1.0

if (mode == "all")
{
    allGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL<=minPVal & abs(indata$ROB_log2FC) > minFC,]$id)
    bgGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL>minPVal | (indata$ROB_ADJ.PVAL<=minPVal & abs(indata$ROB_log2FC) <= minFC),]$id)

    print(length(allGeneIDs))
    print(length(bgGeneIDs))
} else if (mode == "up")
{
    allGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL<=minPVal & indata$ROB_log2FC > minFC,]$id)
    bgGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL>minPVal | (indata$ROB_ADJ.PVAL<=minPVal & indata$ROB_log2FC <= minFC),]$id)


    print(length(allGeneIDs))
    print(length(bgGeneIDs))


} else if (mode == "down")
{
    allGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL<=minPVal & indata$ROB_log2FC < -minFC,]$id)
    bgGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL>minPVal | (indata$ROB_ADJ.PVAL<=minPVal & indata$ROB_log2FC >= -minFC),]$id)


    print(length(allGeneIDs))
    print(length(bgGeneIDs))
}

allGeneIDs = allGeneIDs[!is.na(allGeneIDs)]
bgGeneIDs = bgGeneIDs[!is.na(bgGeneIDs)]

if (length(allGeneIDs) == 0)
{
    print("NO GENE ID SELECTED. TERMINATING")
    quit(status=0, save='no')
}

# just in case we got stupid ensembl stable IDs ... ENSG00012817271.x => ENSG00012817271
allGeneIDs = gsub("\\..*","",allGeneIDs)
bgGeneIDs = gsub("\\..*","",bgGeneIDs)

annotTable = NULL;

if (organism == "mmu")
{
    annotTable = grcm38
    conversionTable = org.Mm.eg.db
    egid = annotTable %>% dplyr::filter(ensgene %in% allGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()


} else if (organism == "hsa")
{
    annotTable = grch38
    conversionTable = org.Hs.eg.db

    egid = annotTable %>% dplyr::filter(ensgene %in% allGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()

} else if (organism == "sce")
{

    suppressMessages(suppressWarnings(require(org.Sc.sgd.db)))
    conversionTable = org.Sc.sgd.db

    readableState =F
    allEntrez = allGeneIDs#sapply(allGeneIDs, function(x) {get(x, org.Sc.sgdENTREZID)} )
    egid = data.frame("ensgene"=allGeneIDs, "entrez"=allEntrez)

}

entrezGenes = egid[!is.na(egid$entrez),]
entrezGenes = as.vector(entrezGenes$entrez)

annotTable = NULL;

if (organism == "mmu")
{
    annotTable = grcm38
    egid = annotTable %>% dplyr::filter(ensgene %in% bgGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()


} else if (organism == "hsa")
{
    annotTable = grch38
    egid = annotTable %>% dplyr::filter(ensgene %in% bgGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()

} else if (organism == "sce")
{

    suppressMessages(suppressWarnings(require(org.Sc.sgd.db)))

    readableState =F
    allEntrez = bgGeneIDs#sapply(bgGeneIDs, function(x) {get(x, org.Sc.sgdENTREZID)} )
    egid = data.frame("ensgene"=bgGeneIDs, "entrez"=allEntrez)

}

#head(egid)

entrezBG = egid[!is.na(egid$entrez),]
entrezBG = as.vector(entrezBG$entrez)

entrezGenesC = unlist(lapply(entrezGenes, as.character))
entrezBGC = unlist(lapply(entrezBG, as.character))

head(entrezGenesC)
head(entrezBGC)

#like
#https://www.r-bloggers.com/kegg-enrichment-analysis-with-latest-online-data-using-clusterprofiler/
 
kk <- enrichKEGG(entrezGenesC, organism=organism,pAdjustMethod="BH",universe=c(entrezGenesC, entrezBGC), pvalueCutoff=0.5, qvalueCutoff=0.5)

if (is.null(kk))
{
    print("NO DATA RETURNED. TERMINATING")
    quit(status=0, save='no')
}

kk <- setReadable(kk, OrgDb = conversionTable, keyType="ENTREZID")

rs = as.data.frame(kk)

if (nrow(rs) == 0)
{
    print("NO DATA in DF. TERMINATING")
     quit(status=0, save='no')
}

rsc = colnames(rs)
rsc[1] = "KEGG ID"
colnames(rs) = rsc

write.table(rs, file=paste(filename,"kegg",mode,"tsv", sep="."), sep="\t", quote=F, row.names=FALSE)

print("finished")
quit(status=0, save='no')