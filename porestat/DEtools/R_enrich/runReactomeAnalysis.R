suppressMessages(suppressWarnings(require(ReactomePA)))
#suppressMessages(suppressWarnings(require(clusterProfiler)))
suppressMessages(suppressWarnings(require(annotables)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(qvalue)))
suppressMessages(suppressWarnings(require(writexl)))

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

if (organism == "mouse")
{
    annotTable = grcm38
    egid = annotTable %>% dplyr::filter(ensgene %in% allGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()


} else if (organism == "human")
{
    annotTable = grch38
    egid = annotTable %>% dplyr::filter(ensgene %in% allGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()

} else if (organism == "yeast")
{

    suppressMessages(suppressWarnings(require(org.Sc.sgd.db)))

    readableState =F
    allEntrez = sapply(allGeneIDs, function(x) {get(x, org.Sc.sgdENTREZID)} )
    egid = data.frame("ensgene"=allGeneIDs, "entrez"=allEntrez)

}

#print("got egid")
#head(egid)

entrezGenes = egid[!is.na(egid$entrez),]
entrezGenes = as.vector(entrezGenes$entrez)


annotTable = NULL;
readableState = T

if (organism == "mouse")
{
    annotTable = grcm38
    egid = annotTable %>% dplyr::filter(ensgene %in% bgGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()


} else if (organism == "human")
{
    annotTable = grch38
    egid = annotTable %>% dplyr::filter(ensgene %in% bgGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()

} else if (organism == "yeast")
{

    suppressMessages(suppressWarnings(require(org.Sc.sgd.db)))
    readableState =F
    allEntrez = sapply(bgGeneIDs, function(x) {get(x, org.Sc.sgdENTREZID)} )
    egid = data.frame("ensgene"=bgGeneIDs, "entrez"=allEntrez)
}

print("got egid")
head(egid)

#egid = grcm38 %>% dplyr::filter(ensgene %in% bgGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()
#head(egid)

entrezBG = egid[!is.na(egid$entrez),]
entrezBG = as.vector(entrezBG$entrez)

entrezGenesC = unlist(lapply(entrezGenes, as.character))
entrezBGC = unlist(lapply(entrezBG, as.character))

head(entrezGenesC)
head(entrezBGC)

#like 
#https://www.r-bloggers.com/kegg-enrichment-analysis-with-latest-online-data-using-clusterprofiler/
 
kk <- enrichPathway(entrezGenesC, organism=organism, pAdjustMethod="BH", universe=c(entrezGenesC, entrezBGC), pvalueCutoff=0.5, qvalueCutoff=0.5, readable=readableState)

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
rsc[1] = "Reactome ID"
colnames(rs) = rsc

write.table(rs, file=paste(filename,"reactome",mode,"tsv", sep="."), sep="\t", quote=F, row.names=FALSE)
write_xlsx(rs, path=paste(filename,"reactome", mode,"xlsx", sep="."))

print("finished")
quit(status=0, save='no')