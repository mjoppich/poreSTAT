
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
allGeneIDs = allGenes$id

allGenes = allGenes[allGenes$ROB_ADJ.PVAL != 1.0 & allGenes$ROB_log2FC != 0.0,]

if (nrow(allGenes) == 0)
{
    print("NO GENE ID SELECTED. TERMINATING")
    quit(status=0, save='no')
}

egid = grcm38 %>% dplyr::filter(ensgene %in% allGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()
head(egid)

egid = merge(x=egid, y=allGenes, by.x="ensgene", by.y="id")
head(egid)

entrezGenes = egid[!is.na(egid$entrez) & !is.null(egid$ROB_log2FC)& !is.nan(egid$ROB_log2FC) &!is.na(egid$ROB_log2FC),]
head(entrezGenes)
print(sum(is.na(entrezGenes$ROB_log2FC)))

#entrezGenesC = entrezGenes$ROB_ADJ.PVAL
entrezGenesC = as.vector(entrezGenes$ROB_log2FC)
names(entrezGenesC) = entrezGenes$entrez

entrezGenesC = sort(entrezGenesC, decreasing = TRUE)

head(entrezGenesC)


#like
#https://www.r-bloggers.com/kegg-enrichment-analysis-with-latest-online-data-using-clusterprofiler/

for (GODB in c("BP", "MF", "CC")) { #


    kk <- gseGO(entrezGenesC, organism, keyType="ENTREZID", ont=GODB, pvalueCutoff=0.5, pAdjustMethod="BH")

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

    mode="all"
    write.table(rs, file=paste(filename,"GeneOntology", GODB, mode,"gsea.tsv", sep="."), sep="\t", quote=F, row.names=FALSE)

    mode="down"
    rsd = rs[rs$NES < 0,]
    write.table(rsd, file=paste(filename,"GeneOntology", GODB, mode,"gsea.tsv", sep="."), sep="\t", quote=F, row.names=FALSE)

    mode="up"
    rsd = rs[rs$NES > 0,]
    write.table(rsd, file=paste(filename,"GeneOntology", GODB, mode,"gsea.tsv", sep="."), sep="\t", quote=F, row.names=FALSE)
}