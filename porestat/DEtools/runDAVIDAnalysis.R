suppressMessages(suppressWarnings(require(clusterProfiler)))
suppressMessages(suppressWarnings(require(annotables)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(RDAVIDWebService)))
suppressMessages(suppressWarnings(require(qvalue)))

args = commandArgs(trailingOnly=TRUE)



callDAVIDEnrichment <- function(gene,
                        idType        = "ENTREZ_GENE_ID",
                        universe,
                        minGSSize     = 10,
                        maxGSSize     = 500,
                        annotation    = "GOTERM_BP_FAT",
                        pvalueCutoff  = 0.05,
                        pAdjustMethod = "BH",
                        qvalueCutoff  = 0.2,
                        species       = NA,
                        david.user){

    Count <- List.Total <- Pop.Hits <- Pop.Total <- NULL

    pAdjustMethod <- match.arg(pAdjustMethod, c("bonferroni", "BH"))

    david.pkg <- "RDAVIDWebService"
    pkgs <- installed.packages()[,1]
    if (! david.pkg %in% pkgs) {
        stop("You should have RDAVIDWebService package installed before using enrichDAVID...")
    }

    require(david.pkg, character.only=TRUE)
    DAVIDWebService <- eval(parse(text="DAVIDWebService"))
    addList <- eval(parse(text="addList"))
    setAnnotationCategories <- eval(parse(text="setAnnotationCategories"))
    getFunctionalAnnotationChart <- eval(parse(text="getFunctionalAnnotationChart"))
    getSpecieNames <- eval(parse(text="getSpecieNames"))
    getIdTypes <- eval(parse(text="getIdTypes"))

    david <- DAVIDWebService$new(email=david.user,
                                 url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

    setTimeOut(david, 120000)

    ## addList will throw error if idType is not match.
    ## use match.arg to check before addList make it more readable
    
    idType <- match.arg(idType, getIdTypes(david))
    
    ##     getIdTypes(david)
    ##  [1] "AFFYMETRIX_3PRIME_IVT_ID" "AFFYMETRIX_EXON_ID"      
    ##  [3] "AGILENT_CHIP_ID"          "AGILENT_ID"              
    ##  [5] "AGILENT_OLIGO_ID"         "APHIDBASE_ID"            
    ##  [7] "BEEBASE_ID"               "BEETLEBASE_ID"           
    ##  [9] "BGD_ID"                   "CGNC_ID"                 
    ## [11] "CRYPTODB_ID"              "DICTYBASE_ID"            
    ## [13] "ENSEMBL_GENE_ID"          "ENSEMBL_TRANSCRIPT_ID"   
    ## [15] "ENTREZ_GENE_ID"           "FLYBASE_GENE_ID"         
    ## [17] "GENBANK_ACCESSION"        "GENOMIC_GI_ACCESSION"    
    ## [19] "GENPEPT_ACCESSION"        "LOCUS_TAG"               
    ## [21] "MGI_ID"                   "MIRBASE_ID"              
    ## [23] "MRNA_GI_ACCESSION"        "NASONIABASE_ID"          
    ## [25] "PROTEIN_GI_ACCESSION"     "PSEUDOCAP_ID"            
    ## [27] "REFSEQ_MRNA"              "REFSEQ_PROTEIN"          
    ## [29] "RGD_ID"                   "SGD_ID"                  
    ## [31] "TAIR_ID"                  "UNIGENE"                 
    ## [33] "UNIPROT_ACCESSION"        "UNIPROT_ID"              
    ## [35] "VECTORBASE_ID"            "WORMBASE_GENE_ID"        
    ## [37] "XENBASE_ID"               "ZFIN_ID"
    
    david.res <- addList(david, gene, idType=idType,
                         listName="clusterProfiler",
                         listType="Gene")


    if (david.res$inDavid == 0) {
        stop("All id can not be mapped. Please check 'idType' parameter...")
    }

    if (!missing(universe)) {
        david.res <- addList(david, universe, idType=idType,
                             listName="universe",
                             listType="Background")
    }

    setAnnotationCategories(david, annotation)
    x <- getFunctionalAnnotationChart(david, threshold=1, count=minGSSize)

    if (length(x@.Data) == 0) {
        warning("No significant enrichment found...")
        return(NULL)
    }

    term <- x$Term
    sep <- "~"

    term.list <- sapply(term, function(y) strsplit(y, split=sep))
    term.df <- do.call("rbind", term.list)


    ID <- term.df[,1]
    Description <- term.df[,2]
    GeneRatio <- with(x, paste(Count, List.Total, sep="/"))
    BgRatio <- with(x, paste(Pop.Hits, Pop.Total, sep="/"))
    Over <- data.frame(ID          = ID,
                       Description = Description,
                       GeneRatio   = GeneRatio,
                       BgRatio     = BgRatio,
                       pvalue      = x$PValue,
                       Category = x$Category,
                       FoldEnrichment = x$Fold.Enrichment,
                       stringsAsFactors = FALSE)


    row.names(Over) <- term

    if (pAdjustMethod == "bonferroni") {
        Over$p.adjust <- x$Bonferroni
    } else {
        Over$p.adjust <- x$Benjamini
    }

    qobj <- tryCatch(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"),error=function(e) print(e))
    if (class(qobj) == "qvalue") {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
    }
    Over$qvalue <- qvalues
    Over$geneID <- gsub(",\\s*", "/", x$Genes)
    Over$Count <- x$Count

    Over <- Over[ Over$pvalue <= pvalueCutoff, ]
    Over <- Over[ Over$p.adjust <= pvalueCutoff, ]
    if (! any(is.na(Over$qvalue))) {
        Over <- Over[ Over$qvalue <= qvalueCutoff, ]
    }

    org <- getSpecieNames(david)
    org <- gsub("\\(.*\\)", "", org)

    ## gc <- strsplit(Over$geneID, "/")
    ## names(gc) <- Over$ID

    if (!is.na(maxGSSize) && !is.null(maxGSSize)) {
        idx <- as.numeric(sub("/\\d+", "", Over$BgRatio)) <= maxGSSize
        Over <- Over[idx,]
    }

    new("enrichResult",
        result         = Over,
        pvalueCutoff   = pvalueCutoff,
        pAdjustMethod  = pAdjustMethod,
        organism       = org,
        ontology       = annotation, ## as.character(x$Category[1]),
        gene           = as.character(gene),
        keytype        = idType)
}






filename = args[1]
organism = args[2]

mode=args[3]

print(paste("Running in mode", mode, sep=" "))

allGeneIDs = NULL
bgGeneIDs = NULL
indata = read.table(filename, header=TRUE, sep="\t")

minFC = 1.0

if (mode == "all")
{
    allGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL<=0.05 & abs(indata$ROB_log2FC) > minFC,]$id)
    bgGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL>0.05 | (indata$ROB_ADJ.PVAL<=0.05 & abs(indata$ROB_log2FC) <= minFC),]$id)

    print(length(allGeneIDs))
    print(length(bgGeneIDs))
} else if (mode == "up")
{
    allGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL<=0.05 & indata$ROB_log2FC > minFC,]$id)
    bgGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL>0.05 | (indata$ROB_ADJ.PVAL<=0.05 & indata$ROB_log2FC <= minFC),]$id)


    print(length(allGeneIDs))
    print(length(bgGeneIDs))


} else if (mode == "down")
{
    allGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL<=0.05 & indata$ROB_log2FC < -minFC,]$id)
    bgGeneIDs = as.vector(indata[indata$ROB_ADJ.PVAL>0.05 | (indata$ROB_ADJ.PVAL<=0.05 & indata$ROB_log2FC >= -minFC),]$id)


    print(length(allGeneIDs))
    print(length(bgGeneIDs))
}


if (length(allGeneIDs) == 0)
{
    print("NO GENE ID SELECTED. TERMINATING")
    quit(status=0, save='no')
}


egid = grcm38 %>% dplyr::filter(ensgene %in% allGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()
head(egid)

entrezGenes = egid[!is.na(egid$entrez),]
entrezGenes = as.vector(entrezGenes$entrez)



egid = grcm38 %>% dplyr::filter(ensgene %in% bgGeneIDs) %>% dplyr::select(ensgene, entrez) %>% as.data.frame()
head(egid)

entrezBG = egid[!is.na(egid$entrez),]
entrezBG = as.vector(entrezBG$entrez)

entrezGenesC = unlist(lapply(entrezGenes, as.character))
entrezBGC = unlist(lapply(entrezBG, as.character))

head(entrezGenesC)
head(entrezBGC)

#like 
#https://www.r-bloggers.com/kegg-enrichment-analysis-with-latest-online-data-using-clusterprofiler/

#source("https://bioconductor.org/biocLite.R")
#biocLite("RDAVIDWebService")

annotations = c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT","REACTOME_PATHWAY","KEGG_PATHWAY", "UP_KEYWORDS")

kk = callDAVIDEnrichment(gene=entrezGenesC, universe=c(entrezGenesC, entrezBGC), species=organism, idType="ENTREZ_GENE_ID", annotation = annotations, david.user="joppich@bio.ifi.lmu.de", pAdjustMethod="BH", pvalueCutoff=0.5, qvalueCutoff=0.5)

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
rsc[1] = "TERM.ID"
colnames(rs) = rsc

write.table(rs, file=paste(filename,"david",mode,".tsv", sep="."), sep="\t", quote=F, row.names=FALSE)