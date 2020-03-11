from collections import defaultdict, Counter

from matplotlib_venn._venn3 import compute_venn3_subsets

from porestat.utils.Numbers import toNumber, toFloat

from ..analysis.similarity_analysis import SimilarityAnalysis
from .DataFrame import DataFrame, DataRow, ExportTYPE
import os
import mpld3
from ..plots.plotconfig import PlotConfig, PlotSaveTYPE
from ..plots.poreplot import PorePlot
import pickle
from ..tools.PTToolInterface import PSToolException


class EnrichmentDF(DataFrame):

    def __init__(self, otherDF=None):
        super(EnrichmentDF, self).__init__()
        self.idCol = self.addColumn('id')

        if otherDF:
            self.data = otherDF.data
            self.column2idx = otherDF.column2idx


    def addConditions(self, condDatas, condName):

        allCols = set()
        for x in condDatas:
            for y in x:
                allCols.add(y)

        self.addColumns(sorted(allCols), default=None, ignoreDuplicates=True)

        self.updateRowIndexed("id", condDatas, ignoreMissingCols=True, addIfNotFound=True)


    def addCondition(self, condData, condName):

        isDataRow = isinstance(condData, DataRow)
        isDict = isinstance(condData, dict)

        if not (isDataRow or isDict):
            raise ValueError("parameter not of type datarow or dict")

        if isDataRow:
            return self._addConditionFromDataRow(condData, condName)

        if isDict:
            return self._addConditionFromDict(condData, condName)


    def _addConditionFromDict(self, condData, condName):
        dr = DataRow.fromDict(condData)
        return self._addConditionFromDataRow(dr, condName)


    def _addConditionFromDataRow(self, condData, condName):

        newColIdx = self.addColumn(condName, None)

        setFoundKeys = set()
        cDataHeader = condData.getHeader()

        for i in range(0, len(self.data)):
            row = self.data[i]
            rowID = row[ self.idCol ]

            if rowID in cDataHeader:
                lrow = list(row)
                lrow[newColIdx] = condData[rowID]
                self.data[i] = tuple(lrow)
                setFoundKeys.add(rowID)

        for geneID in cDataHeader:
            if not geneID in setFoundKeys:

                lrow = [self.idx2default[i] for i in range(0, len(self.column2idx))]
                lrow[0] = geneID
                lrow[newColIdx] = condData[geneID]

                self.data.append(tuple(lrow))


    def writeEnrichmentBrowserFiles(self, prepData, replicates, exprFile, pdataFile, fdataFile):

        self.writeExpressionFile(exprFile, prepData)
        self.writepDataFile(pdataFile, replicates)
        self.writefDataFile(fdataFile, prepData)



    def prepareEBData(self, cond1Samples, cond2Samples):

        assert (type(cond1Samples) == list)
        assert (type(cond2Samples) == list)

        print(cond1Samples)
        print(cond2Samples)

        header = ['gene'] + cond1Samples + cond2Samples
        prepData = [tuple(header)]


        condSample2Idx = {}

        for cond1Sample in cond1Samples:
            condSample2Idx[cond1Sample] = self.getColumnIndex(cond1Sample)

        for cond2Sample in cond2Samples:
            condSample2Idx[cond2Sample] = self.getColumnIndex(cond2Sample)


        for row in self:

            gene = row[0]
            baseRow = [gene]

            for sampleName in cond1Samples+cond2Samples:

                count = int(row[condSample2Idx[sampleName]])

                if count != None:
                    baseRow.append(count)

            if len(baseRow) == len(cond1Samples) + len(cond2Samples) + 1:
                prepData.append( tuple(baseRow) )

        return prepData


    def writepDataFile(self, file, replicates):

        with open(file, 'w') as pdataFile:

            for ridx, replicate in enumerate(replicates):
                for sample in replicates[replicate]:
                    pdataFile.write(sample + "\t" + str(ridx) + "\n")


    def writefDataFile(self, file, prepData):

        with open(file, 'w') as fdataFile:
            sampleCount = len(prepData[0])
            for line in prepData[1:]:
                fdataFile.write( "\t".join( [line[0]] * sampleCount) + "\n")

    def writeExpressionFile(self, file, prepData):

        with open(file, 'w') as exprFile:

            for line in prepData[1:]:
                exprFile.write( "\t".join([str(x) for x in line[1:]]) + "\n")



    @property
    def supported_de_methods(self):
        return ['NOISeq', 'DESeq2', 'DirectDESeq2', 'msEmpiRe', 'limma', 'edgeR']


    def runDEanalysis(self, outputFolder, replicates, prefix= "", methods=['NOISeq', 'msEmpiRe', 'DESeq2', "DirectDESeq2"], rscriptPath="/usr/bin/Rscript", noDErun=False, enhanceSymbol=None, geneLengths=None, norRNA=False):


        filePrefix = prefix
        if prefix != None and prefix != "" and prefix[len(prefix)-1] != "_":
            filePrefix = prefix + "_"

        print("Running DE analysis with prefix " + str(prefix) + " and methods " + str(methods))

        basePath = outputFolder

        if not basePath[-1] == '/':
            basePath += "/"

        base = basePath + filePrefix

        exprFile = base + "expr"
        pdataFile = base + "p_data"
        fdataFile = base + "f_data"
        outFileBase = base + "out_data"

        noiseqFile = base + "noiseq"

        for x in methods:
            if not x in self.supported_de_methods:
                raise PSToolException("Unsupported DE Method in runDEAnalysis: " + str(x) + "\nAvailable Methods: " + ", ".join(self.supported_de_methods))


        createdComparisons = []
        conditions = [x for x in replicates]

        for cond in replicates:
            print(cond, replicates[cond])

        for i in range(0, len(conditions)):
            cond1 = conditions[i]
            for j in range(i+1, len(conditions)):
                cond2 = conditions[j]

                prepData = self.prepareEBData(replicates[cond1], replicates[cond2])

                with open(noiseqFile, 'w') as fout:

                    for elem in prepData:
                        fout.write("\t".join([elem[0]] + [str(x) for x in reversed(elem[1:])])  +"\n")

                self.writeEnrichmentBrowserFiles(prepData, replicates, exprFile, pdataFile, fdataFile)

                condResult = {}

                for method in methods: # 'limma' 'edgeR'
                    outFile = outFileBase + "_" + method

                    execStr = None
                    if method in ['DESeq2', 'limma', 'edgeR']:
                        scriptPath = os.path.dirname(os.path.abspath(__file__)) + "/../data/de_rseq.R"
                        execStr = rscriptPath+" "+scriptPath+" "+exprFile+" "+pdataFile+" "+fdataFile+" "+method+" " + outFile
                            #raise PSToolException("Unable to run enrichment analsyis. " + str(execStr))

                    elif method in ['DirectDESeq2']:
                        scriptPath = os.path.dirname(os.path.abspath(__file__)) + "/../data/deseq_direct.R"
                        execStr = rscriptPath+" "+scriptPath+" "+exprFile+" "+pdataFile+" "+fdataFile+" " + outFile

                    elif method in ['msEmpiRe']:
                        scriptPath = os.path.dirname(os.path.abspath(__file__)) + "/../data/empire_diffreg.R"

                        execStr = rscriptPath+" "+scriptPath+" "+exprFile+" "+pdataFile+" "+noiseqFile+" "+method+" " + outFile

                    elif method in ['NOISeq']:

                        sample4condition = [conditions[i]] * len(replicates[cond2]) + [conditions[j]] * len(replicates[cond1])
                        print(sample4condition)

                        scriptPath = os.path.dirname(os.path.abspath(__file__)) + "/../data/noiseq_diffreg.R"

                        execStr = rscriptPath+" "+scriptPath+" "+noiseqFile + " "+ outFile + " " + " ".join(sample4condition)

                    print(execStr)

                    if not noDErun:
                        sysret = os.system(execStr)
                    sysret = 1

                    if sysret != 0:
                        pass

                    methDF = DataFrame.parseFromFile(outFile)
                    condResult[method] = methDF


                # TODO this must be done for all samples!
                compDF = EnrichmentDF()
                #id column already there!
                evCol = compDF.addColumn('evidence', prefix) # TODO this is just a quick hack!


                geneNames = self.getColumnIndex('id')

                updateRows = []

                for row in self:
                    baseDict = {"id": row["id"]}

                    for sampleName in replicates[cond1] + replicates[cond2]:
                        baseDict[sampleName] = row[sampleName]

                        if sampleName + ".FPKM" in self.column2idx:
                            baseDict[sampleName + ".FPKM"] = row[sampleName + ".FPKM"]

                        if sampleName + ".TPM" in self.column2idx:
                            baseDict[sampleName + ".TPM"] = row[sampleName + ".TPM"]

                    updateRows.append(baseDict)

                    #condSampleCounts = self.toDataRow(geneNames, self.getColumnIndex(sampleName))

                allCols = set()
                for x in updateRows:
                    for y in x:
                        allCols.add(y)

                compDF.addColumns(sorted(allCols), default=0, ignoreDuplicates=True)
                compDF.updateRowIndexed("id", updateRows, ignoreMissingCols=True, addIfNotFound=True)
                    #compDF.addCondition(condSampleCounts, sampleName)



                condVPData = {}
                usedMethod = []
                for method in condResult:

                    methDF = condResult[method]

                    if methDF == None:
                        continue

                    print("Combining Method", method)

                    usedMethod.append(method)

                    l2FCTitle = method + "_log2FC"
                    rawpTitle = method + "_RAW.PVAL"
                    adjpTitle = method + "_ADJ.PVAL"

                    def recognizeColumnIndex( df, colNames):

                        for colName in colNames:
                            if df.columnExists(colName):
                                return df.getColumnIndex(colName)

                        return -1


                    geneIDidx = recognizeColumnIndex(methDF, ["GENE.ID", "PROBEID"])#methDF.getColumnIndex('GENE.ID') if methDF.columnExists('GENE.ID') else methDF.getColumnIndex('PROBEID')
                    log2FCidx = recognizeColumnIndex(methDF, ["log2FC", "FC"]) #methDF.getColumnIndex('log2FC') if methDF.columnExists('log2FC') else methDF.getColumnIndex('FC')
                    rawPValidx = recognizeColumnIndex(methDF, ["RAW.PVAL", "PVAL", "ADJ.PVAL"]) # methDF.getColumnIndex('RAW.PVAL') if methDF.columnExists('RAW.PVAL') else methDF.getColumnIndex('ADJ.PVAL')
                    adjPValidx = methDF.getColumnIndex('ADJ.PVAL')

                    print("Combining Method to data row", method)
                    l2FCdata = methDF.toDataRow( geneIDidx, log2FCidx)
                    rawPdata = methDF.toDataRow( geneIDidx, rawPValidx)
                    adjPdata = methDF.toDataRow( geneIDidx, adjPValidx)

                    noiSeqProbIdx = None
                    noiSeqProbTitle = None
                    if method in ['NOISeq']:
                        noiSeqProbTitle = method + "_prob"
                        if methDF.columnExists("prob"):
                            noiSeqProbIdx = methDF.getColumnIndex('prob')


                    methRows = []
                    for row in methDF:

                        rowDict = {
                            "id": row[geneIDidx],
                            l2FCTitle: row.getIndex(log2FCidx),
                            rawpTitle: row.getIndex(rawPValidx),
                            adjpTitle: row.getIndex(adjPValidx),
                        }

                        if noiSeqProbIdx != None:
                            rowDict[noiSeqProbTitle] = row[noiSeqProbIdx]


                        methRows.append(rowDict)


                    print("Combining Method adding cond", method)
                    print("Adding total rows:", len(methRows))
                    allCols = set()
                    for x in methRows:
                        for y in x:
                            allCols.add(y)

                    compDF.addColumns(sorted(allCols), default=None, ignoreDuplicates=True)
                    compDF.updateRowIndexed("id", methRows, ignoreMissingCols=True, addIfNotFound=True)

                    #compDF.addCondition(l2FCdata, l2FCTitle)
                    #compDF.addCondition(rawPdata, rawpTitle)
                    #compDF.addCondition(adjPdata, adjpTitle)


                    print("cond VP Data ready")
                    condVPData[method] = ( l2FCdata.getHeader(), l2FCdata.to_list(), rawPdata.to_list(), adjPdata.to_list() )


                def addInfoFunc(x):
                    gene = x[0]
                    x[evCol] = prefix

                    link = "<a target='_blank' href='http://bacteria.ensembl.org/Helicobacter_pylori_p12/Gene/Summary?g="+gene+"'>EnsemblBacteria</a><br/>" \
                           "<a target='_blank' href='http://www.uniprot.org/uniprot/?query="+gene+"&sort=score'>Uniprot</a>"

                    if not gene.startswith("HP"):
                        link = ""

                    x[addInfo] = link
                    return tuple(x)

                if enhanceSymbol != None:

                    def addGeneSymbol(x, gscol):
                        gene = x[0]
                        genesym = enhanceSymbol.get(gene, (gene, ""))

                        geneSymbol = genesym[0]

                        if geneSymbol == None or geneSymbol == "":
                            geneSymbol = gene

                        x[gscol] = geneSymbol

                        return tuple(x)


                    gsCol = compDF.addColumn("gene_symbol")
                    compDF.applyToRow( lambda x: addGeneSymbol(x, gsCol) )

                    def addGeneBiotype(x, gscol):
                        gene = x[0]
                        genesym = enhanceSymbol.get(gene, ("", ""))

                        x[gscol] = genesym[1]

                        return tuple(x)


                    gsCol = compDF.addColumn("gene_biotype")
                    compDF.applyToRow( lambda x: addGeneBiotype(x, gsCol) )

                if geneLengths != None:

                    def addGeneLengths(x, gscol):
                        gene = x[0]
                        genesym = geneLengths.get(gene, "")

                        x[gscol] = genesym

                        return tuple(x)


                    gsCol = compDF.addColumn("ue_gene_length")
                    compDF.applyToRow( lambda x: addGeneLengths(x, gsCol) )


                #print("Applying add info")
                #addInfo = compDF.addColumn("EnsemblBacteria")
                #compDF.applyToRow( addInfoFunc )
                
                cond1FName = cond1.replace("/", "_").replace("\\", "_").replace(".", "_")
                cond2FName = cond2.replace("/", "_").replace("\\", "_").replace(".", "_")

                outname = base + cond1FName + "_" + cond2FName + ".tsv"

                createdComparisons.append( (cond1, cond2, outname) )

                compDF.export(outFile=outname, exType=ExportTYPE.TSV)

        return createdComparisons


    def printResult(self, outputFolder, prefix, conditionPair2File, replicates):

        filePrefix = prefix
        if prefix != None and prefix != "" and prefix[len(prefix)-1] != "_":
            filePrefix = prefix + "_"

        print("Running DE analysis with prefix " + str(prefix))

        basePath = outputFolder

        if not basePath[-1] == '/':
            basePath += "/"

        base = basePath + filePrefix

        import mpld3, shutil
        d3js_path = mpld3.getD3js()
        mpld3_path = mpld3.getmpld3js(True)

        d3js_dest = outputFolder + os.path.split(d3js_path)[1]
        mpld3_dest = outputFolder + os.path.split(mpld3_path)[1]

        shutil.copyfile(d3js_path, d3js_dest)
        shutil.copyfile(mpld3_path, mpld3_dest)

        fileOutPuts = {}
        allMethods = set()

        for cond1, cond2 in conditionPair2File:

            file = conditionPair2File[(cond1, cond2)]

            compDF = DataFrame.parseFromFile(file)
            # 0 evidence, 1 gene, 2 cond1, 3 cond2, 4+ methods

            compHeader = compDF.getHeader()

            methods = []
            for i in range(4, len(compHeader)-1):
                colVal = compHeader[i]

                if colVal.endswith("_log2FC"):
                    method = colVal.replace('_log2FC', '')
                    methods.append(method)
                    allMethods.add(method)

            (headHTML, bodyHTML) = compDF.export(outFile=None, exType=ExportTYPE.HTML_STRING)

            if 'NOISeq' in methods:
                bodyHTML = bodyHTML + "<h3>Warning: the NOISeq calculated p-values are transformed z-scores from probabilities of differential expression</h3>"

            pltcfg = PlotConfig()
            pltcfg.setOutputType(PlotSaveTYPE.HTML_STRING)

            pltcfg.d3js = os.path.relpath(d3js_dest, outputFolder)
            pltcfg.mpld3js = os.path.relpath(mpld3_dest, outputFolder)

            geneNames = compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex('id')).to_list()

            l2FCs = []
            l2fcSuffix = "_log2FC"

            def parseList(x):
                ret = [None] * len(x)
                for i in range(0, len(x)):
                    if x[i] != 'None' and x[i] != 'NA':
                        ret[i] = float(x[i])

                return ret

            for method in methods:
                l2FCTitle = method + l2fcSuffix
                rawpTitle = method + "_RAW.PVAL"
                adjpTitle = method + "_ADJ.PVAL"

                l2FCs.append(l2FCTitle)

                l2FCdata = parseList(compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(l2FCTitle)).to_list())
                rawPdata = parseList(compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(rawpTitle)).to_list())
                adjPdata = parseList(compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(adjpTitle)).to_list())

                PorePlot.volcanoPlot(geneNames, l2FCdata, rawPdata, "Volcano Plot " + cond1 + " vs " + cond2 + "\n ("+method+")", "log2 FC", "raw pValue", pltcfg)

                if method in ['NOISeq']:
                    pltcfg.addHTMLPlot("<h3>Warning: the calculated p-values are transformed z-scores from probabilities of differential expression</h3>")

                PorePlot.volcanoPlot(geneNames, l2FCdata, adjPdata, "Volcano Plot " + cond1 + " vs " + cond2 + "\n ("+method+")", "log2 FC", "adj pValue", pltcfg)

                if method in ['NOISeq']:
                    pltcfg.addHTMLPlot("<h3>Warning: the calculated p-values are transformed z-scores from probabilities of differential expression</h3>")




            for methodI in range(0,len(methods)):
                for methodJ in range(methodI+1,len(methods)):


                    il2FCTitle = methods[methodI] + l2fcSuffix
                    jl2FCTitle = methods[methodJ] + l2fcSuffix

                    il2FCdata = compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(il2FCTitle)).to_dict()
                    jl2FCdata = compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(jl2FCTitle)).to_dict()

                    idata = []
                    jdata = []

                    unionIdx = set([x for x in jl2FCdata] + [x for x in il2FCdata])

                    for val in unionIdx:

                        if val in il2FCdata and val in jl2FCdata:

                            ielem = il2FCdata[val]
                            jelem = jl2FCdata[val]

                            if ielem == None or ielem == "None" or jelem == None or jelem=="None":
                                continue

                            idata.append(float(ielem))
                            jdata.append(float(jelem))


                    print(method, "i/jdata", len(idata), len(jdata))

                    PorePlot.plotscatter(idata, jdata, "log2FC comparison " + methods[methodI] + " vs " + methods[methodJ], "log2 FC ("+methods[methodI] + ")", "log2 FC ("+methods[methodJ] + ")", None, pltcfg)

            if all([x in methods for x in ['NOISeq', 'msEmpiRe', 'DESeq2']]):

                method2calls = {}#defaultdict(set)
                for method in ['NOISeq', 'msEmpiRe', 'DESeq2']:
                    adjpTitle = method + "_ADJ.PVAL"
                    adjPdata = compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(adjpTitle)).to_pairs()
                    method2calls[method] = set()
                    for geneid, pval in adjPdata:

                        fpval = toFloat(pval)

                        if fpval != None and fpval < 0.05:

                            #if method == 'msEmpiRe' and not geneid in method2calls["DESeq"] and not geneid in method2calls["NOISeq"]:
                            #    print(geneid, pval, fpval)

                            method2calls[method].add(geneid)


                def compute_venn3_sets(a,b,c):

                    def set_size(x):
                        return x

                    if not (type(a) == type(b) == type(c)):
                        raise ValueError("All arguments must be of the same type")

                    return (set_size(a - (b | c)),  # TODO: This is certainly not the most efficient way to compute.
                            set_size(b - (a | c)),
                            set_size((a & b) - c),
                            set_size(c - (a | b)),
                            set_size((a & c) - b),
                            set_size((b & c) - a),
                            set_size(a & b & c))


                allsets = []
                alllabels = []

                for method in sorted([x for x in method2calls]):
                    allsets.append(method2calls[method])
                    alllabels.append(method2calls[method])

                allss = compute_venn3_sets(*allsets)
                print(alllabels)
                print("c-(b|a)", len(allss[3]))

                PorePlot.plotVennDiagram(method2calls, title="Overlap of DE methods for adj. pval 0.05", pltcfg=pltcfg)

            outFilePath = base + cond1.replace("/", "_").replace(".", "_") + "_" + cond2.replace("/", "_").replace(".", "_") + ".html"
            with open( outFilePath, 'w') as resultHTML:

                mpld3js = "<script src=" + pltcfg.mpld3js + "></script>\n"
                d3js = "<script src=" + pltcfg.d3js + "></script>\n"

                outStr = "<html><head>" + headHTML +"\n" + d3js + mpld3js + "</head><body><p>" + "\n".join(pltcfg.getCreatedPlots()) + "</p>" + bodyHTML + "</body></html>"
                resultHTML.write( outStr )


                cond2genecount = {}
                for cond in [cond1, cond2]:
                    condSamples = replicates[cond]
                    gene2count = defaultdict(lambda: 0.0)

                    for sample in condSamples:

                        dr = compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(sample))

                        genes = dr.getHeader()
                        print(genes[:10])

                        for gene in genes:
                            gene2count[gene] += dr[gene]

                    if len(condSamples)>0:
                        for gene in gene2count:
                            gene2count[gene] = gene2count[gene] / len(condSamples)

                    cond2genecount[cond] = DataRow.fromDict(gene2count)

                simScore = SimilarityAnalysis.calcSimilarity(cond2genecount[cond1], cond2genecount[cond2])

                retStuff = {
                    "path": outFilePath,
                    "sim_score": simScore
                }
                for method in methods:

                    methodL2FC = method + l2fcSuffix

                    if not methodL2FC in l2FCs:
                        continue

                    l2FCCol = compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(methodL2FC))
                    avgL2FC = self.calcAvgL2FC(l2FCCol)

                    retStuff[method] = avgL2FC

                fileOutPuts[(cond1, cond2)] = retStuff


        # making overview!
        overviewDF = DataFrame()

        overviewCols = ['Condition1', 'Condition2', 'Similarity', 'Comparison HTML'] + ['Similarity (avg log2 FC) ' + method for method in allMethods]

        overviewDF.addColumns( overviewCols )

        allMethods = sorted(allMethods)
        avgL2FCByMethod = defaultdict(lambda : defaultdict())

        for condPair in fileOutPuts:

            condData = fileOutPuts[condPair]
            relHTMLFile = os.path.relpath(condData["path"], basePath)


            condPairData = {
                'Condition1': condPair[0],
                'Condition2': condPair[1],
                'Similarity': condData["sim_score"],
                'Comparison HTML': "<a target='_blank' href='"+relHTMLFile+"'>Comparison</a>"
            }

            for midx, method in enumerate(allMethods):
                condPairData['Similarity (avg log2 FC) ' + method] = condData.get(method, 0)

            overviewDF.addRow( DataRow.fromDict(condPairData) )


        (headHTML, bodyHTML) = overviewDF.export(outFile=None, exType=ExportTYPE.HTML_STRING)

        outFilePath = base + "overview.html"
        with open(outFilePath, 'w') as resultHTML:
            outStr = "<html><head>" + headHTML + "</head><body>" + bodyHTML + "</body></html>"
            resultHTML.write(outStr)

        return fileOutPuts

    def calcAvgL2FC(self, datarow):

        cnt = 0
        sum = 0

        for x in datarow.to_list():

            if x == 'None' or x == 'NA':
                xn = None
            else:
                xn = toNumber(x)

            if xn == None:
                #print("Not a number in calcAvgL2FC", str(xn))
                pass


            if xn != None:
                sum += abs(xn)
                cnt += 1

        return sum / cnt

    def getConditions(self):

        ownHeader = self.getHeader()
        return [ownHeader[2], ownHeader[3]]