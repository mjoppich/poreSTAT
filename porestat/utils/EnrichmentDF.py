from ..analysis.similarity_analysis import SimilarityAnalysis
from .DataFrame import DataFrame, DataRow, ExportTYPE
import os
import mpld3
from ..plots.plotconfig import PlotConfig, PlotSaveTYPE
from ..plots.poreplot import PorePlot
import pickle


class EnrichmentDF(DataFrame):

    def __init__(self, otherDF=None):
        super(EnrichmentDF, self).__init__()
        self.idCol = self.addColumn('id')

        if otherDF:
            self.data = otherDF.data
            self.column2idx = otherDF.column2idx

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

        for i in range(0, len(self.data)):
            row = self.data[i]
            rowID = row[ self.idCol ]

            if rowID in condData.getHeader():
                lrow = list(row)
                lrow[newColIdx] = condData[rowID]
                self.data[i] = tuple(lrow)
                setFoundKeys.add(rowID)

        for geneID in condData.getHeader():
            if not geneID in setFoundKeys:

                lrow = [None] * len(self.column2idx)
                lrow[0] = geneID
                lrow[newColIdx] = condData[geneID]

                self.data.append(tuple(lrow))


    def writeEnrichmentBrowserFiles(self, prepData, exprFile, pdataFile, fdataFile):

        self.writeExpressionFile(exprFile, prepData)
        self.writepDataFile(pdataFile, prepData)
        self.writefDataFile(fdataFile, prepData)



    def prepareEBData(self, cond1, cond2):

        header = ['gene', cond1, cond2]
        prepData = [tuple(header)]

        cond1Idx = self.getColumnIndex(cond1)
        cond2Idx = self.getColumnIndex(cond2)

        for row in self:

            count1 = int(row[cond1Idx])
            count2 = int(row[cond2Idx])

            gene = row[0]

            if count1 == None or count2 == None:
                continue

            prepData.append( (gene, count1, count2) )

        return prepData


    def writepDataFile(self, file, prepData):

        with open(file, 'w') as pdataFile:

            conds = prepData[0][1:]
            cnt = 0
            for cond in conds:
                pdataFile.write(cond + "\t" + str(cnt) + "\n")
                cnt += 1


    def writefDataFile(self, file, prepData):

        with open(file, 'w') as fdataFile:

            for line in prepData[1:]:
                fdataFile.write( line[0] + "\t" + line[0] + "\n")

    def writeExpressionFile(self, file, prepData):

        with open(file, 'w') as exprFile:

            for line in prepData[1:]:
                exprFile.write( "\t".join([str(x) for x in line[1:]]) + "\n")



    def runDEanalysis(self, outputFolder, rscriptPath="/usr/bin/Rscript", prefix= "", conditions=None, methods=['DESeq', 'edgeR']):

        if conditions == None:
            conditions = self.getHeader()[1:]

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

        scriptPath = os.path.dirname(os.path.abspath(__file__)) + "/../data/de_rseq.R"

        createdComparisons = []

        for i in range(0, len(conditions)):
            cond1 = conditions[i]
            for j in range(i+1, len(conditions)):
                cond2 = conditions[j]

                prepData = self.prepareEBData(cond1, cond2)
                self.writeEnrichmentBrowserFiles(prepData, exprFile, pdataFile, fdataFile)

                condResult = {}

                for method in methods: # 'limma' 'edgeR'
                    outFile = outFileBase + "_" + method

                    execStr = rscriptPath+" "+scriptPath+" "+exprFile+" "+pdataFile+" "+fdataFile+" "+method+" " + outFile
                    print(execStr)

                    os.system(execStr)

                    methDF = DataFrame.parseFromFile(outFile)
                    condResult[method] = methDF

                geneNames = self.getColumnIndex('id')
                cond1Counts = self.toDataRow(geneNames, self.getColumnIndex(cond1))
                cond2Counts = self.toDataRow(geneNames, self.getColumnIndex(cond2))


                compDF = EnrichmentDF()
                evCol = compDF.addColumn('evidence', prefix) # TODO this is just a quick hack!
                compDF.addCondition(cond1Counts, cond1)
                compDF.addCondition(cond2Counts, cond2)

                condVPData = {}
                usedMethod = []
                for method in condResult:

                    methDF = condResult[method]

                    if methDF == None:
                        continue

                    usedMethod.append(method)

                    l2FCTitle = method + "_log2FC"
                    rawpTitle = method + "_RAW.PVAL"
                    adjpTitle = method + "_ADJ.PVAL"

                    l2FCdata = methDF.toDataRow( methDF.getColumnIndex('GENE.ID'), methDF.getColumnIndex('log2FC') )
                    rawPdata = methDF.toDataRow(methDF.getColumnIndex('GENE.ID'), methDF.getColumnIndex('RAW.PVAL'))
                    adjPdata = methDF.toDataRow(methDF.getColumnIndex('GENE.ID'), methDF.getColumnIndex('ADJ.PVAL'))

                    compDF.addCondition(l2FCdata, l2FCTitle)
                    compDF.addCondition(rawPdata, rawpTitle)
                    compDF.addCondition(adjPdata, adjpTitle)

                    condVPData[method] = ( l2FCdata.getHeader(), l2FCdata.to_list(), rawPdata.to_list(), adjPdata.to_list() )


                addInfo = compDF.addColumn("EnsemblBacteria")

                def addInfoFunc(x):
                    gene = x[0]
                    x[evCol] = prefix

                    link = "<a target='_blank' href='http://bacteria.ensembl.org/Helicobacter_pylori_p12/Gene/Summary?g="+gene+"'>EnsemblBacteria</a><br/>" \
                           "<a target='_blank' href='http://www.uniprot.org/uniprot/?query="+gene+"&sort=score'>Uniprot</a>"

                    if not gene.startswith("HP"):
                        link = ""

                    x[addInfo] = link
                    return tuple(x)

                compDF.applyToRow( addInfoFunc )

                outname = base + cond1 + "_" + cond2 + ".tsv"

                createdComparisons.append(outname)

                compDF.export(outFile=outname, exType=ExportTYPE.TSV)

        return createdComparisons


    def printResult(self, outputFolder, prefix, conditions, files):

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

        for file in files:

            compDF = DataFrame.parseFromFile(file)
            # 0 evidence, 1 gene, 2 cond1, 3 cond2, 4+ methods

            compHeader = compDF.getHeader()

            cond1 = compHeader[2]
            cond2 = compHeader[3]

            if not (cond1 in conditions and cond2 in conditions):
                continue # comparison not needed

            methods = []
            for i in range(4, len(compHeader)-1, 3):
                colVal = compHeader[i]
                method = colVal.replace('_log2FC', '')
                methods.append(method)

            (headHTML, bodyHTML) = compDF.export(outFile=None, exType=ExportTYPE.HTML_STRING)

            pltcfg = PlotConfig()
            pltcfg.setOutputType(PlotSaveTYPE.HTML_STRING)

            pltcfg.d3js = os.path.relpath(d3js_dest, outputFolder)
            pltcfg.mpld3js = os.path.relpath(mpld3_dest, outputFolder)

            geneNames = compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex('id')).to_list()

            l2FCs = []
            l2fcSuffix = "_log2FC"

            for method in methods:
                l2FCTitle = method + l2fcSuffix
                rawpTitle = method + "_RAW.PVA"
                adjpTitle = method + "_ADJ.PVAL"

                def parseList(x):
                    ret = [None] * len(x)
                    for i in range(0, len(x)):
                        if x[i] != 'None':
                            ret[i] = float(x[i])

                    return ret

                l2FCs.append(l2FCTitle)

                l2FCdata = parseList(compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(l2FCTitle)).to_list())
                rawPdata = parseList(compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(rawpTitle)).to_list())
                adjPdata = parseList(compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(adjpTitle)).to_list())

                PorePlot.volcanoPlot(geneNames, l2FCdata, rawPdata, "Volcano Plot " + cond1 + " vs " + cond2 + "\n ("+method+")", "log2 FC", "raw pValue", pltcfg)
                PorePlot.volcanoPlot(geneNames, l2FCdata, adjPdata, "Volcano Plot " + cond1 + " vs " + cond2 + "\n ("+method+")", "log2 FC", "adj pValue", pltcfg)

            outFilePath = base + cond1 + "_" + cond2 + ".html"
            with open( outFilePath, 'w') as resultHTML:

                mpld3js = "<script src=" + pltcfg.mpld3js + "></script>\n"
                d3js = "<script src=" + pltcfg.d3js + "></script>\n"

                outStr = "<html><head>" + headHTML +"\n" + d3js + mpld3js + "</head><body><p>" + "\n".join(pltcfg.getCreatedPlots()) + "</p>" + bodyHTML + "</body></html>"
                resultHTML.write( outStr )

                fileOutPuts[(cond1, cond2)] = tuple([compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(cond1)),
                                               compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(cond2)),
                                               outFilePath] + [compDF.toDataRow(compDF.getColumnIndex('id'), compDF.getColumnIndex(lfcm)) for lfcm in l2FCs],)


        # making overview!
        overviewDF = DataFrame()
        overviewDF.addColumns(['Condition1', 'Condition2', 'Similarity', 'Comparison HTML'])

        avgL2FCByMethod = defaultdict(lambda : defaultdict)

        for condPair in fileOutPuts:

            condData = fileOutPuts[condPair]
            relHTMLFile = os.path.relpath(condData[2], basePath)

            simScore = SimilarityAnalysis.calcSimilarity(condData[0], condData[1])

            avgL2FCByMethod[methods[0]][(condPair[0], condPair[1])] = self.calcAvgL2FC(condPair[3])
            avgL2FCByMethod[methods[1]][(condPair[0], condPair[1])] = self.calcAvgL2FC(condPair[4])

            condPairData = {
                'Condition1': condPair[0],
                'Condition2': condPair[1],
                'Similarity': simScore,
                'Similarity (avg log2 FC) ' + methods[0]: avgL2FCByMethod[methods[0]][(condPair[0], condPair[1])],
                'Similarity (avg log2 FC) ' + methods[1]: avgL2FCByMethod[methods[1]][(condPair[0], condPair[1])],
                'Comparison HTML': "<a target='_blank' href='"+relHTMLFile+"'>Comparison</a>"
            }

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
            if x != None:
                sum += abs(x)
                cnt += 1

        return sum / cnt

    def getConditions(self):

        ownHeader = self.getHeader()
        return [ownHeader[2], ownHeader[3]]