import argparse
from collections import defaultdict
from typing import OrderedDict

import HTSeq
import matplotlib
from porestat.utils.ArgParseExt import FolderType

from porestat.plots.poreplot import PorePlot, MultiAxesPointHTMLTooltip

from porestat.plots.plotconfig import PlotConfig, PlotSaveTYPE
from porestat.utils.EnrichmentDF import EnrichmentDF
from ..utils.DataFrame import DataFrame, DataRow, ExportTYPE
from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory, PSToolException, PSToolInterface

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Files import fileExists, eprint
from scipy import stats
import sys, math
import numpy as np
from matplotlib import pyplot as plt
import random
import os

from scipy.optimize import curve_fit


class FoldChangeFeatureCountsDistributionFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):
        super(FoldChangeFeatureCountsDistributionFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which + ' help')
        parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), default=None,
                            help='counts summary file', required=False)
        parser.add_argument('-p', '--prefixes', nargs='+', type=str, default=None,
                            help='counts summary file', required=False)

        parser.add_argument('-co', '--conditions', nargs='+', type=str, action='append', default=None,
                            help='counts summary file', required=False)

        parser.add_argument('-v', '--no-analysis', dest='noanalysis', action='store_true', default=False)
        parser.add_argument('-m', '--methods', type=str, nargs='+', default=['DESeq2', 'DirectDESeq2', 'edgeR', 'limma', 'msEmpiRe', "nlEmpiRe"])

        parser.add_argument('-o', '--output', type=FolderType('w'), help='output location, default: std out',
                            required=True)
        parser.add_argument('-r', '--rscript', type=argparse.FileType('r'), help='path to Rscript',
                            default='/usr/bin/Rscript')

        parser.add_argument('-e', '--enhanced', type=argparse.FileType('r'), default=None)
        parser.add_argument('-l', '--lengths', type=argparse.FileType('r'), default=None)

        parser.add_argument('-rrna', '--no-rrna', dest='norrna', action='store_true', default=False)
        parser.add_argument('-rmtrna', '--remove-mtrna', dest='removemtrna', action='store_true', default=False)
        parser.add_argument('-opc', '--only-protein-coding', dest='only_protein_coding', action='store_true', default=False)

        parser.add_argument('-fpkm', '--fpkm', dest='fpkm', action='store_true', default=False)
        parser.add_argument('-tpm', '--tpm', dest='tpm', action='store_true', default=False)
        parser.add_argument('-libsize', '--libsize', dest='libsize', action='store_true', default=False)


        parser.add_argument('-anc', '--allow-nonexistant-cond', dest='allow_nonexistant_cond', action='store_true', default=False)

    

        parser = PlotConfig.addParserArgs(parser)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):
        if (not PSToolInterface.hasArgument('counts', args) and not PSToolInterface.hasArgument('diffreg', args)) or (
                args.counts == None and args.diffreg == None):
            raise argparse.ArgumentParser().error(
                "error: Either counts [--counts] or diffreg results [--diffreg] must be set!")

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return FoldChangeFeatureCountsAnalysis(simArgs)


class FoldChangeFeatureCountsAnalysis(ParallelPSTInterface):

    def __init__(self, args):

        super(FoldChangeFeatureCountsAnalysis, self).__init__(args)

        self.counts = None
        self.condData = EnrichmentDF()

    def _makePropDict(self):

        return None

    def readCounts(self, args, biotypes=None, gene2length=None):

        if args.norrna and biotypes == None:
            raise argparse.ArgumentParser().error("removal of rRNA requires --enhanced!")
        if args.removemtrna and biotypes == None:
            raise argparse.ArgumentParser().error("removal of mtRNA requires --enhanced!")
        if args.only_protein_coding and biotypes == None:
            raise argparse.ArgumentParser().error("--only-protein-coding requires --enhanced!")

        if args.fpkm and gene2length == None:
            raise argparse.ArgumentParser().error("calculation of FPKM requires --lengths!")
        if args.tpm and gene2length == None:
            raise argparse.ArgumentParser().error("calculation of TPM requires --lengths!")

        featureCountsColumns = ["Geneid",	"Chr",	"Start",	"End",	"Strand",	"Length"]

        counts = defaultdict(lambda: list())

        condition2samples = defaultdict(list)

        for idx, countFile in enumerate(args.counts):

            countFilePrefix = args.prefixes[idx]
            df = DataFrame.parseFromFile(countFile.name, skipChar='#')

            allheaders = df.getHeader()
            sampleHeaders = [x for x in allheaders if not x in featureCountsColumns]

            for sample in sampleHeaders:
                condition2samples[sample].append(countFilePrefix + sample)


            for condGroup in args.conditions:
                condName = condGroup[0]
                for condElement in condGroup:

                    print(condName, condElement, condElement in allheaders)
                    if args.allow_nonexistant_cond and not condElement in allheaders:
                        continue

                    subDf = df.selectColumns({"Geneid": "gene", condElement: "count"})

                    if biotypes != None and args.norrna:
                        geneColIdx = subDf.getColumnIndex("gene")
                        subDf.filterRows(lambda x: x[geneColIdx] in biotypes and not "rRNA" in biotypes[x[geneColIdx]][1] )

                    if biotypes != None and args.removemtrna:
                        geneColIdx = subDf.getColumnIndex("gene")
                        subDf.filterRows(lambda x: x[geneColIdx] in biotypes and not "Mt_" in biotypes[x[geneColIdx]][1] )

                    if biotypes != None and args.only_protein_coding:
                        subDf.filterRows(lambda x: x[geneColIdx] in biotypes and "protein_coding" in biotypes[x[geneColIdx]][1] )

                    if os.path.isdir(condElement) or os.path.isfile(condElement):
                        subDf.setFilepath(os.path.abspath(condElement))
                    else:
                        subDf.setFilepath(condElement)

                    if args.libsize:
                        countCol = subDf.getColumnIndex("count")
                        geneCol = subDf.getColumnIndex("gene")

                        totalCounts = sum([x[countCol] for x in subDf.data])

                        libSizeIdx = subDf.addColumn("LS", 0)

                        def addLibSize(x):
                            x[libSizeIdx] = (x[countCol]/totalCounts) * 10000

                            return tuple(x)

                        subDf.applyToRow(addLibSize)     


                    if args.fpkm:

                        countCol = subDf.getColumnIndex("count")
                        geneCol = subDf.getColumnIndex("gene")

                        totalCounts = sum([x[countCol] for x in subDf.data])

                        fpkmIdx = subDf.addColumn("FPKM", 0)

                        def addFPKM(x):

                            geneID = x[geneCol]
                            geneLength = gene2length.get(geneID, 0)

                            if geneLength == 0:
                                x[fpkmIdx] = 0
                            else:
                                x[fpkmIdx] = x[countCol]/(totalCounts*geneLength) * pow(10,9)

                            return tuple(x)

                        subDf.applyToRow(addFPKM)

                    if args.tpm:

                        countCol = subDf.getColumnIndex("count")
                        geneCol = subDf.getColumnIndex("gene")

                        totalCounts = sum([x[countCol] for x in subDf.data])
                        totalRatio = 0

                        for row in subDf:
                            geneID = row["gene"]
                            geneCount = row["count"]
                            geneLength = gene2length.get(geneID, 0)

                            if geneLength == 0:
                                pass
                            else:
                                totalRatio += geneCount/geneLength

                        tpmIdx = subDf.addColumn("TPM", 0)

                        def addTPM(x):

                            geneID = x[geneCol]
                            geneLength = gene2length.get(geneID, 0)

                            if geneLength == 0:
                                x[fpkmIdx] = 0
                            else:
                                x[tpmIdx] = x[countCol]/(geneLength * totalRatio) * pow(10,6)

                            return tuple(x)

                        subDf.applyToRow(addTPM)

                    counts[condName].append(subDf)

        return counts, condition2samples

    def readDiffRegs(self, args):

        # TODO what did I want to do with this argument?
        for diffFile in args.diffreg:
            df = EnrichmentDF.parseFromFile(diffFile)

    def prepareInputs(self, args):
        return []

    def execParallel(self, data, environment):

        return None

    def joinParallel(self, existResult, newResult, oEnvironment):

        return None

    def loadEnhancement(self, fileE):

        if fileE == None:
            print("Not loading gene name enhancements")
            return {}

        print("Loading gene name enhancements", fileE.name)

        ens2sym = {}

        for lidx, line in enumerate(fileE):
            line = line.strip().split("\t")

            if lidx == 0 or line[0].startswith("#"):
                continue

            ensemblID = line[0]
            geneSymbol = line[1]
            biotype = line[2]

            #if len(geneSymbol) == 0:
            #    continue

            ens2sym[ensemblID] = (geneSymbol, biotype)

        return ens2sym

    def loadGeneLengths(self, fileE):

        if fileE == None:
            print("Not loading gene lengths")
            return None

        print("Loading gene lengths", fileE.name)

        """
            Ensembl_gene_identifier GeneID  length
            ENSMUSG00000000001      14679   3262
            ENSMUSG00000000003      54192   902
            ENSMUSG00000000028      12544   2252
        """

        ens2gl = {}
        for lidx, line in enumerate(fileE):
            line = line.strip().split("\t")

            if lidx == 0:
                try:
                    int(line[1])
                except:
                    continue

            ensemblID = line[0]
            geneLength = line[1]

            if len(ensemblID) == 0 or len(geneLength) == 0:
                continue

            geneLength = int(geneLength)
            ens2gl[ensemblID] = geneLength

        return ens2gl


    def makeResults(self, parallelResult, oEnvironment, args):

        if not args.counts == None:

            """
            counts is a defaultdict(list) for each condition name with maybe multiple samples
            """

            geneEnhancement = self.loadEnhancement(args.enhanced)
            geneLengths = self.loadGeneLengths(args.lengths)

            counts, cond2samples = self.readCounts(args, biotypes=geneEnhancement, gene2length=geneLengths)

            #vConds = sorted([x for x in counts])
            vConds = [x for x in counts]

            createdComparisons = defaultdict(list)
            conditions = []

            for valueSource in ['count']:
                self.condData = EnrichmentDF()
                replicates = OrderedDict()

                for condition in vConds:

                    condData = counts[condition]

                    condReplicates = []
                    for condDataSample in condData:
                        geneNames = condDataSample.getColumnIndex('gene')
                        geneCounts = condDataSample.getColumnIndex(valueSource)

                        rowUpdates = []
                        sampleName = condDataSample.filepath

                        print(sampleName, len(condDataSample))
                        for row in condDataSample:

                            rowData = {
                                    "id": row["gene"],
                                    sampleName: row[valueSource]
                                }

                            if args.libsize:
                                rowData[sampleName+".LS"] = row["LS"]

                            if args.fpkm:
                                rowData[sampleName+".FPKM"] = row["FPKM"]

                            if args.tpm:
                                rowData[sampleName+".TPM"] = row["TPM"]

                            rowUpdates.append(rowData)

                        #condRows = condDataSample.namedRows(geneNames, interestCols)
                        #condRow = condDataSample.toDataRow(geneNames, geneCounts)

                        conditions.append(sampleName)
                        condReplicates.append(sampleName)

                        print("Add Condition", sampleName, rowUpdates[0])
                        self.condData.addConditions(rowUpdates, sampleName)

                    replicates[condition] = condReplicates

                print("Running for conditions: " + str(vConds))

                createdComparisons[valueSource] += self.condData.runDEanalysis(args.output, prefix=valueSource,
                                                                               rscriptPath=args.rscript.name,
                                                                               methods=args.methods,
                                                                               replicates=replicates,
                                                                               noDErun=args.noanalysis,
                                                                               enhanceSymbol=geneEnhancement,
                                                                               geneLengths=geneLengths
                                                                               )

            self.prepareHTMLOut(createdComparisons, replicates, args)


    def getValueSource(self, df):
        return df.data[0][1]

    def prepareHTMLOut(self, createdComparisons, replicates, args):

        for valueSource in createdComparisons:

            allComparisons = createdComparisons[valueSource]

            condPair2File = {}
            for x in allComparisons:
                condPair2File[(x[0], x[1])] = x[2]

            print("Comparisons")
            for x in condPair2File:
                print(x, condPair2File[x])

            self.condData.printResult(args.output, prefix=valueSource, conditionPair2File=condPair2File,
                                      replicates=replicates)  # conditions=conditions, files=createdComparisons[valueSource])
