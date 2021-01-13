import argparse
import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")


import pickle
import scipy

from scipy import stats
from collections import Counter, defaultdict, OrderedDict

import networkx
import requests
import json
import datetime

from statsmodels.stats.multitest import multipletests
from scipy.stats import hypergeom

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
from porestat.DEtools.mirtools.miRNAUtils import miRNA, miRNAPART, isNumber
from porestat.plots.GraphPlotter import GraphPlot


class DataBaseAccessor:
    serverAddress = "http://localhost"
    serverPort = "65500"
    serverPath = "/"

    @classmethod
    def makeServerAddress(cls, address, port, path, page):

        ret = address

        if port != None:
            ret += ":" + str(port) + "/"

        if path != None:

            if port == None:
                ret += "/"

            ret += path + "/"

        if page != None:
            if port == None and path == None:
                ret += "/"

            ret += page

        return ret

    @classmethod
    def checkContext(cls, requestDict):

        r = requests.post(cls.makeServerAddress(cls.serverAddress, cls.serverPort, cls.serverPath, "check_context"),
                          data=json.dumps(requestDict))

        jsonRes = r.json()

        for x in jsonRes:
            jsonRes[x] = [{"termid": y} for y in jsonRes[x]]

        return jsonRes

    @classmethod
    def fetchSimple(cls, requestDict):

        r = requests.post(cls.makeServerAddress(cls.serverAddress, cls.serverPort, cls.serverPath, "find_interactions"),
                          data=json.dumps(requestDict))

        jsonRes = r.json()

        return jsonRes

    @classmethod
    def fetch_mirna_interactions(cls, requestDict, MIRNASTRPARTS=[miRNAPART.MATURE, miRNAPART.ID], verbose=False):

        if not "sentences" in requestDict:
            requestDict["sentences"] = "FALSE"

        if verbose:
            print(json.dumps(requestDict))

        jsonRes = cls.fetchSimple(requestDict)

        if not "rels" in jsonRes:
            print(jsonRes)

        foundInteractions = set()

        for rel in jsonRes['rels']:

            source = rel['lid']
            target = rel['rid']

            try:
                target = miRNA(target)
                target = target.getStringFromParts(MIRNASTRPARTS, normalized=True)

            except:
                pass

            hasCanonicalRegEvidence = False

            for ev in rel["evidences"]:

                if ev["data_source"] == "pmid":

                    if ev["rel_interaction"] == "MIR_GENE" and ev["rel_category"] in ["DOWN", "NEU"]:
                        hasCanonicalRegEvidence = True

                else:

                    if not ev.get("data_source", "") in ["DIANA"]:
                        hasCanonicalRegEvidence = True


            if hasCanonicalRegEvidence:
                edge = (source, target)
                foundInteractions.add(edge)

        return foundInteractions


class miRGeneGraph:

    def logFC2Size(self, logFC, minS=20, maxS=40):

        logFCr = maxLogFC + (-minLogFC)

        logFCx = logFC + -minLogFC

        logFCx = logFCx / logFCr

        return minS + logFCx * (maxS - minS)

    def exprColor(self, logfc, nodeDirection):

        if not str(logfc).isnumeric():

            if logfc == "up":
                logfc = 1
            elif logfc == "down":
                logfc = -1

        if logfc == "down" or logfc < 0:
            return "#FF0000"

        if logfc == 'up' or logfc > 0:
            return "#00FF00"

        if logfc == 0 and nodeDirection.upper() in ["UNCHANGED"]:
            return "#c7c7c7"

        return "#555555"

    def findDataForNode(self, nodeName, isDefinitelyGene):
        un = nodeName.upper()

        if not isDefinitelyGene and ("MIR" in un or "LET" in un):
            nodeData = None

            try:
                objMir = miRNA(nodeName)

                for mir in allMIRs:
                    if objMir.accept(mir):
                        nodeData = allMIRs.get(mir, None)
                        #print("Matching", nodeName, "with", mir)
                        break

                nodeType = "mirna"
                nodeShape = "triangle"

            except:
                print("Could not parse", nodeName, "as miRNA", file=sys.stderr)
                return None, False

        else:
            nodeData = allGenes.get(un, None)
            nodeType = "gene"
            nodeShape = "square"

        if nodeData == None:
            return {"type": nodeType,
                    "shape": nodeShape,
                    "border_style": "dashed",
                    "log2FC": 0,
                    "adjPval": 1,
                    "de_measured": "false",
                    "node_expr_detection": "imputed0",
                    "node_expr_direction": "N/A"
                    }, False

        isMeasured = True if (nodeData[0], nodeData[1]) != (0.0, 1.0) else False

        returnDict = {
            "de_measured": str(isMeasured).lower(),
            "border_style": "dashed" if nodeData[1] > sigThreshold else "solid",
            "type": nodeType,
            "log2FC": nodeData[0],
            "adjPval": nodeData[1],
            "shape": nodeShape,
            "size": self.logFC2Size(nodeData[0]),
            "node_expr_detection": "data" if isMeasured else "imputed0",
        }

        if nodeType == "mirna":
            dataOrigin = ";".join(nodeData[2])
            returnDict["mirna_data_origin"] = dataOrigin

        if isMeasured:

            if nodeData[0] < 0:
                returnDict["node_expr_direction"] = "down"
            elif nodeData[0] > 0:
                returnDict["node_expr_direction"] = "up"
            else:
                returnDict["node_expr_direction"] = "unchanged"
        else:
            returnDict["node_expr_direction"] = "N/A"

        returnDict["color"] = self.exprColor(nodeData[0], returnDict["node_expr_direction"])

        return returnDict, True

    def __init__(self, props):
        self.priorityEdgeTypes = ["expected", "unexpected"]
        self.dirToOpposite = {"up": "down", "down": "up", "N/A": "N/A"}

        self.gp = GraphPlot()

        self.minLogFC = props["minLogFC"] if "minLogFC" in props else -5
        self.maxLogFC = props["maxLogFC"] if "maxLogFC" in props else 5

    def createGraph(self, mirnaHits, geneHits, genename2mirs, graph=None):

        if graph is None:
            graph = networkx.Graph()

        deNodes = set()

        print("Starting mirnaHits")

        if not mirnaHits is None:
            for edge in mirnaHits:

                if edge[0].upper() in genename2mirs and edge[1].upper() in genename2mirs:
                    print("skip edge", edge, "for genename2mirs")
                    continue

                srcData, srcDE = self.findDataForNode(edge[0], edge[0] in genename2mirs)
                tgtData, tgtDE = self.findDataForNode(edge[1], edge[1] in genename2mirs)

                if srcData != None and tgtData != None:

                    if srcDE:
                        deNodes.add(edge[0])
                    if tgtDE:
                        deNodes.add(edge[1])

                    graph.add_node(edge[0], attr_dict=srcData)
                    graph.add_node(edge[1], attr_dict=tgtData)
                    graph.add_edge(edge[0], edge[1])

                    graph.edges[edge]['edge_creation'] = "original_mirhit"

        print("Starting geneHits")

        if not geneHits is None:
            for edge in geneHits:

                if edge[0].upper() in genename2mirs and edge[1].upper() in genename2mirs:
                    print("skip edge", edge, "for genename2mirs")
                    continue

                srcData, srcDE = self.findDataForNode(edge[0], edge[0] in genename2mirs)
                tgtData, tgtDE = self.findDataForNode(edge[1], edge[1] in genename2mirs) # was edge[1], edge[0]

                if srcData != None and tgtData != None:

                    if srcDE:
                        deNodes.add(edge[0])
                    if tgtDE:
                        deNodes.add(edge[1])

                    graph.add_node(edge[0], attr_dict=srcData)
                    graph.add_node(edge[1], attr_dict=tgtData)

                    graph.add_edge(edge[0], edge[1])

                    graph.edges[edge]['edge_creation'] = "original_genehit"

        print("Finished geneHits")

        # number of genes for each mirna
        for node in graph.nodes():
            nodeData = graph.node[node]["attr_dict"]

            if nodeData["type"] != "mirna":
                continue

            miRNeighbours = networkx.all_neighbors(graph, node)
            nonMirNeighbours = set()

            for y in miRNeighbours:
                neighData = graph.node[y]["attr_dict"]

                if neighData["type"] != "mirna":
                    nonMirNeighbours.add(y)

            graph.node[node]["attr_dict"]["all_targets"] = nonMirNeighbours



        # let's remove genes which have not been measured at all ...
        delNode = set()
        for node in graph.nodes():

            nodeData = graph.node[node]["attr_dict"]

            if nodeData["type"] != "gene":
                continue

            if nodeData["de_measured"] == "false":
                delNode.add(node)

            if nodeData.get("adjPval", 1.0) > 0.05:
                delNode.add(node)

        print("Removing Genes:", delNode)

        for x in delNode:
            graph.remove_node(x)

        return graph

    def imputeGraph(self, graph):

        self.removeSingletons(graph)

        self.colorInitial(graph)
        self.__impute1(graph)
        self.__impute1_1(graph)
        self.__impute2(graph)
        self.__impute3(graph)
        self.__impute4(graph)

        self.colorEdges(graph)
        self.makeEdgeCounters(graph)

        return graph

    def removeSingletons(self, graph):

        removeNodes = set()

        for node in graph.nodes():
            nodeData = graph.node[node]["attr_dict"]

            nodeNN = [x for x in graph.neighbors(node)]

            if len(nodeNN) == 0:

                removeNodes.add(node)

        for x in removeNodes:
            graph.remove_node(x)
            print("Removed from graph for singleton", x)


    def makeEdgeCounters(self, graph):

        for node in graph.nodes():

            nodeData = graph.node[node]["attr_dict"]

            nodeNN = [x for x in graph.neighbors(node)]

            correctDirection = 0
            incorrectDirection = 0
            unknownDirection = 0

            directionCounter = Counter()

            for nn in nodeNN:

                edge = (node, nn)
                src = node
                tgt = nn

                srcData = graph.node[src]["attr_dict"]
                tgtData = graph.node[tgt]["attr_dict"]

                srcDir = srcData.get("node_expr_direction", "N/A")
                tgtDir = tgtData.get("node_expr_direction", "N/A")

                directionCounter[tgtDir] += 1

                srcDataOrigin = srcData.get("node_expr_detection", "N/A")
                tgtDataOrigin = tgtData.get("node_expr_detection", "N/A")

                if srcDir == self.dirToOpposite[tgtDir] and srcDir != "N/A":
                    correctDirection += 1

                elif srcDir == tgtDir and srcDir != "N/A":

                    incorrectDirection += 1

                else:
                    unknownDirection += 1

            graph.node[node]['attr_dict']["edge_elements"] = {
                'correct': correctDirection,
                'incorrect': incorrectDirection,
                'unknown': unknownDirection
            }
            graph.node[node]['attr_dict']["node_elements"] = dict(directionCounter)

            if nodeData["node_expr_direction"] in ["N/A"]:
                print(node, nodeData["node_expr_direction"], nodeData["node_expr_detection"], directionCounter)

        return graph




    def to_print_format(self, dictElem):

        outlistH = [("Category", "Count")]
        outlist = []

        if dictElem == None:
            return outlist

        for x in dictElem:

            keyElem = x

            if type(keyElem) in [list, tuple]:
                keyElem = ", ".join([str(y) for y in x])

            outlist.append((keyElem, str(dictElem[x])))

        return outlistH + sorted(outlist)

    def getUnexplainedGenes(self, graph):

        returnResult = [("Gene", "Number of inconsistent Edges", "Number of Edges for Gene")]

        totalNodeCount = 0
        nodeConditionCount = 0
        for node in sorted([x for x in graph.nodes()]):

            nodeData = graph.node[node]["attr_dict"]

            if nodeData["type"] == "mirna":
                continue

            totalNodeCount += 1
            nodeNN = [x for x in graph.neighbors(node)]

            mirDirection = nodeData.get("node_expr_direction", "N/A")

            inconsistCount = 0
            for nn in nodeNN:
                nnData = graph.node[nn]["attr_dict"]

                nnDirection = nnData["node_expr_direction"]
                edgeType = graph.edges[(node, nn)]["edge_type"]
                edgeCreation = graph.edges[(node, nn)]["edge_creation"]

                if mirDirection == nnDirection and mirDirection != None:
                    inconsistCount += 1

            if inconsistCount == len(nodeNN):
                nodeConditionCount += 1
                print(node, inconsistCount, len(nodeNN))

                returnResult.append((node, inconsistCount, len(nodeNN)))

        print(nodeConditionCount, totalNodeCount)

        return returnResult

    def getMeasuredInconsistencies(self, graph):

        returnResult = [("Source", "Target")]

        for edge in graph.edges():

            if graph.edges[edge]["edge_type"] == "unexpected":
                returnResult.append((edge[0], edge[1]))

        return returnResult

    def getImputedInconsistencies(self, graph):

        returnResult = [("Source", "Target", "Has Other Explanation")]

        for edge in graph.edges():

            if graph.edges[edge]["dir_type"] == "unexpected_imputed":

                hasOtherExplain = False

                src = edge[0]
                tgt = edge[1]

                srcData = graph.node[src]["attr_dict"]
                tgtData = graph.node[tgt]["attr_dict"]

                if srcData["type"] == "gene":

                    srcNN = [x for x in graph.neighbors(src)]

                    for nn in srcNN:
                        if graph.edges[(src, nn)]["dir_type"] in ["expected", "expected_imputed"]:
                            hasOtherExplain = True
                            break
                elif tgtData["type"] == "gene":
                    tgtNN = [x for x in graph.neighbors(tgt)]

                    for nn in tgtNN:
                        if graph.edges[(tgt, nn)]["dir_type"] in ["expected", "expected_imputed"]:
                            hasOtherExplain = True
                            break

                if not hasOtherExplain:
                    returnResult.append((edge[0], edge[1], hasOtherExplain))

        return returnResult


    def scoreGenes(self, graph):

        for node in graph.nodes():

            graph.nodes[node]["attr_dict"]["rs"] = 0.0

            nodeData = graph.nodes[node]["attr_dict"]

            if nodeData["type"] == "mirna":
                continue

            nodeDirection = nodeData.get("node_expr_direction", "N/A")

            nodeNN = [x for x in graph.neighbors(node)] #nodeNN is miRNAs

            if nodeDirection in ["N/A", "UNKNOWN"]:
                continue


            unknownCount = 0
            consistentCount = 0
            inconsistentCount = 0

            for nn in nodeNN:
                tgt = nn
                tgtData = graph.node[tgt]["attr_dict"]
                tgtDir = tgtData.get("node_expr_direction", "N/A")

                if tgtDir == "N/A":
                    unknownCount+= 1
                    continue

                if tgtDir == self.dirToOpposite[nodeDirection]:
                    consistentCount += 1
                else:
                    inconsistentCount += 1


            effectSize = abs(nodeData.get("log2FC", 0))
            rs = 0

            if consistentCount != 0:
                rs = effectSize / consistentCount

            graph.nodes[node]["attr_dict"]["rs"] = rs

    def scoreMIRs(self, graph):

        for node in graph.nodes():

            graph.nodes[node]["attr_dict"]["ns"] = 0.0

            nodeData = graph.nodes[node]["attr_dict"]

            if not nodeData["type"] == "mirna":
                continue

            nodeDirection = nodeData.get("node_expr_direction", "N/A")

            nodeNN = [x for x in graph.neighbors(node)]

            if nodeDirection in ["N/A", "UNKNOWN"]:
                continue

            ns = abs(nodeData.get("log2FC", 0.0))
            nnns = 0

            for nn in nodeNN:
                tgt = nn
                tgtData = graph.node[tgt]["attr_dict"]
                tgtDir = tgtData.get("node_expr_direction", "N/A")

                if tgtDir == "N/A":
                    continue

                assert(tgtData.get("rs", 0.0) >= 0.0)

                if tgtDir == self.dirToOpposite[nodeDirection]:
                    nnns += tgtData.get("rs", 0.0)
                else:
                    nnns -= tgtData.get("rs", 0.0)

            ns += (nnns / len(nodeNN))
            graph.nodes[node]["attr_dict"]["ns"] = ns

    def getNSScores(self, graph):

        returnResult = [("miRNA", "Score", "z-value", "pval", "adj_pval")]


        foundScores = []

        for node in sorted([x for x in graph.nodes()]):

            nodeData = graph.node[node]["attr_dict"]

            if not nodeData["type"] == "mirna":
                continue

            foundScores.append((node, nodeData["ns"]))

        zScores = stats.zscore([x[1] for x in foundScores])

        allpvals = [scipy.stats.norm.sf(abs(x)) for x in zScores]
        rej, allAdjPvals, _, _ = multipletests(allpvals, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

        for eidx, elem in enumerate(foundScores):

            returnResult.append((elem[0], elem[1], zScores[eidx], allpvals[eidx], allAdjPvals[eidx]))

        return returnResult

    def getTopRegulatedMIRNAS(self, graph, direction, maxNodes=50, measured=False, significant=False, numDeGenes=0, numAllGenes=0):

        node2ratio = {}

        for node in graph.nodes():

            nodeData = graph.nodes[node]["attr_dict"]

            if not nodeData["type"] == "mirna":
                continue

            nodeDirection = nodeData.get("node_expr_direction", "N/A")

            if nodeDirection == "N/A":
                continue

            if nodeDirection.upper() != direction.upper():
                continue

            nodeNN = [x for x in graph.neighbors(node)]

            totalDirectionCounter = Counter()
            directionCounter = Counter()

            for nn in nodeNN:
                tgt = nn
                tgtData = graph.node[tgt]["attr_dict"]
                tgtDir = tgtData.get("node_expr_direction", "N/A")

                totalDirectionCounter[self.dirToOpposite[tgtDir]] += 1  # opposite, because this is evidence for opp dir

                if measured and tgtData.get("de_measured", "false") == "false":
                    continue

                if significant and tgtData.get("adjPval", 1.0) > sigThreshold:
                    continue

                directionCounter[self.dirToOpposite[tgtDir]] += 1 # opposite, because this is evidence for opp dir


            divValue = directionCounter[self.dirToOpposite[direction]]
            directionCount = directionCounter[direction]

            ratio = float(directionCount) / float(divValue if divValue > 0 else 1)

            if ratio < 1:
                ratio = 0

            node2ratio[node] = (node, nodeData["de_measured"], ratio*directionCount, ratio, directionCount , divValue)


        allNodeNames = [node2ratio[x][0] for x in node2ratio]
        allNodeScores = [node2ratio[x][2] for x in allNodeNames]


        if len(allNodeScores) == 0:
            rej = []
            allAdjPvals = []
            allNodeScoresZScores = []
            allpvals = []

        else:
            allNodeScoresZScores = stats.zscore(allNodeScores)
            allpvals = [scipy.stats.norm.sf(abs(z)) for z in allNodeScoresZScores]
            rej, allAdjPvals, _, _ = multipletests(allpvals, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

        ovaPval = []
        ovaAdjPval = []
        rejOva = []

        for x in node2ratio:
            data = node2ratio[x]
            populationSize = numAllGenes
            numSuccInPopulation = numDeGenes

            drawnSuccesses = data[4]
            sampleSize = len(graph.node[data[0]]["attr_dict"]["all_targets"])

            pval = hypergeom.sf(drawnSuccesses - 1, populationSize, numSuccInPopulation, sampleSize)
            #print(x, drawnSuccesses - 1, populationSize, numSuccInPopulation, sampleSize, pval)
            ovaPval.append(pval)

        if len(node2ratio) > 0:
            rejOva, ovaAdjPval, _, _ = multipletests(ovaPval, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)



        header = [("miRNA", "Measured?", "Score", "Measured Ratio", "Measured {} elements".format(direction), "Measured {} elements".format(self.dirToOpposite[direction]),
                     "zscore", "pval", "adj_pval", "genes", "strict_opp_dir", "ova_pval", "ova_adj_pval", "mir_target_count")]
        retRatio = []
        for x in node2ratio:
            data = node2ratio[x]
            savedElem = list(node2ratio[x])

            nodeIdx = allNodeNames.index(savedElem[0])

            savedElem.append(allNodeScoresZScores[nodeIdx])
            savedElem.append(allpvals[nodeIdx])
            savedElem.append(allAdjPvals[nodeIdx])

            genes = self.getConsistentNeighbors(graph, savedElem[0], measuredOnly=measured, significantOnly=significant)

            strictGenes = self.getConsistentNeighbors(graph, savedElem[0], measuredOnly=measured, significantOnly=significant, strict=True)

            savedElem.append(";".join(genes))
            
            strict_opp_dir = len(strictGenes)
            savedElem.append( str(strict_opp_dir) )

            savedElem.append( ovaPval[nodeIdx] )
            savedElem.append( ovaAdjPval[nodeIdx] )

            sampleSize = len(graph.node[data[0]]["attr_dict"]["all_targets"])
            savedElem.append(sampleSize)

            retRatio.append(tuple(savedElem))


            assert(len(savedElem) == len(header[0]))

        bestNodes = sorted(retRatio, key=lambda x: x[12], reverse=False)

        if len(bestNodes) > maxNodes and maxNodes != -1:
            bestNodes = bestNodes[0:maxNodes]

        return header + bestNodes

    def writeTopRegulatedMIRNA(self, graph, outfilename, upregs, downregs):
        """

        :param basename:
        :param upregs:
        :param downregs:
        :return:
        """

        outfile = open(outfilename, 'w')

        """
        elem_id	population_size	success_population	sample_size	success_samples	sample_success_fraction	pval	adj_pval	direction   genes
        ubiquitination	58253	3941	9418	1777	0.1886812486727543	0.0	0.0	UP H1-3;CENPK;CLCN2;CCNB1;KANK2;BCAT1;PLAAT3;MYO5B;MPP2;ARHGAP9;SLC41A3;SYNE1;EGFLAM;CNDP2;RNFT1;SORT1;TNPO2;USP46;RASGRP2;MID2;POLR1E;MAP1B;MXRA5;SLC38A1;NEK7;KNL1;CENPN;CORO2A;CRYAB;ATF5;SRD5A3;TMEM9B;MTBP;MXRA7;SLFN11;MCF2L2;IQGAP2;SCCPDH;SLC7A11;TPM4;MSL3;ARPC3;NCEH1;SP140;FANCD2;MTMR2;PSMB9;ANTXR1;PYCARD;SPIN4;BCL7A;ZNF367;HOOK3;SEC23A;UTP18;CKS2;ARRB2;HDAC11;SUSD1;TSR2;KIF20A;SMOX;PARP12;H2AC8;VCL;GSTO1;GRAP;LRRC8C;KYNU;COLGALT1;MALT1;BAG2;FAM83D;RILPL2;DENND5B;TESPA1;TJP2;MAN1A1;SLC2A3;YIPF1;HIP1;FADS2;LRP6;DHCR24;MYLIP;SPINT1;GFPT1;FLNA;EZR;MDFIC;TBC1D4;TANK;ANGEL1;DDHD1;CD48;TMPO;KANK1;CBX8;PNMA2;CLDN7;BTBD1;MGME1;ROCK2;GJA1;SLC41A1;NUDT6;PIGU;SPAST;HINT3;TNFRSF11A;ADRB2;HLA-E;IL2RG;ADGRB2;DGLUCY;DENND1C;RAB9B;FGFR1;ADD1;HIVEP2;OXSR1;GPRASP2;KLHL7;FLT4;APBB1;LGALS8;GK;CAPZA1;PKN1;ATP11A;LFNG;DIP2C;LBR;CEP128;C5AR1;CTTN;MRPL49;MICB;RIMS1;LSS;PTPN11;DYSF;HLA-F;LDLRAD4;EFCAB13;RBPMS;GBA2;SGO2;SPDL1;FXYD5;PLAUR;PIK3R3;ZSWIM6;COA6;DCP2;ACTR1B;MED11;GNA11;IVNS1ABP;NFE2L3;LCP2;ARL4C;MBNL1;TJP1;PWWP3A;ERAP2;LAYN;ATP2B1;STARD7;SLC2A12;HLF;ID1;LRIG1;SKAP2;TMEM94;PPARGC1B;BAIAP2;PGK1;HIF3A;KPNA2;AOPEP;MAP4;TEAD2;PDK3;NCS1;SPINT2;UBE2D4;FSTL3;DYNC1LI2;KLHL21;METTL14;ZNRF2;ANKS6;GRM5;BID;TSPYL5;CRLF3;UBXN11;TYW1;SLC35F6;PPY;SVIL;PHF6;RGPD8;IRS1;TXNDC12;BUB1B;GOLGA2;ITGA7;STRN3;ANLN;DNM1;MTM1;MIOS;SPTB;TVP23C;ANO6;ZEB1;LDLRAP1;AKAP1;TRIM24;ABCC4;GPRC5C;VPS8;SHCBP1;SYTL3;ITGAL;CADM1;TRIP10;SH3BP2;ARHGAP19;DNAH17;TULP3;IKBKE;VAC14;PCDHB15;CD83;PRDX3;ARHGAP1;PDE5A;DMPK;MET;CNNM2;HELLS;NPRL2;RAB34;PLCG1;P3H2;ZCCHC18;MAGED2;TALDO1;PELO;PSMA6;SLC22A5;GMFG;CHST11;MFNG;PRKAR2B;TMEM33;MX2;NUDT10;TDP2;ARHGEF39;PACC1;BLVRB;ITPR1;ZDHHC7;ZDBF2;INPP5D;SGK1;EFHD1;ASB1;DSTYK;VCAN;ANKRD42;SPRYD3;PLEKHA1;ATP9A;IL27RA;BLVRA;H2BC12;GADD45A;CMTM6;TOP2A;TBC1D1;DIP2B;RBMS2;KIFC1;CD151;PRKRA;COL18A1;SMAD6;TMEM120A;SLC12A7;EGLN3;PSMB8;DENND4A;GPAT3;KITLG;WFS1;TREX1;F11R;SCRN1;NPTN;SNX25;TRPV2;STK4;CDK5RAP2;SFMBT2;ANKRD37;TUB;SELENOF;TRERF1;ATAD5;HS3ST1;MFAP3;ISG15;ATP1A2;ARNTL;NACAD;BEX3;STK17B;CDK1;NFIX;SCARF2;PLS3;ARVCF;HHAT;STK19B;KIAA0355;MACROH2A1;MYO7B;GNB3;MCM7;VAMP8;ZSCAN18;LPXN;RAPGEF1;MTX3;RPS6KA1;VDR;BRCA1;RFTN1;PLEKHO2;ATP6V1A;LYN;CTSL;EHD3;RALGAPA2;FEN1;PLA2G7;SMG7;WDR91;NUSAP1;RRAS;USP11;UNC93B1;FCHO2;KCTD10;SELPLG;CCND2;TEX264;BCAT2;ADCY8;CD180;KIF5B;HLA-A;MPP1;ARRB1;RAVER2;LGALS9;EXOC6B;PPTC7;PMEPA1;FADS3;CSRNP2;NUDCD3;ENAH;H1-5;CBX6;CEP170;PLAGL1;BMPR2;HPCAL1;MAPKAPK3;RAP2B;B3GALNT2;MED13;ORMDL3;NETO2;ST7;DBN1;PHF13;RPGR;ARHGAP15;TSPYL1;ACE;HIP1R;LYST;SLA;PM20D2;ABCB1;RPL28;TAF10;APBB2;KNTC1;CREB3L2;ATP6V1C1;NAP1L2;PEX19;SMIM15;TSPAN33;IL17RA;NRSN2;GAMT;POLR2B;SYNPO;PCDHB16;ITM2C;GNAO1;PLCD3;ADAM15;TMEM263;BMP2K;SNTB1;PLIN2;TCEAL2;EPB41L5;LAMC1;GLUL;PPP1R3C;NAB2;PLBD2;TMEM51;SPECC1;RGS19;DOP1B;ATP9B;FRY;RGS14;MAPK7;LAIR1;NABP1;ALDH2;SUDS3;ARMT1;PNPLA2;MYO1D;LAP3;MTURN;RBBP8;BNC2;CRY2;CRMP1;TPM1;FLVCR2;CDC7;GGT7;SLAMF6;LRRC70;CD84;SUMO3;COL15A1;TNRC6A;GINS2;GIMAP2;SLC43A3;LMNB1;KCNA3;TTLL10;MIEF2;HPGD;TARBP1;KRT17;GNPTAB;CEP55;MYH10;ING5;CUTC;SMARCA1;NOP16;SLC39A11;MYO5C;PEPD;RBL1;CYFIP2;IFI44;IFT43;FOXN2;PHLDB1;KIF2A;LIMK1;KCTD12;TRAF5;SMTN;DCLRE1C;SHKBP1;HMMR;GATB;FSCN1;RNASET2;NCAPG;PARK7;RB1;EHBP1L1;CDCP1;GUSB;HM13;PCK2;OLFM2;TUT7;ZNF217;NCSTN;POLR2L;UNC13C;DCK;NDN;PLCL2;CASP7;FUCA2;ZBTB24;NAT8;HPLH1;TCEA2;PAK3;UTP14C;STOML1;EML2;SBDS;H2BC11;CLTB;EGFL7;EMC3;NAP1L3;SMARCD3;STARD4;CAMK4;MAGEE1;BIN2;NBEAL1;CD4;RPS6KA6;LARP1B;LRP5;FKBP1A;NCKIPSD;FAM53B;TCF7L1;RTL8C;ZNF470;ATP6AP1;GLT8D2;NIBAN1;TOM1L2;JAM3;RNASEH2B;GAB2;MICAL2;LRP8;TNS2;SYTL4;FTL;WWTR1;CPED1;ITGA8;HIGD1A;SERPINA2;MMP15;DOCK8;AMPD3;NDRG1;ASAH1;UGP2;RAB39A;HIGD2A;CPVL;CNOT6;KIF20B;IL17RC;ARMCX3;CCDC190;SKI;SPTBN1;CYBB;PPP1R13B;EPM2A;MFSD6;TLE4;PDLIM1;NAA25;KIF2C;STEAP4;RETREG1;SSR1;SLC9A7;CLSPN;NLK;RBM47;TRIM32;HCST;KIF24;MLX;CSF1R;ALG6;FDXACB1;CACFD1;TIGD1;HOXB7;KRT16;TMEM209;CAVIN3;HMOX1;ARL6IP1;WDR76;DEPDC7;AMOTL2;CALHM5;INPP4B;CXorf56;STARD3NL;REV3L;SPSB1;MAN2B1;UBE2J1;HOMER1;TNC;ANKRD28;SNRPN;CYP1B1;ADAM17;PARP14;METAP1;ESCO2;ADCY9;RACGAP1;HENMT1;ITGA4;GDI1;SMURF1;PIK3CG;RRAS2;PPP1R3B;CLCN3;B3GNT2;STIL;RUNX2;IGF2R;GRK3;SLC12A4;KLF5;FCER1G;SAMD9;PTGER4;ATG9B;C19orf47;CDS2;CYB5R3;TSHZ3;SCGB1D1;TOX4;CRIP1;NOTCH3;ACTN1;BEX5;KRT86;PRKCH;GAS8;PLA2G4A;KIF14;DAG1;TCF12;TMEM245;ANKS1B;ERI1;RCN2;PLEK2;EIF4EBP1;CFL2;CLASP2;PTPRE;IQGAP3;UBFD1;GLS2;CISD2;MAP2K6;MOAP1;BLOC1S2;CTSB;ELF4;ATG3;RPLP0;FYB1;IL18;TMEM150C;ANO5;RYR2;TMBIM6;SH3BP5;CYTH2;BMP6;BCL11A;C16orf54;CYB5R2;CSK;MGAT5;TIMELESS;TPMT;PLEKHG3;GM2A;POFUT1;BCAS3;NAGA;ZNF267;INTS7;MCM9;ELOVL4;FKBP15;SGPL1;SYNPO2;SFT2D1;KCNAB2;MYL6;NAPG;CLIC4;SYNJ2;DPF2;SLMAP;ZNF135;ADAMTS1;ABCE1;FILIP1;EHD2;BAMBI;LACTB;FAM114A1;NRBF2;SEL1L3;RMI1;WDR36;HHEX;MAFK;HSD3B7;C12orf45;OTUD7B;RAD51B;FANCC;DENND2D;GPRC5A;CERS6;CYP7B1;CALCOCO1;TPST2;ABHD3;PARPBP;VRK1;HTATIP2;ANKEF1;PDLIM7;LHFPL2;IRF4;CTSS;RPS6KA3;EFEMP1;PDE7A;IFIH1;OASL;SUN2;PAQR3;MED7;USP20;RIPK2;CASP10;ITGA3;SULF2;VAMP2;HRH1;TWSG1;SLC30A7;ACP2;ATPAF1;PGD;KRT8;SOAT1;GOLGA3;FABP5P1;BRIP1;PTPRJ;BTN3A2;CSPG4;HES4;SLC40A1;UPP1;PACS1;MORC4;CHST15;SKA2;THRA;ZGRF1;WWP2;CARMIL1;CDK19;LYVE1;CENPF;TACC3;PGR;OAS3;AGFG1;OSBPL10;USO1;RCAN3;NLN;MAP3K1;PLAAT4;OXTR;CCNE2;CTSD;ACVR2B;RAC2;F8;PRMT2;ZNF92;PHB2;SH3BGRL2;PLK4;HACD2;GNS;TNFSF10;COLQ;BCR;ITGA6;MCU;UBA2;SLC20A2;NPR1;CENPO;SAMSN1;MCM4;TRIM3;PSAP;SLC16A10;TRIM25;KIRREL1;HEATR3;MTARC2;PARVA;CCDC18;ILK;SELENON;WDFY4;TUBA4A;KRT7;AIF1;PTBP3;MPP7;CXorf38;PARP15;EFNA1;MPRIP;MFGE8;PAWR;FDX1;BLM;PDE2A;PTPN7;ETS1;GCLM;ANXA2;ARHGAP18;STIM2;ENDOD1;RGS5;CLEC2B;RANBP17;SORBS2;CRYBG1;FBH1;MPHOSPH9;TPP1;PPP3CA;ACOT13;WWP1;H1-10;TMED5;CIAO2A;EDEM1;DAPK3;PPP1R12A;DNAJB5;PAQR4;EPB41L3;PIAS3;UBE2D1;PAM;COMMD10;DNMT1;CYTH1;RESF1;SEPTIN10;SH3BGRL3;MAGEF1;SULT1A1;TCEAL9;BPNT2;ALDH5A1;RRM2;RAB23;PIGW;RNF217;SUGCT;NFKBIE;FECH;LIN7A;ORAI3;BMPR1A;TLN1;NCKAP1L;QPCT;ACTA2;PDK2;AMIGO2;ANTXR2;TRIP13;PLXDC1;MACF1;MYO5A;ZMYND11;PYGL;ACTN4;ZBTB16;MYH7B;CAVIN1;AVEN;ECT2;GIMAP6;NUP50;CPSF3;NSMAF;TCEAL4;KIFBP;HMGB3;PRDM16;SH3GLB2;ITPRIPL1;ARHGAP26;OBSL1;EMILIN2;VRK2;SNAP25;RAD51AP1;RNF185;GUF1;C20orf27;MAPK13;PCDH7;WEE1;ZNF551;TCEAL3;RAPGEF4;UBE2H;BCL10;HOMER2;RHOBTB3;TPX2;FLT1;INPP5A;SYDE2;ZFP36L2;DDAH1;DPYD;TMEM237;FLI1;NDE1;NCOA4;CNN1;PURA;COMMD9;RNF115;NCAPH;FANCH;BCAR1;ITGB1BP1;PIK3AP1;CKAP2;RIPK4;HSPBAP1;MYH11;RNF138;ZFPL1;GAS2L3;ERBB2;EVI2A;CNKSR3;JAG1;BTK;TRAF3IP3;PTK2B;IRF8;ARPC1A;CD276;OAS2;RIC3;HSPB8;TRIM56;NDUFB9;MICAL3;SINHCAF;RALA;CYTH3;ST14;SACS;SERPINB8;POLR2C;IMPA2;HYAL2;PRCP;ERGIC1;PRC1;ATP6V1F;EGFR;SRPK3;NT5DC1;TSPAN14;MFAP3L;GALNT18;GIMAP7;CHPF;SLC43A2;PIK3R5;PHF1;MOB1B;ATP6V0D1;FIGN;RAB5IF;CALD1;LIX1L;IRAK2;IDH1;P3H3;RASSF2;COBLL1;CNPY3;CD109;CCSAP;DDX50;BIRC3;WASHC5;ARHGEF3;GNPTG;BMERB1;FIGNL1;MPDU1;NCAPG2;SLC10A7;TRABD2B;SIRPA;STK38L;AP1S3;C9orf72;STEAP3;PLA2G15;MAD2L1;ASRGL1;TMEM43;ALG11;BAZ1A;TMEM159;GPNMB;SLC35D1;CMIP;NCDN;C1orf198;RDH10;RNPEP;DCLRE1B;GINS4;PPP2CB;TNS1;SMURF2;DDX24;HEXIM1;NDEL1;PRKAB2;HAS2;NEXN;DCLK1;EPN2;EPHX2;RUSF1;RBMXL1;TOM1L1;TMEM102;SMOC2;S1PR3;NR2F6;ZEB2;WAS;TERF2IP;CDC42SE2;SNURF;LRR1;AGPAT5;VASN;FAM81A;ASPH;GCLC;RAB3GAP1;CREB3L4;GOLM1;OPTN;ATP6V1B2;LIMCH1;OSBPL5;BSDC1;SEMA4C;COQ8A;CTSC;NDUFA4;PTPN6;FAR2;LYRM2;HSPG2;ARHGAP21;RNF130;TUBGCP3;RIMKLB;SLC25A1;CSRP1;CLIP3;ZWINT;MTR;OAT;OTULINL;FXYD6;DTX3L;OGFRL1;NRP1;GMNN;ID3;ABHD6;CYRIB;FZD6;CEMIP2;ARL4A;B2M;NUP210;RNF125;RNF13;TSPYL2;ABLIM1;NKIRAS1;TSPO;FYCO1;HSPA2;TOPBP1;BBS1;EIF5A2;PECAM1;KIAA1794;PELI3;PINK1;PTPRU;NUP62;FBXW4;CARD8;GNG2;H3C8;NOS3;PKN3;FKBP9;AMPH;GSR;ARHGAP11A;SNX5;PTGES3L-AARSD1;RAB8B;RAB21;SGCB;KCTD9;PTTG1;TAL1;DPYSL3;PHKA1;SCAMP3;C19orf54;CRIM1;LNPK;PRIM1;VAV3;SMARCC2;KLHL42;HCK;LONRF2;UNC13D;DMXL2;CAPG;NADK2;SLIT3;GLRB;RORA;BTG1;AP2S1;METAP2;P3H4;TMSB10;DUSP19;IGFBP2;HLX;FTH1;CD55;ICAM2;LSM11;KLHL5;AKAP2;PLXNC1;PIM2;CRTAP;AVL9;SECTM1;TTC39C;ADCY6;DIRAS1;GLI3;ANO1;SIRPB1;PRXL2B;POLQ;LYPLA1;CHST14;MAP3K5;MTMR14;MTHFD1L;GKAP1;MYL6B;SDC4;ERCC6;SGK3;YAP1;DNAJB4;SYK;ADAMTS5;ANK1;COL17A1;SPAG5;H1-4;RASIP1;LPP;RREB1;RNF25;HACD1;CCDC8;CHAC1;LAPTM5;SLC19A1;RBMS3;POMGNT1;ANXA11;SNRK;ESYT2;TAGLN;EFNB1;CTPS1;CAPN5;ZNRF3;SLC6A6;CERT1;TTYH3;SPRYD7;DUSP4;STN1;PACSIN3;GPRASP1;TRAPPC5;MFSD1;HSPA4L;CBR4;SKOR1;DAB2;PGRMC1;NDC80;NFIL3;TRIP6;NMI;ACOX3;DAB2IP;PSEN1;SOX4;TYMS;CCDC102A;TBC1D31;ANXA3;UBAC1;DCAF12;PYGB;MATK;C12orf75;ZNF710;MAGEH1;TADA2A;ELF2;GMFB;ELOF1;CHD7;ANK3;PIP4K2B;PKIB;NQO2;ELMO1;ZNF711;DNAJC21;PXDN;GIT2;AEBP1;PCNA;HERC6;SPRED1;TRAT1;HLA-DRA;LRRFIP2;RHOB;ADGRG1;TRAPPC6A;CCDC88A;RASSF3;URGCP;MAP4K1;DERA;SLA2;AGPS;ADI1;TNFAIP8L2;LZTS2;DDR2;DYNLL2;RPS6KA5;NAP1L5;ADAM10;RNASEK;AUNIP;DMTN;SGO1;ITGB4;NOD2;NEDD4L;CA2;LMTK2;SLC26A6;APBB1IP;CCDC138;PHTF2;FERMT3;STK32C;SCA12;HAUS7;BTN3A3;MPP6;C2CD2;PATJ;FYN;MCM5;RNF8;ALAD;PGAP4;RUSC2;HSPB1;IRF7;ERBB4;EPHB4;FKBP10;C9;PDLIM4;TCIRG1;SMAD9;B3GNTL1;SNX10;N4BP2;SH3TC1;DOCK11;SNCG;ITGB1;ZNF888;MTRF1;ARHGAP44;ST8SIA4;MREG;HLA-B;KCNN4;PRKD2;MSH5;BUB1;PCSK6;MSH2;CLEC7A;ELK3;TIMP1;HECW2;SQOR;LRRC8B;ADAM9;ERMP1;CRYBG3;RTL8B;C5orf24;ARL11;TXN;GK5;CLIC2;ATL3;JMJD6;TAOK3;CDC42EP1;ZMIZ1;PCYOX1;ACSL5;LSM14A;BRI3BP;NRCAM;ZNF397;MUS81;SCN3A;IFI16;RUBCNL;SLC25A23;TCAF1;CD5;CLIC1;MPDZ;IL1RAP;HS2ST1;MAP4K2;ZBTB33;SV2A;EVA1C;HVCN1;DDR1;OSTF1;NRROS;ARHGAP30;CIP2A;CYBA;TUBB3;TPM3;GINS3;DUSP14;NECTIN3;SIRPG;GMCL1;LUC7L2;LTBP4;DDI2;MYH9;RFC3;ACAP2;NPC1;TPM2;MBNL3;RBPJ;SELENOI;KBTBD8;PFN2;FAM111A;DOCK2;CEP85;SNX18;SEPTIN6;CENPE;TENT2;SLFN5;FOXM1;GSDME;KIF4A;TRIM5;SERPINB1;WRNIP1;COL4A3;PHLPP2;PCBP4;SCARA3;ARL15;POGLUT3;NME3;BCAP31;SNX33;SLC25A4;NCKAP1;LAMA5;RNASE1;WDR11;AMOT;S100A11;MESD;GNG12;GRK5;PFKM;IKZF1;PDXK;APOBEC3G;CSTA;ZBTB5;PTPRC;NPRL3;GDPD1;IRF5;PASK;ANKRD6;CYREN;CSTF2T;PLEKHB2;FAM111B;GGT1;SLC4A8;POLR3D;CTSZ;RAI2;SHTN1;IPMK;C1GALT1;ITGB5;MACROD2;TXNIP;MX1;LPL;TMX1;GIMAP4;GFPT2;FRAT2;S1PR1;GLRX;ASAP2;ASXL2;CHST3;TANC1;VGLL4;CPNE8;IFI30;LRRC8D;TRO;FERMT2;PDLIM5;UGCG;GNL1;CASP8;CARD10;B9D1;AEBP2;MORF4L2;SLC7A8;PRKAA2;SP110;MAOB;NGLY1;AMOTL1;TES;DSTN;NSF;AKR1A1;PPIF;STAMBPL1;NPC2;CIRBP;ITGB7;OGN;ARHGAP23;SAMD9L;LEF1;TACC2;ZFP36L1;TWF2;LGALS3;OSBPL8;CAMK2G;ABCG2;PELI2;KIF11;OSMR;SPOP;TEC;IPP;TDP1;PARP4;ADH5;ENTPD7;CELIAC3;PEBP1;IFNGR1;ODF2L;EIF4G3;GALNT2;UBASH3B;MRPL36;CATSPERB;CDH2;CHST7;TFDP2;EEPD1;MPG;MTMR1;SDCBP;FAM177A1;RUSC1;NFE2L1;TREM1;DUSP5;GAB3;CDCA4;UHRF1;GSTK1;GNA13;HMGA1;COTL1;IFIT3;CHEK2;CEP152;CLDN1;SLC7A2;CORO7;PAK2;MERTK;C11orf96;PPT1;SAMHD1;GPD2;RTN4RL1;EPHB2;CDK17;TBC1D9;PLD3;LRRC8A;TNFRSF10A;ITK;GALNT1;CEP135;ETNK1;PTPN3;QKI;ESPL1;R3HCC1L;FLNC;DMWD;PTDSS1;NRAS;SLC16A4;PCDHB11;OXCT1;RUNX1;DLG5;APOL6;MAP3K20;CLASP1;UAP1L1;CD33;PI4K2B;SEPSECS;OSBPL11;TAF1;OSBPL3;DYNLL1;MYLK;SLC45A1;SERGEF;EIF4A1;NXT2;LINS1;DECR1;ALS2CL;EMD;PLEKHH3;SLC37A2;MYL9;MTFR1L;SH2D3C;CD99L2;MKI67;ALDH1A2;FBXO32;HEYL;SEMA4A;ANP32E;SEZ6L2;REEP3;EDIL3;PTER;PRR11;CNN3;ATP10D;B4GALT2;AKAP12;KIAA0930;CMPK2;IFIT1;PIP5K1A;SH3RF1;NCAM2;DOCK4;SLC39A10;HIBCH;PDGFRA;PDE8B;PPP1R13L;AP3M2;SLC16A6;PNMA8A;LCORL;DNAJB2;SEPTIN8;TNFRSF12A;NFXL1;KLC1;SLC25A19;ATP10A;CDC6;SLC2A10;IGF1R;ZNF200;PARP9;CDCA8;LAGE3;DTL;RHBDF1;CCDC103;DEPDC1;H2BC5;DOK2;BTN3A1;SNX6;PLXND1;KIF7;PBXIP1;ZFYVE16;DSC3;MIB1;ISYNA1;AP4M1;CAND2;TMEM98;NRIP1;ERG;AK1;PHF10;KIF15;S100A10;DNAH14;MYO1E;SH2B3;HABP4;STK17A;UHMK1;TMEM132A;NEURL1B;SGCE;PLEKHG2;EID1;H2BC3;SOWAHC;TLN2;PKMYT1;MIPEP;ASPM;NACC2;TGFB1I1;NIPA2;ERLEC1;FMNL1;SERINC3;EPPK1;MPHOSPH8;SHMT1;NEO1;PTPMT1;CSTB;PLCD1;GINS1;IRF6;NDRG3;CD82;CCNA2;RNF41;SUN1
        """

        header = ["elem_id","measured","score","ratio","dir_elements","opp_dir_elements","pval","adj_pval", "direction", "genes", "strict_targets", "target_gene_count", "all_target_genes", "ova_pval", "ova_adj_pval"]

        print("\t".join(header), file=outfile)

        #reg elems
        #miRNA, Measured, Score, Ratio, dir_elements, opp_dir_elements, total_ratio, total_dir_elements, total_opp_dir_elements, zscore, pval, adjpval, genes


        #allpvals = [scipy.stats.norm.sf(abs(elem[9])) for elem in upregs[1:] + downregs[1:]]
        #rej, allAdjPvals, _, _ = multipletests(allpvals, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

        directions = ["UP"] * len(upregs[1:]) + ["DOWN"] * len(downregs[1:])

        for idx, elem in enumerate(upregs[1:]+downregs[1:]):
            miRNA = elem[0]
            pop_size = elem[1]

            success_pop = elem[2]
            sample_size = elem[3]
            success_sample = elem[4]
            success_fraction  = elem[5]

            pval = elem[7]#allpvals[idx]
            adj_pval = elem[8]#allAdjPvals[idx]
            direction = graph.node[miRNA]["attr_dict"].get("node_expr_direction", "N/A").upper()

            all_targets = ";".join([str(x) for x in graph.node[miRNA]["attr_dict"]["all_targets"]])
            target_gene_count = len(graph.node[miRNA]["attr_dict"]["all_targets"])

            genes = elem[9]
            strict_count = elem[10]

            ovaPval = elem[11]
            ovaAdjPval = elem[12]

            if type(genes) == list:
                genes = ";".join(genes)

            print(  miRNA, pop_size, success_pop, sample_size, success_sample, success_fraction, pval, adj_pval,
                    direction, genes, strict_count, target_gene_count, all_targets,
                    ovaPval, ovaAdjPval,
                     sep="\t", file=outfile)

        outfile.close()

    def getConsistentNeighbors(self, graph, mirna, measuredOnly=False, significantOnly=False, strict=False):

        nodeData = graph.nodes[mirna]["attr_dict"]

        nodeDir = nodeData.get("node_expr_direction", "N/A")

        if nodeDir == "N/A":
            return []

        nodeNN = [x for x in graph.neighbors(mirna)]

        consistentNeighbors = []
        isExplainedByOtherNode = False

        for nn in nodeNN:
            tgt = nn
            tgtData = graph.node[tgt]["attr_dict"]
            tgtDir = tgtData.get("node_expr_direction", "N/A")

            if measuredOnly:
                if tgtData.get("de_measured", "false") == "false":
                    continue

            if significantOnly:
                if tgtData.get("adjPval", 1.0) > sigThreshold:
                    continue

            nnNeighbors = [x for x in graph.neighbors(nn)]

            nnExplained = False
            for nnn in nnNeighbors:
                if nnn == mirna:
                    continue

                nnnDir = graph.node[nnn]["attr_dict"].get("node_expr_direction", "N/A")
                if tgtDir != "N/A" and tgtDir == self.dirToOpposite[nnnDir]:
                    nnExplained = True
                    break
            
            isExplainedByOtherNode = isExplainedByOtherNode or nnExplained


            if tgtDir != "N/A" and tgtDir == self.dirToOpposite[nodeDir]:

                if not strict or isExplainedByOtherNode:
                    consistentNeighbors.append(nn)

        return consistentNeighbors

    def saveGraph(self, outpath, outhtml, graph, numDeGenes, numAllGenes):

        graphStats = {}

        graphStats["Node Stats"] = self.to_print_format(self.countAttributeNodes(graph, ["type", "node_expr_direction", "node_expr_detection"]))
        graphStats["Edge Stats"] = self.to_print_format(self.countAttributeEdges(graph, ["edge_type", "edge_creation"]))
        graphStats["Edge Stats (2)"] = self.to_print_format(self.countAttributeEdges(graph, ["dir_type", "edge_creation"]))

        graphStats["Unexplained Genes"] = self.getUnexplainedGenes(graph)
        graphStats["Measured Inconsistencies"] = self.getMeasuredInconsistencies(graph)
        graphStats["Remaining Imputed Inconsistencies"] = self.getImputedInconsistencies(graph)

        self.scoreGenes(graph)
        self.scoreMIRs(graph)
        graphStats["miRNA NS Scores"] = self.getNSScores(graph)

        graphStats["UP regulated miRNAs"] = self.getTopRegulatedMIRNAS(graph, "up", -1, numDeGenes=numDeGenes, numAllGenes=numAllGenes)
        graphStats["DOWN regulated miRNAs"] = self.getTopRegulatedMIRNAS(graph, "down", -1, numDeGenes=numDeGenes, numAllGenes=numAllGenes)

        """
        write out TSVs
        """

        basename = os.path.abspath(outpath.name)
        basename = os.path.splitext(basename)[0]

        upregMirs = self.getTopRegulatedMIRNAS(graph, "up", -1, measured=True, numDeGenes=numDeGenes, numAllGenes=numAllGenes)
        downregMirs = self.getTopRegulatedMIRNAS(graph, "down", -1, measured=True, numDeGenes=numDeGenes, numAllGenes=numAllGenes)
        self.writeTopRegulatedMIRNA(graph, basename + ".mirs.tsv", upregMirs, downregMirs)

        upregMirs = self.getTopRegulatedMIRNAS(graph, "up", -1, measured=False, numDeGenes=numDeGenes, numAllGenes=numAllGenes)
        downregMirs = self.getTopRegulatedMIRNAS(graph, "down", -1, measured=False, numDeGenes=numDeGenes, numAllGenes=numAllGenes)
        self.writeTopRegulatedMIRNA(graph, basename + ".mirs.all.tsv", upregMirs, downregMirs)


        self.gp.saveGraph(graph, outpath, name=outhtml, stats=graphStats, title=outhtml + " {}".format(datetime.datetime.now()))

    def colorInitial(self, graph):

        for edge in graph.edges():

            srcNode = edge[0]
            tgtNode = edge[1]

            srcNodeData = graph.nodes[srcNode]["attr_dict"]
            tgtNodeData = graph.nodes[tgtNode]["attr_dict"]

            graph.edges[edge]["linestyle"] = "solid"
            graph.edges[edge]["edge_type"] = "unknown"
            graph.edges[edge]["dir_type"] = "unexpected"

            if (srcNodeData["type"] == "gene" and tgtNodeData["type"] == "mirna") or (srcNodeData["type"] == "mirna" and tgtNodeData["type"] == "gene"):

                reg0 = srcNodeData.get("log2FC", None)
                reg1 = tgtNodeData.get("log2FC", None)

                if reg0 == None or reg1 == None:
                    continue

                graph.edges[edge]["linestyle"] = "solid"

                if (reg0 > 0 and reg1 < 0) or (reg0 < 0 and reg1 > 0):
                    # expected regulation
                    graph.edges[edge]["color"] = "g"
                    graph.edges[edge]["edge_type"] = "expected"
                    graph.edges[edge]['dir_type'] = "expected"

                elif (reg0 > 0 and reg1 > 0) or (reg0 < 0 and reg1 < 0):
                    graph.edges[edge]["color"] = "r"
                    graph.edges[edge]["edge_type"] = "unexpected"
                    graph.edges[edge]['dir_type'] = "unexpected"



    def __impute1(self, graph):
        """

        given a miRNA -> if all targets are in one direction => miRNA opposite direction

        """
        # impute regulations
        for node in graph.nodes():
            nodeData = graph.node[node]["attr_dict"]

            if nodeData["type"] != "mirna":
                continue

            if nodeData["de_measured"] == "true":
                # we cannot change such a node ...
                continue

            nodeNN = [x for x in graph.neighbors(node)]

            nnRegs = {}
            for nn in nodeNN:

                nndata = graph.node[nn]["attr_dict"]
                if nndata["log2FC"] < 0:
                    nnRegs[nn] = "down"
                elif nndata["log2FC"] > 0:
                    nnRegs[nn] = "up"
                else:
                    nnRegs[nn] = "N/A"

            mirDir = set([nnRegs[x] for x in nnRegs if nnRegs[x] != "N/A"])
            sameDir = len(mirDir) == 1

            if sameDir:

                exprVal = 0
                nodeDirection = "N/A"
                if 'up' in mirDir:
                    exprVal = -1
                    nodeDirection = "down"
                else:
                    exprVal = 1
                    nodeDirection = "up"

                graph.nodes[node]["attr_dict"]["node_expr_direction"] = nodeDirection
                graph.nodes[node]["attr_dict"]["color"] = self.exprColor(exprVal, nodeDirection)
                graph.nodes[node]["attr_dict"]['node_expr_detection'] = 'imputed1'
                print(node, "imputed")

                for nn in nodeNN:

                    if graph.edges[(node, nn)]["edge_type"] in self.priorityEdgeTypes:
                        continue

                    graph.edges[(node, nn)]["color"] = "g"
                    graph.edges[(node, nn)]["linestyle"] = "dashed"
                    graph.edges[(node, nn)]["edge_type"] = "imputed"
                    graph.edges[(node, nn)]["dir_type"] = "expected"
                    graph.edges[(node, nn)]["edge_creation"] = "imputed1"

    def __impute1_1(self, graph):
        """

        given a miRNA -> if all targets are in one direction => miRNA opposite direction

        here a neighbor direction may also be N/A
        => if a gene has no direction yet, it is not excluded

        """

        # impute regulations
        for node in graph.nodes():

            nodeData = graph.node[node]["attr_dict"]

            if nodeData["type"] == "mirna":
                continue

            if nodeData["de_measured"] == "true":
                # we do not want to impute miRs here
                continue

            nodeNN = [x for x in graph.neighbors(node)]

            nnDirs = set()
            acceptNode = True
            for nn in nodeNN:

                nndata = graph.node[nn]["attr_dict"]
                if nndata["log2FC"] < 0:
                    nnDir = "down"
                elif nndata["log2FC"] > 0:
                    nnDir = "up"
                else:
                    nnDir = "N/A"

                if nndata["de_measured"] == "true":
                    nnDirs.add(nnDir)
                else:
                    acceptNode = False
                    break

            if "N/A" in nnDirs or len(nnDirs) != 1:
                acceptNode = False

            if acceptNode:
                geneDir = self.dirToOpposite[list(nnDirs)[0]]

                graph.nodes[node]["attr_dict"]["node_expr_direction"] = geneDir
                graph.nodes[node]["attr_dict"]["node_expr_detection"] = "imputed1.1"
                graph.nodes[node]["attr_dict"]["color"] = self.exprColor(geneDir, geneDir)

                for nn in nodeNN:

                    if graph.edges[(node, nn)]["edge_type"] in self.priorityEdgeTypes:
                        continue

                    graph.edges[(node, nn)]["color"] = "g"
                    graph.edges[(node, nn)]["linestyle"] = "dashed"
                    graph.edges[(node, nn)]["edge_type"] = "imputed"
                    graph.edges[(node, nn)]["dir_type"] = "expected"
                    graph.edges[(node, nn)]["edge_creation"] = "imputed1.1"

                    print(node, nn, "imputed1.1")

        # impute regulations

    def nnHasEvidence(self, graph, node):
        nodeData = graph.node[node]["attr_dict"]
        nodeNN = [x for x in graph.neighbors(node)]
        for nn in nodeNN:

            nndata = graph.node[nn]["attr_dict"]
            #edgeType = graph.edges[(node, nn)]["edge_type"]
            #if edgeType in ['imputed', 'expected']:
            #    return True

            edgeDir = graph.edges[(node, nn)].get("dir_type", "unexpected").split("_")[0]
            if edgeDir in ['expected']:
                return True


        return False

    def __impute2(self, graph):

        """

        for a miRNA, locally optimize the number of edges

        impute in direction with most evidence-based targets

        => all evidences must have the same direction!

        """

        mirNodes = 0

        for node in graph.nodes():

            nodeData = graph.node[node]["attr_dict"]

            if nodeData["type"] != "mirna":
                continue

            if nodeData["de_measured"] == "true":
                continue

            mirNodes += 1

            nodeNN = [x for x in graph.neighbors(node)]

            nnRegs = {}
            nnEdgeTypes = Counter()
            nnWithEvidence = set()
            nnWithNoEvidence = set()

            if nodeData["log2FC"] < 0:
                mirDir = "down"
            elif nodeData["log2FC"] > 0:
                mirDir = "up"
            else:
                mirDir = nodeData.get("node_expr_direction", "N/A")

            for nn in nodeNN:

                nndata = graph.node[nn]["attr_dict"]
                if nndata["log2FC"] < 0:
                    nnRegs[nn] = "down"
                elif nndata["log2FC"] > 0:
                    nnRegs[nn] = "up"
                else:
                    nnRegs[nn] = nndata.get("node_expr_direction", "N/A")

                edgeData = graph.edges[(node, nn)]

                edgeType = edgeData["edge_type"]
                edgeDir = edgeData.get("dir_type", "unexpected").split("_")[0]

                nnEdgeTypes[edgeType] += 1

                if edgeDir in ["expected"] or self.nnHasEvidence(graph, nn):
                    nnWithEvidence.add((nn, nnRegs[nn]))
                else:
                    nnWithNoEvidence.add((nn, nnRegs[nn]))

            regCounter = Counter()
            for gene in nodeNN:
                regCounter[nnRegs[gene]] += 1

            # are there unassigned edges?
            if sum([regCounter[x] for x in regCounter]) != nnEdgeTypes["imputed"] + nnEdgeTypes["expected"]:

                # remove all elems with other evidence
                ncounter = {}
                for x in regCounter:
                    ncounter[x] = regCounter[x]

                for e in nnWithEvidence:
                    ncounter[e[1]] -= 1

                # check whether remainings can be imputed
                imputeDirs = [x for x in ncounter if ncounter[x] > 0]
                sameDir = len(imputeDirs) == 1

                undetEdges = set()
                newImputEdges = set()

                if node == "miR-23":
                    print("impute2", node)
                    print(ncounter)
                    print(regCounter)
                    print(nnWithEvidence)
                    print(nnEdgeTypes)

                if sameDir:

                    for e in nnWithEvidence:
                        undetEdges.add((node, e[0]))

                    for e in nnWithNoEvidence:
                        newImputEdges.add((node, e[0]))

                    undetEdges = list(undetEdges)
                    imputeDir = imputeDirs[0]
                    imputeDir = self.dirToOpposite[imputeDir]

                    print("impute2 samedir", node, sameDir, imputeDir, ncounter)

                    if imputeDir == 'up':
                        graph.nodes[node]["attr_dict"]['color'] = 'g'
                    else:
                        graph.nodes[node]["attr_dict"]['color'] = 'r'

                    graph.nodes[node]["attr_dict"]['node_expr_detection'] = 'imputed2'
                    graph.nodes[node]["attr_dict"]["node_expr_direction"] = imputeDir

                    for e in nnWithEvidence:
                        if e[1] == imputeDir:
                            undetEdges.remove((node, e[0]))
                            newImputEdges.add((node, e[0]))

                    for nn in nodeNN:

                        if graph.edges[(node, nn)]["edge_type"] in self.priorityEdgeTypes:
                            continue

                        if imputeDir == nnRegs[nn]:
                            graph.edges[(node, nn)]["color"] = "#D50000"
                            graph.edges[(node, nn)]["linestyle"] = "dotted"
                            graph.edges[(node, nn)]["edge_type"] = "imputed"
                            graph.edges[(node, nn)]["dir_type"] = "unexpected"
                            graph.edges[(node, nn)]["edge_creation"] = "imputed2"

                        elif imputeDir == self.dirToOpposite[nnRegs[nn]]:
                            graph.edges[(node, nn)]["color"] = "g"
                            graph.edges[(node, nn)]["linestyle"] = "dashed"
                            graph.edges[(node, nn)]["edge_type"] = "imputed"
                            graph.edges[(node, nn)]["dir_type"] = "expected"
                            graph.edges[(node, nn)]["edge_creation"] = "imputed2"
                        else:
                            print("error, no direction", imputeDir, nnRegs[nn], node, nn)

                else:  # not sameDir:
                    # if graph.nodes[node]["attr_dict"]["de_measured"] == "false":
                    #    graph.nodes[node]["attr_dict"]['color'] = '#FF8F00'
                    # else:
                    #    graph.nodes[node]["attr_dict"]['color'] = exprColor(graph.nodes[node]["attr_dict"]["log2FC"])
                    print("impute2 not samedir", node, regCounter, nnEdgeTypes, nnWithEvidence, sameDir, undetEdges)

    def nnGetEvidences(self, graph, node):
        nodeData = graph.node[node]["attr_dict"]
        nodeNN = [x for x in graph.neighbors(node)]

        evidenceNodes = set()
        for nn in nodeNN:

            nndata = graph.node[nn]["attr_dict"]
            edgeType = graph.edges[(node, nn)]["edge_type"]

            if graph.edges[(node, nn)].get("linestyle", None) == "dotted":
                continue

            if edgeType in ['imputed', 'expected']:
                evidenceNodes.add(nn)

        return evidenceNodes

    def __impute3(self, graph):

        """
            for a miRNA, locally optimize the number of edges

            may not create inconsistency with measured gene!
        """

        for node in graph.nodes():

            nodeData = graph.node[node]["attr_dict"]

            if nodeData["type"] != "mirna":
                continue

            if nodeData.get("node_expr_direction", None) not in [None, "N/A", "unknown"]:
                continue

            nodeNN = [x for x in graph.neighbors(node)]

            if len(nodeNN) == 0:
                print("Singlet Node", node)
                continue

            nnRegs = {}
            nnEdgeTypes = Counter()
            nnDirDirections = Counter()
            nnWithEvidence = set()
            nnWithNoEvidence = set()

            dir2Singletons = set()

            for nn in nodeNN:

                nndata = graph.node[nn]["attr_dict"]
                if nndata["log2FC"] < 0:
                    childDir = "down"
                elif nndata["log2FC"] > 0:
                    childDir = "up"
                else:
                    childDir = "N/A"

                nnRegs[nn] = childDir
                nnDirDirections[childDir] += 1

                edgeType = graph.edges[(node, nn)]["edge_type"]
                nnEdgeTypes[edgeType] += 1

                edgeDir = graph.edges[(node, nn)].get("dir_type", "unexpected").split("_")[0]

                if edgeDir in ["expected"] or self.nnHasEvidence(graph, nn):
                    nnWithEvidence.add((nn, nnRegs[nn]))
                else:
                    nnWithNoEvidence.add((nn, nnRegs[nn]))

            if len(nnWithNoEvidence) == 0:
                print("can be explained:", node, nnDirDirections, )

                childDir2Effect = defaultdict(list)

                for nn in nodeNN:

                    allEvs = self.nnGetEvidences(graph, nn)
                    nndata = graph.node[nn]["attr_dict"]

                    if nndata["log2FC"] < 0:
                        childDir = "down"
                    elif nndata["log2FC"] > 0:
                        childDir = "up"
                    else:
                        childDir = "N/A"

                    averageEffect = abs(nndata["log2FC"]) / len(allEvs) if len(allEvs) > 0 else 0
                    childDir2Effect[childDir].append((nn, averageEffect))
                    print("impute3, nnWithNoEv", node, nn, childDir, allEvs, averageEffect)

                etToEffect = [(x, sum([y[1] for y in childDir2Effect[x]])) for x in childDir2Effect]

                etToEffect = sorted(etToEffect, key=lambda x: x[1], reverse=True)

                finalTarget = etToEffect[0]

                targetDir = finalTarget[0]

                if targetDir == "up":
                    targetDir = "down"
                elif targetDir == "down":
                    targetDir = "up"

                print("impute3 intermediate", node, finalTarget, targetDir)

                if targetDir == 'up':
                    graph.nodes[node]["attr_dict"]['color'] = 'g'
                    graph.nodes[node]["attr_dict"]['node_expr_detection'] = 'imputed3'
                    graph.nodes[node]["attr_dict"]["node_expr_direction"] = targetDir
                else:
                    graph.nodes[node]["attr_dict"]['color'] = 'r'
                    graph.nodes[node]["attr_dict"]['node_expr_detection'] = 'imputed3'
                    graph.nodes[node]["attr_dict"]["node_expr_direction"] = targetDir

                for nn in nodeNN:

                    allEvs = self.nnGetEvidences(graph, nn)
                    nndata = graph.node[nn]["attr_dict"]

                    childDir = graph.nodes[nn]["attr_dict"].get("node_expr_direction", "N/A")

                    if childDir == targetDir:
                        if graph.edges[(node, nn)]["edge_type"] in self.priorityEdgeTypes:
                            continue
                        graph.edges[(node, nn)]["color"] = "#D50000"
                        graph.edges[(node, nn)]["linestyle"] = "dotted"
                        graph.edges[(node, nn)]["edge_type"] = "imputed"
                        graph.edges[(node, nn)]["dir_type"] = "unexpected"
                        graph.edges[(node, nn)]["edge_creation"] = "imputed3"

                    elif childDir == self.dirToOpposite[targetDir]:
                        if graph.edges[(node, nn)]["edge_type"] in self.priorityEdgeTypes:
                            continue
                        graph.edges[(node, nn)]["color"] = "g"
                        graph.edges[(node, nn)]["linestyle"] = "dashed"
                        graph.edges[(node, nn)]["edge_type"] = "imputed"
                        graph.edges[(node, nn)]["dir_type"] = "expected"
                        graph.edges[(node, nn)]["edge_creation"] = "imputed3"

    def changeMirDirection(self, graph, gnode, newDir):
        if newDir == 'up':
            graph.nodes[gnode]["attr_dict"]['color'] = 'g'
            graph.nodes[gnode]["attr_dict"]['node_expr_detection'] = 'imputed4'
            graph.nodes[gnode]["attr_dict"]["node_expr_direction"] = newDir
        elif newDir == "down":
            graph.nodes[gnode]["attr_dict"]['color'] = 'r'
            graph.nodes[gnode]["attr_dict"]['node_expr_detection'] = 'imputed4'
            graph.nodes[gnode]["attr_dict"]["node_expr_direction"] = newDir
        else:
            return

        nodeData = graph.node[gnode]["attr_dict"]
        nodeNN = [x for x in graph.neighbors(gnode)]

        for nn in nodeNN:

            mirEdge = (gnode, nn)
            ndata = graph.node[nn]["attr_dict"]

            childDir = graph.nodes[nn]["attr_dict"].get("node_expr_direction", "N/A")

            if childDir == newDir:
                if graph.edges[(gnode, nn)]["edge_type"] in self.priorityEdgeTypes:
                    continue
                graph.edges[mirEdge]["color"] = "#D50000"
                graph.edges[mirEdge]["linestyle"] = "dotted"
                graph.edges[mirEdge]["edge_type"] = "imputed"
                graph.edges[mirEdge]["edge_creation"] = "imputed4"
                graph.edges[mirEdge]["dir_type"] = "unexpected"

            elif childDir in self.dirToOpposite and self.dirToOpposite[childDir] == newDir:
                if graph.edges[(gnode, nn)]["edge_type"] in self.priorityEdgeTypes:
                    continue
                graph.edges[mirEdge]["color"] = "g"
                graph.edges[mirEdge]["linestyle"] = "dashed"
                graph.edges[mirEdge]["edge_type"] = "imputed"
                graph.edges[mirEdge]["edge_creation"] = "imputed4"
                graph.edges[mirEdge]["dir_type"] = "expected"

            else:
                print("Warning: undetermined directions")
                graph.edges[mirEdge]["color"] = "#546E7A"
                graph.edges[mirEdge]["linestyle"] = "dashed"

                if childDir in ["N/A", "unchanged"]:
                    graph.edges[mirEdge]["edge_type"] = "unknown"
                else:
                    graph.edges[mirEdge]["edge_type"] = "unexplained"

                graph.edges[mirEdge]["edge_creation"] = "imputed4"

    def colorEdges(self, graph):

        for edge in graph.edges():

            src = edge[0]
            tgt = edge[1]

            srcData = graph.node[src]["attr_dict"]
            tgtData = graph.node[tgt]["attr_dict"]

            srcDir = srcData.get("node_expr_direction", "N/A")
            tgtDir = tgtData.get("node_expr_direction", "N/A")

            srcDataOrigin = srcData.get("node_expr_detection", "N/A")
            tgtDataOrigin = tgtData.get("node_expr_detection", "N/A")

            srcType = srcData.get("type", None)
            tgtType = tgtData.get("type", None)

            edgeType = graph.edges[edge]["edge_type"]


            # here we want to color edges with data, but not imputed to highlight correct/incorrect edges
            if srcType == "mirna":
                if srcDataOrigin == "data":
                    edgeType = "imputed"
            else:
                if tgtType == "mirna":
                    if tgtDataOrigin == "data":
                        edgeType = "imputed"


            if edgeType == "expected":
                graph.edges[edge]["color"] = "#4CAF50"
                graph.edges[edge]["linestyle"] = "solid"
                graph.edges[edge]["dir_type"] = "expected"

            elif edgeType == "unexpected":

                graph.edges[edge]["color"] = "#F44336"
                graph.edges[edge]["linestyle"] = "solid"
                graph.edges[edge]["dir_type"] = "unexpected"

            elif edgeType == "imputed":

                if srcDir == self.dirToOpposite[tgtDir]:
                    graph.edges[edge]["color"] = "#4CAF50"  # green
                    graph.edges[edge]["linestyle"] = "dashed"
                    graph.edges[edge]["dir_type"] = "expected_imputed"

                elif srcDir == tgtDir and srcDir != "N/A":

                    graph.edges[edge]["color"] = "#F44336"  # red
                    graph.edges[edge]["linestyle"] = "dashed"
                    graph.edges[edge]["dir_type"] = "unexpected_imputed"

                else:
                    graph.edges[edge]["color"] = "#607D8B"  # gray
                    graph.edges[edge]["linestyle"] = "dashed"
                    graph.edges[edge]["dir_type"] = "unknown_imputed"


            else:
                graph.edges[edge]["color"] = "#607D8B" # gray
                graph.edges[edge]["linestyle"] = "dotted"
                graph.edges[edge]["dir_type"] = "unknown"

    def __impute4(self, graph):

        """

        imputes all remaining miRNA, where useful. this could actually be parametrized

        len(nnWithNoEvidence) == 0 // previous round could have provided evidences for neighbours // set miR direction such that fewest inconsistencies are added

        len(noEvidenceDir) == 1 // there are neighbours with no evidence, but all in one direction => set miR direction such that these can be explained
        minimze unexplained relations over inconsistencies, but this hopefully coincidences with minHarm ... most times

        downCount > upCount // upCount > downCount => change in appropriate direction

        downCount == upCount => no change

        :param graph:
        :return: None
        """

        mirNodes = 0

        for node in graph.nodes():

            nodeData = graph.node[node]["attr_dict"]

            if node in ["miR-181"]:
                print(node)

            if nodeData["type"] != "mirna":
                continue

            if nodeData.get("node_expr_direction", None) not in [None, "N/A", "unknown"]:
                continue

            mirNodes += 1

            nodeNN = [x for x in graph.neighbors(node)]

            if len(nodeNN) == 0:
                print("Singlet Node", node)
                continue

            nnRegs = {}
            nnEdgeTypes = Counter()
            nnDirDirections = Counter()
            nnWithEvidence = set()
            nnWithNoEvidence = set()

            dir2Singletons = set()

            for nn in nodeNN:

                nndata = graph.node[nn]["attr_dict"]
                childDir = graph.nodes[nn]["attr_dict"].get("node_expr_direction", "N/A")

                nnRegs[nn] = childDir
                nnDirDirections[childDir] += 1

                edgeType = graph.edges[(node, nn)].get("edge_type", "N/A")
                edgeDir = graph.edges[(node, nn)].get("dir_type", "unexpected").split("_")[0]

                nnEdgeTypes[edgeType] += 1

                if edgeDir in ["expected"] or self.nnHasEvidence(graph, nn):
                    nnWithEvidence.add((nn, nnRegs[nn]))
                else:
                    nnWithNoEvidence.add((nn, nnRegs[nn]))

            noEvidenceDir = set([x[1] for x in nnWithNoEvidence])
            print("impute4", node, nnDirDirections, noEvidenceDir, nnWithNoEvidence, nnWithEvidence)

            ndirDirs = {}
            for x in nnDirDirections:
                ndirDirs[x] = nnDirDirections[x]

            for x in nnWithNoEvidence:
                ndirDirs[x[1]] -= 1

            if len(nnWithNoEvidence) == 0:
                minHarm = sorted([(x, nnDirDirections[x]) for x in nnDirDirections], key=lambda x: x[1], reverse=True)
                newMirDir = self.dirToOpposite[minHarm[0][0]]

                print("impute4", node,"no missing evidence rule", dir2count, ndirDirs, newMirDir)

            elif len(noEvidenceDir) == 1:

                oppDir = list(noEvidenceDir)[0]
                newMirDir = self.dirToOpposite[oppDir]

                minHarm = sorted([(x, nnDirDirections[x]) for x in nnDirDirections], key=lambda x: x[1], reverse=True)


                print("impute4", node,"evidence dir rule", newMirDir)
                print("impute4", node,"evidence dir minharm", minHarm)
            else:

                dir2count = Counter()
                for x in nnWithNoEvidence:
                    dir2count[x[1]] += 1

                upCount = dir2count["up"]
                downCount = dir2count["down"]

                if downCount > upCount:
                    print("impute4", node,"lesser rule up", dir2count)
                    newMirDir = "up"

                elif upCount > downCount:
                    print("impute4", node,"lesser rule down", dir2count)
                    newMirDir = "down"


                else:

                    allNN = set(nnWithNoEvidence).union(nnWithEvidence)
                    allNNDirs = Counter()
                    for x in allNN:
                        allNNDirs[x[1]] += 1

                    allNNDirsBest = allNNDirs.most_common(1)[0]
                    minHarm = sorted([(x, dir2count[x]) for x in dir2count], key=lambda x: x[1], reverse=True)

                    minHarmEq = False
                    if len(minHarm) > 1:
                        if minHarm[0][1] == minHarm[1][1]:
                            minHarmEq = True

                    print("impute4", node,"nndirbest minharm test", allNNDirsBest, minHarm, minHarmEq)

                    if minHarmEq or allNNDirsBest == minHarm[0][0]:
                        newMirDir = self.dirToOpposite[minHarm[0][0]]

                        print("impute4", node,"nndirbest minharm rule", newMirDir, dir2count, allNNDirs)

                    else:

                        newMirDir = "N/A"
                        print("impute4", node,"no rule", newMirDir, dir2count, allNNDirs)



            self.changeMirDirection(graph, node, newMirDir)

    def checkConsistency(self, graph):
        for node in sorted([x for x in graph.nodes()]):

            nodeData = graph.node[node]["attr_dict"]

            if nodeData["type"] != "mirna":
                continue

            nodeNN = [x for x in graph.neighbors(node)]

            mirDirection = nodeData["node_expr_direction"]

            inconsistCount = 0
            for nn in nodeNN:
                nnData = graph.node[nn]["attr_dict"]

                nnDirection = nnData["node_expr_direction"]
                edgeType = graph.edges[(node, nn)]["edge_type"]
                edgeCreation = graph.edges[(node, nn)]["edge_creation"]

                if mirDirection == nnDirection:
                    print(node, nn, "inconsistent", edgeType)
                    inconsistCount += 1

            print(node, inconsistCount, len(nodeNN))

        nodeConditionCount = 0
        totalNodeCount = 0

        for node in sorted([x for x in graph.nodes()]):

            nodeData = graph.node[node]["attr_dict"]

            if nodeData["type"] == "mirna":
                continue

            totalNodeCount += 1
            nodeNN = [x for x in graph.neighbors(node)]

            mirDirection = nodeData["node_expr_direction"]

            inconsistCount = 0
            for nn in nodeNN:
                nnData = graph.node[nn]["attr_dict"]

                nnDirection = nnData["node_expr_direction"]
                edgeType = graph.edges[(node, nn)]["edge_type"]
                edgeCreation = graph.edges[(node, nn)]["edge_creation"]

                if mirDirection == nnDirection:
                    inconsistCount += 1

            if inconsistCount == len(nodeNN):
                nodeConditionCount += 1
                print(node, inconsistCount, len(nodeNN))

        print(nodeConditionCount, totalNodeCount)

    def countAttributeEdges(self, graph, attrNames):
        tCounter = Counter()
        for edge in graph.edges():
            edgeData = graph.edges[edge]

            tCounter[tuple([edgeData.get(x, "") for x in attrNames])] += 1

        print("Counts for attributes: ", ", ".join(attrNames))

        totalTCounter = 0
        for etype in tCounter:
            etypeCount = tCounter[etype]
            totalTCounter += etypeCount
            print(etype, etypeCount)

        print("Total", totalTCounter)
        print()

        return tCounter

    def countAttributeNodes(self, graph, attrNames):

        tCounter = Counter()
        for node in graph.nodes():
            nodeData = graph.node[node]["attr_dict"]

            tCounter[tuple([nodeData.get(x, "") for x in attrNames])] += 1

        print("Counts for attributes: ", ", ".join(attrNames))

        totalTCounter = 0
        for etype in tCounter:
            etypeCount = tCounter[etype]
            totalTCounter += etypeCount
            print(etype, etypeCount)

        print("Total", totalTCounter)
        print()

        return tCounter

global sigThreshold

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-d', '--detable', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True, help="output base")

    parser.add_argument('-p', '--pval', type=float, required=False, default=0.05, help="output base")

    parser.add_argument('--organisms', nargs='+', type=str, default=[], required=False)
    parser.add_argument('--disease', nargs='+', type=str, default=[], required=False)
    parser.add_argument('--go', nargs='+', type=str, default=[], required=False)
    parser.add_argument('--cells', nargs='+', type=str, default=[], required=False)
    parser.add_argument('--ncits', nargs='+', type=str, default=[], required=False)

    args = parser.parse_args()

    sigThreshold = args.pval

    contextDict = defaultdict(list)

    for x in args.organisms:
        contextDict['organism'].append(x)
    for x in args.disease:
        contextDict['disease'].append(x)
    for x in args.go:
        contextDict['go'].append(x)
    for x in args.cells:
        contextDict['cells'].append(x)
    for x in args.ncits:
        contextDict['ncits'].append(x)


    """
    LOAD GENE -> MIRNA HOST GENE
    """
    genes2mirFilename = os.path.dirname(__file__) + "/gene_mir_overlap.tsv"
    genename2mirs = defaultdict(set)

    for line in open(genes2mirFilename):
        aline = line.strip().split()

        aline[2] = aline[2].replace("hsa-", "")

        genename2mirs[aline[1].upper()].add(aline[2])



    """
    FETCH/BUILD INTERACTIONS
    """
    #contextDict = DataBaseAccessor.checkContext(contextDict)

    for didx, detable in enumerate(args.detable):

        defilename = detable.name

        indf = DataFrame.parseFromFile(defilename, skipChar='#', replacements={
            "None": None,
            "": None,
            "NA": None
        })

        geneSymCol = None
        inHeaders = indf.getHeader()

        deMIRs = defaultdict(lambda: (0, 1, ()))
        deGenes = {}

        allMIRs = defaultdict(lambda: (0, 1, ()))
        allGenes = {}

        if "gene_symbol" in inHeaders:
            geneSymCol = "gene_symbol"
        elif "Geneid" in inHeaders:
            geneSymCol = "Geneid"
        elif "id" in inHeaders:
            geneSymCol = "id"

        for ridx, row in enumerate(indf):

            try:
                robL2FC = float(row["ROB_log2FC"])
                robAdjPVal = float(row["ROB_ADJ.PVAL"])
            except:
                continue

            # if robAdjPVal > 0.05:
            #    continue

            geneSymbol = row[geneSymCol]

            #if geneSymbol.upper().startswith("MIR") and isNumber(geneSymbol[3:]):
            if geneSymbol.upper() in genename2mirs:

                assocMirs = genename2mirs[geneSymbol.upper()]

                for amir in assocMirs:

                    geneSymbol = amir#geneSymbol.split("-")[0].upper()

                    newFC, newPV, evs = allMIRs[geneSymbol]

                    if abs(newFC) > abs(robL2FC) and robL2FC != 0:
                        newFC = robL2FC
                        newPV = robAdjPVal

                    #newFC = oldL2FC if abs(oldL2FC) < abs(robL2FC) else robL2FC
                    #newPV = min(oldPVAL, robAdjPVal)

                    evs = set(evs)
                    evs.add(row[geneSymCol].upper())
                    evs = tuple(evs)

                    allMIRs[geneSymbol] = (newFC, newPV, evs)

                    if robAdjPVal < args.pval:

                        newFC, newPV, evs = deMIRs[geneSymbol]

                        if abs(newFC) > abs(robL2FC) and robL2FC != 0:
                            newFC = robL2FC
                            newPV = robAdjPVal

                        evs = set(evs)
                        evs.add(row[geneSymCol].upper())
                        evs = tuple(evs)

                        deMIRs[geneSymbol] = (newFC, newPV, evs)


            if True:
                allGenes[geneSymbol.upper()] = (robL2FC, robAdjPVal)

                if robAdjPVal < args.pval:
                    deGenes[geneSymbol.upper()] = (robL2FC, robAdjPVal)

        print("DE MIRNAS")
        for x in deMIRs:
            print(x, deMIRs[x])

        print("DE GENES", len(deGenes))

        REFETCHDATA = False

        mirnaContextDict = {"mirna": [x for x in deMIRs]}
        for x in contextDict:
            mirnaContextDict[x] = contextDict[x]

        geneContextDict = {"gene": [x for x in deGenes]}
        for x in contextDict:
            geneContextDict[x] = contextDict[x]

        if not REFETCHDATA and os.path.exists(args.output.name + ".mirnafetch.pickle"):
            with open(args.output.name + '.mirnafetch.pickle', 'rb') as f:
                mirnaHits = pickle.load(f)
        else:
            print("Fetching mirnas")
            print(mirnaContextDict)

            mirnaHits = DataBaseAccessor.fetch_mirna_interactions(
                mirnaContextDict,
                MIRNASTRPARTS=[miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]
            )

            with open(args.output.name + '.mirnafetch.pickle', 'wb') as f:
                pickle.dump(mirnaHits, f)
                print("Wrote out mirnafetch")

        if not REFETCHDATA and os.path.exists(args.output.name + ".genefetch.pickle"):

            with open(args.output.name + ".genefetch.pickle", 'rb') as f:
                geneHits = pickle.load(f)
        else:
            print("Fetching genes")
            print(geneContextDict)

            geneHits = DataBaseAccessor.fetch_mirna_interactions(
                geneContextDict,
                MIRNASTRPARTS=[miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]
            )

            with open(args.output.name + '.genefetch.pickle', 'wb') as f:
                pickle.dump(geneHits, f)
                print("Wrote out genefetch")

        print("All Fetched")

        print("Gene Hits", len(geneHits))
        print("miRNA Hits", len(mirnaHits))

        allLogFC = [deMIRs[x][0] for x in deMIRs] + [deGenes[x][0] for x in deGenes]
        minLogFC = min(allLogFC)
        maxLogFC = max(allLogFC)


        mgG = miRGeneGraph({
            "minLogFC": minLogFC,
            "maxLogFC": maxLogFC,
        })

        print("Create Graph")
        mgGraph = mgG.createGraph(mirnaHits, geneHits, genename2mirs)

        print("Impute Graph")
        mgG.imputeGraph(mgGraph)

        mgG.saveGraph(args.output, "graph_" + os.path.basename(defilename) + ".html", mgGraph, len(deGenes), len(allGenes))






