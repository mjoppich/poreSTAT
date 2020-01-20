import argparse
import os,sys
from collections import Counter, defaultdict

import networkx
import requests
import json


sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../")
from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
from porestat.DEtools.miRNAUtils import miRNA, miRNAPART, isNumber
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

        serverAddress = "https://turingwww.bio.ifi.lmu.de"
        serverPort = None
        serverPath = "yancDB"

        r = requests.post(cls.makeServerAddress(cls.serverAddress, cls.serverPort, cls.serverPath, "find_interactions"),
                          data=json.dumps(requestDict))

        jsonRes = r.json()

        return jsonRes

    @classmethod
    def fetch_mirna_interactions(cls, requestDict, MIRNASTRPARTS=[miRNAPART.MATURE, miRNAPART.ID]):

        if not "sentences" in requestDict:
            requestDict["sentences"] = "FALSE"

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

            edge = (source, target)

            foundInteractions.add(edge)

        return foundInteractions


class miRGeneGraph:

    def logFC2Size(self, logFC, minS=20, maxS=40):

        logFCr = maxLogFC + (-minLogFC)

        logFCx = logFC + -minLogFC

        logFCx = logFCx / logFCr

        return minS + logFCx * (maxS - minS)

    def exprColor(self, logfc):

        if not str(logfc).isnumeric():

            if logfc == "up":
                logfc = 1
            elif logfc == "down":
                logfc = -1

        if logfc < 0:
            return "#FF0000"

        if logfc == 'up' or logfc > 0:
            return "#00FF00"

        return "#555555"


    def findDataForNode(self, nodeName):
        un = nodeName.upper()

        if "MIR" in un or "LET" in un:
            nodeData = None
            objMir = miRNA(nodeName)

            for mir in allMIRs:
                if objMir.accept(mir):
                    nodeData = allMIRs.get(mir, None)
                    # print("Matching", nodeName, "with", mir)
                    break

            nodeType = "mirna"
            nodeShape = "triangle"

        else:
            nodeData = allGenes.get(un, None)
            nodeType = "gene"
            nodeShape = "square"

        if nodeData == None:
            return {"type": nodeType, "shape": nodeShape, "border_style": "dashed", "log2FC": 0, "adjPval": 1, "de_measured": "false", }, False

        isMeasured = True if nodeData != (0.0, 1.0) else False

        returnDict = {
            "de_measured": str(isMeasured).lower(),
            "border_style": "dashed" if nodeData[1] > 0.05 else "solid",
            "type": nodeType,
            "log2FC": nodeData[0],
            "adjPval": nodeData[1],
            "shape": nodeShape,
            "color": self.exprColor(nodeData[0]),
            "size": self.logFC2Size(nodeData[0]),
            "node_expr_detection": "data" if nodeData != (0.0, 1.0) else "imputed0",
        }

        if isMeasured:
            returnDict["node_expr_direction"] = "down" if nodeData[0] < 0 else "up"

        return returnDict, True




    def __init__(self, mirnaHits, geneHits, props):
        self.priorityEdgeTypes = ["expected", "unexpected"]
        self.dirToOpposite = {"up": "down", "down": "up"}

        self.gp = GraphPlot()

        self.minLogFC = props["minLogFC"] if "minLogFC" in props else -5
        self.maxLogFC = props["maxLogFC"] if "maxLogFC" in props else 5


        graph = networkx.Graph()
        deNodes = set()

        for edge in mirnaHits:

            srcData, srcDE = self.findDataForNode(edge[0])
            tgtData, tgtDE = self.findDataForNode(edge[1])

            if srcDE:
                deNodes.add(edge[0])
            if tgtDE:
                deNodes.add(edge[1])

            graph.add_node(edge[0], attr_dict=srcData)
            graph.add_node(edge[1], attr_dict=tgtData)
            graph.add_edge(edge[0], edge[1])

            graph.edges[edge]['edge_creation'] = "original_mirhit"

        for edge in geneHits:

            srcData, srcDE = self.findDataForNode(edge[0])
            tgtData, tgtDE = self.findDataForNode(edge[1])

            if srcDE:
                deNodes.add(edge[0])
            if tgtDE:
                deNodes.add(edge[1])

            graph.add_node(edge[0], attr_dict=srcData)
            graph.add_node(edge[1], attr_dict=tgtData)

            graph.add_edge(edge[0], edge[1])

            graph.edges[edge]['edge_creation'] = "original_genehit"


        # let's remove genes which have not been measured at all ...
        delNode = []
        for node in graph.nodes():

            nn = [x for x in networkx.all_neighbors(graph, node)]
            if len(nn) == 1:
                if not node in deNodes:
                    delNode.append(node)

        for x in delNode:
            graph.remove_node(x)


        self.colorInitial(graph)
        self.__impute1(graph)
        self.__impute1_1(graph)
        self.__impute2(graph)
        self.__impute3(graph)
        self.__impute4(graph)

        self.gp.showGraph(graph, location="/mnt/d/dev/data/human_plaque/tmp/")




    def colorInitial(self, graph):

        for edge in graph.edges():

            srcNode = edge[0]
            tgtNode = edge[1]

            srcNodeData = graph.nodes[srcNode]["attr_dict"]
            tgtNodeData = graph.nodes[tgtNode]["attr_dict"]

            graph.edges[edge]["linestyle"] = "solid"
            graph.edges[edge]["edge_type"] = "unknown"

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

                elif (reg0 > 0 and reg1 > 0) or (reg0 < 0 and reg1 < 0):
                    graph.edges[edge]["color"] = "r"
                    graph.edges[edge]["edge_type"] = "unexpected"


    def __impute1(self, graph):
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

                if 'up' in mirDir:
                    exprVal = -1
                    graph.nodes[node]["attr_dict"]["node_expr_direction"] = "down"
                else:
                    exprVal = 1
                    graph.nodes[node]["attr_dict"]["node_expr_direction"] = "up"

                graph.nodes[node]["attr_dict"]["color"] = self.exprColor(exprVal)
                graph.nodes[node]["attr_dict"]['node_expr_detection'] = 'imputed1'
                print(node, "imputed")

                for nn in nodeNN:

                    if graph.edges[(node, nn)]["edge_type"] in self.priorityEdgeTypes:
                        continue

                    graph.edges[(node, nn)]["color"] = "g"
                    graph.edges[(node, nn)]["linestyle"] = "dashed"
                    graph.edges[(node, nn)]["edge_type"] = "imputed"
                    graph.edges[(node, nn)]["edge_creation"] = "imputed1"


    def __impute1_1(self, graph):
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

            if node == "FGF10":
                print(node, nodeNN, nnDirs, acceptNode)

            if acceptNode:
                geneDir = self.dirToOpposite[list(nnDirs)[0]]

                graph.nodes[node]["attr_dict"]["node_expr_direction"] = geneDir
                graph.nodes[node]["attr_dict"]["node_expr_detection"] = "imputed1.1"
                graph.nodes[node]["attr_dict"]["color"] = self.exprColor(geneDir)

                for nn in nodeNN:

                    if graph.edges[(node, nn)]["edge_type"] in self.priorityEdgeTypes:
                        continue

                    graph.edges[(node, nn)]["color"] = "g"
                    graph.edges[(node, nn)]["linestyle"] = "dashed"
                    graph.edges[(node, nn)]["edge_type"] = "imputed"
                    graph.edges[(node, nn)]["edge_creation"] = "imputed1.1"

                    print(node, nn, "imputed1.1")

        # impute regulations


    def nnHasEvidence(self, graph, node):
        nodeData = graph.node[node]["attr_dict"]
        nodeNN = [x for x in graph.neighbors(node)]
        for nn in nodeNN:

            nndata = graph.node[nn]["attr_dict"]
            edgeType = graph.edges[(node, nn)]["edge_type"]

            if edgeType in ['imputed', 'expected']:
                return True

        return False


    def __impute2(self, graph):
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

                edgeType = graph.edges[(node, nn)]["edge_type"]
                nnEdgeTypes[edgeType] += 1

                if edgeType in ["imputed", "expected"] or self.nnHasEvidence(graph, nn):
                    nnWithEvidence.add((nn, nnRegs[nn]))
                else:
                    nnWithNoEvidence.add((nn, nnRegs[nn]))

            regCounter = Counter()
            for gene in nodeNN:
                regCounter[nnRegs[gene]] += 1

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

                if sameDir:

                    for e in nnWithEvidence:
                        undetEdges.add((node, e[0]))

                    for e in nnWithNoEvidence:
                        newImputEdges.add((node, e[0]))

                    undetEdges = list(undetEdges)
                    imputeDir = imputeDirs[0]
                    imputeDir = self.dirToOpposite[imputeDir]

                    print(node, sameDir, imputeDir, ncounter)

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
                            graph.edges[(node, nn)]["color"] = "#ffb733"
                            graph.edges[(node, nn)]["linestyle"] = "dotted"
                            graph.edges[(node, nn)]["edge_type"] = "imputed"
                            graph.edges[(node, nn)]["edge_creation"] = "imputed2"

                        elif imputeDir == self.dirToOpposite[nnRegs[nn]]:
                            graph.edges[(node, nn)]["color"] = "g"
                            graph.edges[(node, nn)]["linestyle"] = "dashed"
                            graph.edges[(node, nn)]["edge_type"] = "imputed"
                            graph.edges[(node, nn)]["edge_creation"] = "imputed2"
                        else:
                            print("error, no direction", imputeDir, nnRegs[nn], node, nn)

                if not sameDir:
                    # if graph.nodes[node]["attr_dict"]["de_measured"] == "false":
                    #    graph.nodes[node]["attr_dict"]['color'] = '#FF8F00'
                    # else:
                    #    graph.nodes[node]["attr_dict"]['color'] = exprColor(graph.nodes[node]["attr_dict"]["log2FC"])
                    print(node, regCounter, nnEdgeTypes, nnWithEvidence, sameDir, undetEdges)


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
        for node in graph.nodes():

            nodeData = graph.node[node]["attr_dict"]

            if nodeData["type"] != "mirna":
                continue

            if nodeData.get("node_expr_direction", None) not in [None, "unknown"]:
                continue


            nodeNN = [x for x in graph.neighbors(node)]

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

                if edgeType in ["imputed", "expected"] or self.nnHasEvidence(graph, nn):
                    nnWithEvidence.add((nn, nnRegs[nn]))
                else:
                    nnWithNoEvidence.add((nn, nnRegs[nn]))

            if len(nnWithNoEvidence) == 0:
                print("can be explained:", node, nnDirDirections, )

                childDir2Effect = defaultdict(list)

                for nn in nodeNN:

                    allEvs = self.nnGetEvidences(graph,nn)
                    nndata = graph.node[nn]["attr_dict"]

                    if nndata["log2FC"] < 0:
                        childDir = "down"
                    elif nndata["log2FC"] > 0:
                        childDir = "up"
                    else:
                        childDir = "N/A"

                    averageEffect = abs(nndata["log2FC"]) / len(allEvs) if len(allEvs) > 0 else 0
                    childDir2Effect[childDir].append((nn, averageEffect))
                    print(node, nn, childDir, allEvs, averageEffect)

                etToEffect = [(x, sum([y[1] for y in childDir2Effect[x]])) for x in childDir2Effect]

                etToEffect = sorted(etToEffect, key=lambda x: x[1], reverse=True)

                finalTarget = etToEffect[0]

                targetDir = finalTarget[0]

                if targetDir == "up":
                    targetDir = "down"
                elif targetDir == "down":
                    targetDir = "up"

                print(node, finalTarget, targetDir)

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
                        graph.edges[(node, nn)]["color"] = "#ffb733"
                        graph.edges[(node, nn)]["linestyle"] = "dotted"
                        graph.edges[(node, nn)]["edge_type"] = "imputed"
                        graph.edges[(node, nn)]["edge_creation"] = "imputed3"

                    elif childDir == self.dirToOpposite[targetDir]:
                        if graph.edges[(node, nn)]["edge_type"] in self.priorityEdgeTypes:
                            continue
                        graph.edges[(node, nn)]["color"] = "g"
                        graph.edges[(node, nn)]["linestyle"] = "dashed"
                        graph.edges[(node, nn)]["edge_type"] = "imputed"
                        graph.edges[(node, nn)]["edge_creation"] = "imputed3"


    def changeMirDirection(self, graph, gnode, newDir):
        if newDir == 'up':
            graph.nodes[gnode]["attr_dict"]['color'] = 'g'
            graph.nodes[gnode]["attr_dict"]['node_expr_detection'] = 'imputed4'
            graph.nodes[gnode]["attr_dict"]["node_expr_direction"] = newDir
        else:
            graph.nodes[gnode]["attr_dict"]['color'] = 'r'
            graph.nodes[gnode]["attr_dict"]['node_expr_detection'] = 'imputed4'
            graph.nodes[gnode]["attr_dict"]["node_expr_direction"] = newDir

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

            elif childDir in self.dirToOpposite and self.dirToOpposite[childDir] == newDir:
                if graph.edges[(gnode, nn)]["edge_type"] in self.priorityEdgeTypes:
                    continue
                graph.edges[mirEdge]["color"] = "g"
                graph.edges[mirEdge]["linestyle"] = "dashed"
                graph.edges[mirEdge]["edge_type"] = "imputed"
                graph.edges[mirEdge]["edge_creation"] = "imputed4"

            else:
                print("Warning: undetermined directions")
                graph.edges[mirEdge]["color"] = "#546E7A"
                graph.edges[mirEdge]["linestyle"] = "dashed"
                graph.edges[mirEdge]["edge_type"] = "errored"
                graph.edges[mirEdge]["edge_creation"] = "imputed4"


    def __impute4(self, graph):

        mirNodes = 0

        for node in graph.nodes():

            nodeData = graph.node[node]["attr_dict"]

            if nodeData["type"] != "mirna":
                continue

            if nodeData.get("node_expr_direction", None) not in [None, "unknown"]:
                continue

            mirNodes += 1

            nodeNN = [x for x in graph.neighbors(node)]

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

                edgeType = graph.edges[(node, nn)]["edge_type"]
                nnEdgeTypes[edgeType] += 1

                if edgeType in ["imputed", "expected"] or self.nnHasEvidence(graph, nn):
                    nnWithEvidence.add((nn, nnRegs[nn]))
                else:
                    nnWithNoEvidence.add((nn, nnRegs[nn]))

            noEvidenceDir = set([x[1] for x in nnWithNoEvidence])
            print(node, nnDirDirections, noEvidenceDir, nnWithNoEvidence)

            if len(nnWithNoEvidence) == 0:
                minHarm = sorted([(x, nnDirDirections[x]) for x in nnDirDirections], key=lambda x: x[1], reverse=True)
                newMirDir = self.dirToOpposite[minHarm[0][0]]

                print("no missing evidence rule", dir2count, ndirDirs, newMirDir)

            elif len(noEvidenceDir) == 1:

                oppDir = list(noEvidenceDir)[0]
                newMirDir = self.dirToOpposite[oppDir]

                print("evidence dir rule", newMirDir)
            else:

                dir2count = Counter()
                for x in nnWithNoEvidence:
                    dir2count[x[1]] += 1

                if dir2count["up"] * 0.5 > dir2count["down"]:
                    print("lesser rule up")
                    newMirDir = "down"
                elif dir2count["down"] * 0.5 > dir2count["up"]:
                    print("lesser rule down")
                    newMirDir = "up"

                else:
                    ndirDirs = {}
                    for x in nnDirDirections:
                        ndirDirs[x] = nnDirDirections[x]

                    for x in nnWithNoEvidence:
                        ndirDirs[x[1]] -= 1

                    minHarm = sorted([(x, dir2count[x]) for x in dir2count], key=lambda x: x[1], reverse=True)
                    newMirDir = self.dirToOpposite[minHarm[0][0]]

                    print("no rule", dir2count, ndirDirs, newMirDir)

            self.changeMirDirection(graph, node, newMirDir)


    def countAttributeEdge(self, graph, attrNames):
        tCounter = Counter()
        for edge in graph.edges():
            edgeData = graph.edges[edge]

            tCounter[tuple([edgeData.get(x, None) for x in attrNames])] += 1

        print("Counts for attributes: ", ", ".join(attrNames))

        totalTCounter = 0
        for etype in tCounter:
            etypeCount = tCounter[etype]
            totalTCounter += etypeCount
            print(etype, etypeCount)

        print("Total", totalTCounter)
        print()

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

    def countAttributeNodes(self, graph, attrNames):

        tCounter = Counter()
        for node in graph.nodes():
            nodeData = graph.node[node]["attr_dict"]

            tCounter[tuple([nodeData[x] for x in attrNames])] += 1

        print("Counts for attributes: ", ", ".join(attrNames))

        totalTCounter = 0
        for etype in tCounter:
            etypeCount = tCounter[etype]
            totalTCounter += etypeCount
            print(etype, etypeCount)

        print("Total", totalTCounter)
        print()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-d', '--detable', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-o', '--output', nargs='+', type=str, required=False, help="output base")


    parser.add_argument('--organisms', nargs='+', type=str, default=[], required=False)
    parser.add_argument('--disease', nargs='+', type=str, default=[], required=False)
    parser.add_argument('--go', nargs='+', type=str, default=[], required=False)
    parser.add_argument('--cells', nargs='+', type=str, default=[], required=False)
    parser.add_argument('--ncits', nargs='+', type=str, default=[], required=False)

    args = parser.parse_args()

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


    contextDict = DataBaseAccessor.checkContext(contextDict)

    defilename = "/mnt/d/dev/data/human_plaque/reports/stable_vs_unstable.star.msEmpiRe_DESeq2.tsv"

    indf = DataFrame.parseFromFile(defilename, skipChar='#', replacements={
        "None": None,
        "": None,
        "NA": None
    })

    geneSymCol = None
    inHeaders = indf.getHeader()

    deMIRs = {}
    deGenes = {}

    allMIRs = {}
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

        if geneSymbol.upper().startswith("MIR") and isNumber(geneSymbol[3:]):

            geneSymbol = geneSymbol.split("-")[0]

            if geneSymbol.upper() in allMIRs:

                oldL2FC, oldPVAL = allMIRs[geneSymbol.upper()]

                newFC = oldL2FC if abs(oldL2FC) < abs(robL2FC) else robL2FC
                newPV = min(oldPVAL, robAdjPVal)

                allMIRs[geneSymbol.upper()] = (newFC, newPV)
            else:
                allMIRs[geneSymbol.upper()] = (robL2FC, robAdjPVal)

            if robAdjPVal < 0.05:
                deMIRs[geneSymbol.upper()] = (robL2FC, robAdjPVal)

        else:
            allGenes[geneSymbol.upper()] = (robL2FC, robAdjPVal)

            if robAdjPVal < 0.05:
                deGenes[geneSymbol.upper()] = (robL2FC, robAdjPVal)

    print("DE MIRNAS", deMIRs)
    print("DE GENES", len(deGenes))

    mirnaContextDict = {"mirna": [x for x in deMIRs]}
    for x in contextDict:
        mirnaContextDict[x] = contextDict[x]

    geneContextDict = {"gene": [x for x in deGenes]}
    for x in contextDict:
        geneContextDict[x] = contextDict[x]

    mirnaHits = DataBaseAccessor.fetch_mirna_interactions(
        mirnaContextDict
    )

    geneHits = DataBaseAccessor.fetch_mirna_interactions(
        geneContextDict
    )

    print("DE MIRNAS", deMIRs, len(mirnaHits))
    print("DE GENES", len(deGenes), len(geneHits))

    allLogFC = [deMIRs[x][0] for x in deMIRs] + [deGenes[x][0] for x in deGenes]
    minLogFC = min(allLogFC)
    maxLogFC = max(allLogFC)

    mgG = miRGeneGraph(mirnaHits, geneHits, {
        "minLogFC": minLogFC,
        "maxLogFC": maxLogFC,
    })


