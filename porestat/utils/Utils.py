import os
from collections import Counter
from collections import defaultdict
from itertools import chain


def mergeDicts( dict1, dict2, resultType=dict):
    dict3 = resultType()

    for k, v in chain(dict1.items(), dict2.items()):

        if k in dict3:

            if not type(v) == type(dict3[k]) and type(dict3[k] == list):

                dict3[k].append(v)
                continue

            if type(v) == list:

                dict3[k] = dict3[k] + v

            elif type(v) == set:

                dict3[k] = dict3[k].union(v)

            elif type(v) == dict:

                dict3[k] = mergeDicts(dict3[k], v)

            elif type(v) == Counter:

                dict3[k] = mergeCounter(dict3[k], v)

            elif type(v) == defaultdict:
                dict3[k] = mergeDefaultDict(dict3[k], v)

            elif type(v) == int or type(v) == float:

                dict3[k] = dict3[k] + v

            elif type(v) == tuple:
                tmp = dict3[k]
                dict3[k] = list()
                dict3[k].append(tmp)
                dict3[k].append(v)


            else:

                retSet = set()
                retSet.add(v)
                retSet.add(dict3[k])

                if len(retSet) != 1:
                    dict3[k] = retSet
        else:

            dict3[k] = v

    return dict3


def readLines(sFileName, encoding = 'utf-8'):
    oFile = open(sFileName, 'r', encoding=encoding)

    return oFile.readlines()

def mergeCounter( counter1, counter2):

    mergedCounter = Counter()

    for x in counter1:
        mergedCounter[x] = counter1[x]

    for x in counter2:
        mergedCounter[x] += counter2[x]

    return mergedCounter

def mergeDefaultDict( dict1, dict2):

    mergedDefaultDict = defaultdict(dict1.default_factory)

    for x in dict1:
        mergedDefaultDict[x] = dict1[x]

    for x in dict2:
        if x in mergedDefaultDict:
            mergedDefaultDict[x] = mergedDefaultDict[x] + dict2[x]
        else:
            mergedDefaultDict[x] = dict2[x]

    return mergedDefaultDict

def readFile(sFileName, sFunc, iSkip = 0, encoding = "utf-8", iMaxLines = -1):

    oFile = open(sFileName, 'r', encoding=encoding)

    iLineCount = 0

    for sLine in oFile:

        if iSkip > 0:
            iSkip = iSkip - 1
            continue

        if iMaxLines >= 0:
            if iLineCount >= iMaxLines:
                break

        sFunc(sLine, iLineCount)
        iLineCount = iLineCount + 1