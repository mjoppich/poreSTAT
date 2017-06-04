import os
import sys
from itertools import chain
from collections import defaultdict, Counter

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def mergeDicts( dict1, dict2):
    dict3 = {}
    for k, v in chain(dict1.items(), dict2.items()):

        if k in dict3:

            if not type(v) == type(dict3[k]):
                raise Exception("You try to merge two different objects!")

            if type(v) == list:

                dict3[k] = dict3[k] + v

            elif type(v) == set:

                dict3[k] = dict3[k].union(v)

            elif type(v) == dict:

                dict3[k] = mergeDicts(dict3[k], v)

            elif type(v) == Counter:

                dict3[k] = mergeCounter(dict3[k], v)

            elif type(v) == int or type(v) == float:

                dict3[k] = dict3[k] + v

            else:

                retSet = set()
                retSet.add(v)
                retSet.add(dict3[k])

                if len(retSet) != 1:
                    dict3[k] = retSet
        else:

            dict3[k] = v

    return dict3

def fileExists(sFileName):

    return os.path.exists(sFileName)

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