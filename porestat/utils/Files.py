import os
import sys


def makePath(path):
    """

    :param path: path to be normalized
    :return: normalizes path by adding path sep to end if needed (/tmp/bla -> /tmp/bla/)
    """

    if os.path.isdir(path):

        if path[len(path)-1] == os.path.sep:
            return path

        return path + os.path.sep

    elif os.path.isfile(path):

        dirname = os.path.dirname(path)

        if path[len(path)-1] == os.path.sep:
            return path

        return path + os.path.sep
    else:
        ValueError("make path can not determine type of: " + str(path))

def printToFile(vVector, sFileName, cInElementSep = None, cElementSep = '\n'):
    """

    :param vVector:
    :param cInElementSep:
    :param cElementSep:
    :return:
    """

    oScoringFile = open(sFileName, 'w')

    iLines = 0
    for elem in vVector:

        if iLines > 0:

            oScoringFile.write( cElementSep )

        if cInElementSep == None:

            oScoringFile.write(str(elem))

        else:

            for i in range(0, len(elem)):

                if (i > 0):
                    oScoringFile.write(cInElementSep)

                oScoringFile.write(str(elem[i]))

        iLines = iLines + 1


    oScoringFile.close()

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

def fileExists(sFileName):

    return os.path.exists(sFileName)


def readLines(sFileName, vLines = None, encoding = "utf-8", iSkipLines = 0):

    vReturn = []

    oFile = open(sFileName, 'r', encoding=encoding)

    iLineCount = 0

    for sLine in oFile:

        if iSkipLines > 0:
            iSkipLines -= 1
            continue

        if (vLines == None) or (iLineCount in vLines):
            vReturn.append(sLine)

            if not vLines is None:
                vLines.remove(iLineCount)

        if (vLines != None) and (len(vLines) == 0):
            break

        iLineCount = iLineCount + 1

    return vReturn


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def pathEmpty(sPathName):
    isDir = os.path.isdir(sPathName)

    isEmpty = False

    if isDir:
        isEmpty = len(os.listdir(sPathName)) == 0

    return isDir and isEmpty


def pathWritable(sPathName):

    isDir = os.path.isdir(sPathName)
    isWritable = False

    if isDir:
        isWritable = os.access(sPathName, os.W_OK)

    return isDir and isWritable