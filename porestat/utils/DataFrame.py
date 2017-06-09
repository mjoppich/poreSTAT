__author__ = 'joppich'

from .Files import readLines, fileExists
from .Numbers import toNumber, isNumber

import operator
from copy import deepcopy

class DataFrameException(Exception):

    def __init__(self, msg):
        super(DataFrameException, self).__init__()

        self.msg = msg

    def __str__(self):
        return self.msg

class DataRowException(Exception):

    def __init__(self, msg):
        super(DataRowException, self).__init__()

        self.msg = msg

class DataRow:

    def __init__(self, elementTuple, dHeader):

        self.dHeader = dHeader
        self.elements = elementTuple

    def __getitem__(self, item):

        if type(item) == int:

            if item < 0 or item > len(self.elements):
                raise DataRowException("Invalid column number: " + str(item))

            return self.elements[item]

        else:

            if not item in self.dHeader:
                raise DataRowException("Invalid column id: " + str(item))

            return self.elements[ self.dHeader[item] ]


class DataFrame:
    def __init__(self, oDefaultEmptyValue=None):

        self.dHeader = {}
        self.oDefaultEmptyValue = oDefaultEmptyValue
        self.vElements = []

    def getHeader(self):

        vReturn = [''] * len(self.dHeader)

        for x in self.dHeader:
            iIdx = self.dHeader[x]
            vReturn[iIdx] = x

        return vReturn

    def addColumns(self, names, default=None):

        for x in names:
            if x in self.dHeader:
                raise DataFrameException("Column already exists: " + x)

        for x in names:
            self.addColumn(x, default)


    def addColumn(self, name, default=None):
        """

        :param name: name of new column (must not be used already)
        :param default: default value to be added to all existing rows
        :return: -1 if not succeeded (name already used), column index otherwise
        """

        if name in self.dHeader:
            return -1

        newColIdx = len(self.dHeader)
        self.dHeader[name] = newColIdx

        def addcol(x):

            newx = []
            for elem in x:
                newx.append(elem)

            newx.append(default)

            return tuple(newx)

        self.applyToRow(addcol)

        return newColIdx

    def removeColumn(self, name):

        if not name in self.dHeader:
            return -1

        colIdx = self.dHeader[name]

        def removeCol(x):

            newx = []
            for i in range(0, len(x)):
                if i != colIdx:
                    newx.append(x[i])

            return tuple(newx)

        self.applyToRow(removeCol)

    def getIdxHeader(self):

        return deepcopy(self.dHeader)

    def toRowDict(self, oTuple):

        if len(oTuple) != len(self.dHeader):
            return None

        dReturn = {}

        for x in self.dHeader:
            dReturn[x] = oTuple[self.dHeader[x]]

        return dReturn

    def __getitem__(self, key):

        if type(key) is tuple:

            if (len(key) != 2):
                return None

            oElem = self.vElements[key[0]]

            try:
                return oElem[self.getColumnIndex(key[1])]

            except KeyError:

                return None

        else:

            try:
                return self.vElements[key]
            except KeyError:
                return None

    def __setitem__(self, key, value):

        if type(key) is tuple:

            if len(key) != 2:
                raise ValueError("invalid length of key")

            iRow = key[0]
            iCol = self.getColumnIndex(key[1])

            vElemLine = list(self.vElements[iRow])

            vElemLine[iCol] = value

            oLineTuple = tuple(vElemLine)

            self.vElements[iRow] = oLineTuple

        else:

            if len(value) != len(self.dHeader):
                raise ValueError("invalid length of value")

            self.vElements[key] = value

    def getColumnIndex(self, oColumn):

        if oColumn in self.dHeader:

            return self.dHeader[oColumn]

        else:

            try:
                return int(oColumn)
            except:
                raise DataFrameException("Invalid column: " + str(oColumn) + "\n\n Available columns: " + str([x for x in self.dHeader]))

    def applyByRow(self, oColumn, oFunc):

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.vElements)):
            vElemLine = list(self.vElements[i])

            self[i, iColumnIndex] = oFunc(vElemLine)

    def applyToRow(self, oFunc):
        """

        :param oFunc: must return tuple/list of length of header
        :return:
        """

        for i in range(0, len(self.vElements)):
            vElemLine = list(self.vElements[i])
            vReturnLine = oFunc(vElemLine)

            self[i] = vReturnLine

    def getColumn(self, oColumn):

        vReturn = []

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.vElements)):
            vReturn.append(self[i, iColumnIndex])

        return vReturn

    def setColumn(self, oColumn, vNewValues):

        if len(vNewValues) != len(self.vElements):
            raise ValueError("elements have different lengths")

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.vElements)):
            self.vElements[i, iColumnIndex] = vNewValues[i]


    def getRow(self, oColumn, oValue, oDefaultValue = None):

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.vElements)):

            if self.vElements[i][iColumnIndex] == oValue:
                return DataRow(self.vElements[i], self.dHeader)

        return oDefaultValue

    def getRows(self, oColumn, vValues):

        iColumnIndex = self.getColumnIndex(oColumn)

        vReturn = []
        for i in range(0, len(self.vElements)):

            if self.vElements[i][iColumnIndex] in vValues:
                vReturn.append( DataRow(self.vElements[i], self.dHeader) )

        return vReturn


    def findRow(self, oColumn, oValue, oDefaultValue=None):

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.vElements)):

            if self.vElements[i][iColumnIndex] == oValue:
                return self.vElements[i]

        return oDefaultValue

    def findRows(self, oColumn, vValues):

        iColumnIndex = self.getColumnIndex(oColumn)

        vReturn = []

        for i in range(0, len(self.vElements)):

            if self.vElements[i][iColumnIndex] in vValues:
                vReturn.append(self.vElements[i])

        return vReturn

    def __str__(self):

        sortedHeader = sorted(self.dHeader.items(), key=operator.itemgetter(1))

        vHeader = [str(x[0]) for x in sortedHeader]

        sStr = "\t".join(vHeader)

        for oLine in self.vElements:
            sStr += "\n"
            sStr += "\t".join([str(x) for x in oLine])

        return sStr

    @classmethod
    def createHeader(cls, oHeaderFrom, cDelim):

        aHeaderLine = None

        if type(oHeaderFrom) == str:
            aHeaderLine = oHeaderFrom.strip().split(cDelim)

        if type(oHeaderFrom) == tuple or type(oHeaderFrom) == list:
            aHeaderLine = oHeaderFrom

        if aHeaderLine == None:
            raise Exception("invalid headerfrom type: must be str, tuple or list")

        dDesc2Idx = {}
        dIdx2Desc = {}

        iIndex = 0
        # oHeader is header
        for sElem in aHeaderLine:
            dDesc2Idx[sElem] = iIndex
            dIdx2Desc[iIndex] = sElem

            iIndex += 1

        return (dDesc2Idx, dIdx2Desc)

    @classmethod
    def createLineTuple(cls, sLine, oHeader, cDelim, vNumberHeader=None, oEmptyValue=None):

        sLine = sLine.strip()

        aLine = sLine.split(cDelim)

        vLine = list(aLine)

        if vNumberHeader is None:
            vNumberHeader = [False] * len(oHeader)

        for e in range(0, len(vNumberHeader)):
            if vNumberHeader[e]:
                oRes = toNumber(vLine[e], vLine[e])
                vLine[e] = oRes

                # if oRes == None:
                #    print("some error " + str(aLine))

        while len(vLine) < len(oHeader):
            vLine.append(oEmptyValue)

        return tuple(vLine)

    @classmethod
    def parseFromFile(cls, sFileName, oHeader=None, cDelim='\t', bConvertTextToNumber=True, encoding="utf-8"):

        if not fileExists(sFileName):
            print("error loading file: " + sFileName)
            print("file does not exist")

            return None

        oNewDataFrame = DataFrame()

        vLines = readLines(sFileName=sFileName, encoding=encoding)

        iStartLine = 0

        if oHeader is None:

            # retrieve header from first line
            iStartLine += 1

            oHeaderRes = cls.createHeader(vLines[0], cDelim)
            oHeader = oHeaderRes[0]

        else:

            oHeaderRes = cls.createHeader(oHeader, cDelim)
            oHeader = oHeaderRes[0]

        oNewDataFrame.dHeader = oHeader

        vNumberHeader = [False] * len(oHeader)

        for i in range(iStartLine, len(vLines)):

            sLine = vLines[i]

            if i == iStartLine and bConvertTextToNumber:

                sLine = vLines[i].strip()
                aLine = sLine.split(cDelim)
                vLine = list(aLine)

                for e in range(0, len(vLine)):
                    if isNumber(vLine[e]):
                        vNumberHeader[e] = True

            lineTuple = cls.createLineTuple(sLine, oHeader, cDelim, vNumberHeader, oNewDataFrame.oDefaultEmptyValue)

            oNewDataFrame.vElements.append(lineTuple)

        return oNewDataFrame