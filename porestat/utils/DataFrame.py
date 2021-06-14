import os
import shutil
from collections import defaultdict
from enum import Enum
import operator
from jinja2 import Template
from openpyxl import Workbook

__author__ = 'joppich'

from .Files import readLines, fileExists
from .Numbers import toNumber, isNumber
import random
import argparse

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

class ExportTYPEAction(argparse.Action):

    def __init__(self, option_strings, dest, nargs=None, const=None, default=None, type=None, choices=None, required=False, help=None, metavar=None):
        super(ExportTYPEAction, self).__init__(option_strings, dest, nargs, const, default, type, choices, required, help, metavar)

        self.help = 'Sets type of file export and must be one of {m}'.format(m=', '.join([str(x.value) for x in ExportTYPE]))

    def __call__(self, parser, args, values, option_string=None):

        try:
            eVal = ExportTYPE[values.upper()]
            args.__dict__[ self.dest ] = eVal

        except:

            raise argparse.ArgumentError(None, 'ExportTYPE can not be {n}, '
                                               'it must be one of {m}'.format(n=values,
                                                                              m=', '.join([str(x.value) for x in ExportTYPE])))


class ExportTYPE(Enum):
    CSV=0
    TSV=1
    XLSX=2
    HTML=3
    HTML_STRING=4
    LATEX=5
    TABLE_FILTER=6


class DataSeries:
    """

    A simple list for data

    """

    def __init__(self, data):
        self.data = data
        self.dataIdx = 0

    def __iter__(self):
        self.dataIdx = 0
        return self

    def __next__(self):
        try:
            item = self.__getitem__(self.dataIdx)
        except (IndexError, DataRowException):
            raise StopIteration()

        self.dataIdx += 1

        return item

    def __getitem__(self, item):
        return self.data[item]

    def __setitem__(self, key, value):

        if type(self.data) == tuple:
            oldData = list(self.data)
            oldData[key] = value
            self.data = tuple(oldData)
        else:
            self.data[key] = value

    def append_data(self, value):
        self.data = [x for x in self.data] + [value]

    def remove_data(self, idx):
        self.data = [ self.data[i] for i in range(0, len(self.data)) if i != idx ]

    def to_tuple(self):
        return tuple(self.data)

    def to_list(self):
        return list(self.data)

    def to_set(self):
        return set(self.data)

    def to_scalar(self):
        if len(self.data) > 1:
            raise DataFrameException("Data is not scalar: " + str(self.data))

        return self.data[0]

class DataColumnAccess:
    """
    A DataSeries with column names
    """

    def __init__(self, column2idx):


        self.column2idx = column2idx
        self.idx2default = []

    def getColumnIndex(self, oColumn):

        if oColumn in self.column2idx:
            return self.column2idx[oColumn]

        else:

            try:
                idx = int(oColumn)

                if idx < 0 or idx > len(self.column2idx):
                    raise DataRowException("Invalid column number: " + str(oColumn))

                return idx

            except:
                raise DataFrameException("Invalid column: " + str(oColumn) + "\n\n Available columns: " + str([x for x in self.column2idx]))




    def copyHeader(self):
        return deepcopy(self.column2idx)

    def getHeader(self):

        vReturn = [''] * len(self.column2idx)

        for x in self.column2idx:
            iIdx = self.column2idx[x]
            vReturn[iIdx] = x

        return vReturn

    def columnExists(self, name):
        return name in self.column2idx

    def addColumn(self, name):
        """

        :param name: name of new column (must not be used already)
        :param default: default value to be added to all existing rows
        :return: -1 if not succeeded (name already used), column index otherwise
        """

        if name in self.column2idx:
            raise DataFrameException("Column already exists: " + name)

        newColIdx = len(self.column2idx)
        self.column2idx[name] = newColIdx

        return newColIdx

    def removeColumn(self, name):

        colIdx = self.getColumnIndex(name)
        DataSeries.remove_data(self, colIdx)

        del self.column2idx[name]

        for x in self.column2idx:
            if self.column2idx[x] > colIdx:
                self.column2idx[x] = self.column2idx[x]-1



class DefaultColumns():

    def __init__(self, **kwargs):

        self.hasDefault = False

        if 'global_default' in kwargs:
            self.hasDefault = True
            self.idx2default = defaultdict(lambda: kwargs['global_default'])
        else:
            self.idx2default = {}

    def add_default(self, colIdx, default):
        self.idx2default[colIdx] = default

    def get_default(self, colIdx):

        if not colIdx in self.idx2default and not self.hasDefault:
            raise DataFrameException("No default value for column: " + str(colIdx))

        return self.idx2default[colIdx]

class DefaultDataColumnAccess(DefaultColumns, DataColumnAccess):

    def __init__(self, dHeader, **kwargs):

        DataColumnAccess.__init__(self, dHeader)
        DefaultColumns.__init__(self, **kwargs)

    def addColumn(self, name, default=None):

        newColIdx = DataColumnAccess.addColumn(self, name)
        DefaultColumns.add_default(self, newColIdx, default)

        return newColIdx


    def addColumns(self, names, ignoreDuplicates=False, default=None):


        if not ignoreDuplicates:
            for x in names:
                if self.columnExists(x):
                    raise DataFrameException("Column already exists: " + x)

        for x in names:
            if not self.columnExists(x):
                self.addColumn(x, default)


class DataRow(DefaultDataColumnAccess, DataSeries):

    def __init__(self, dHeader, data = [], **kwargs):

        DefaultDataColumnAccess.__init__(self, dHeader, **kwargs)
        DataSeries.__init__(self, data)

    def __getitem__(self, item):
        """

        :param item: the key to look for
        :return: if @item is column name, return column name. if @item is int, the item-th element is return.
        """
        try:
            idx = self.getColumnIndex(item)
            return self.data[idx]
        except:
            raise DataRowException("Column not found: " + str(item))

    def getIndex(self, idx):

        assert(type(idx) == int)

        return self.data[idx]

    def addColumn(self, name, default=None):
        DefaultDataColumnAccess.addColumn(self, name, default)
        DataSeries.append_data(self, default)

    def to_tuple(self):
        return DataSeries.to_tuple(self)

    def to_pairs(self):
        return [ (x, self.__getitem__(x)) for x in self.column2idx ]

    def get(self, colname, default=None):

        if self.columnExists(colname):
            return self[colname]

        return default


    def to_dict(self):
        ret = {}

        for x in self.column2idx:
            ret[x] = self.__getitem__(x)

        return ret

    @classmethod
    def fromDict(cls, dictionary):

        allitems = list(dictionary.items())

        dHeader = {}
        velements = []

        for i in range(0, len(allitems)):
            dHeader[ allitems[i][0] ] = i
            velements.append(allitems[i][1])

        return DataRow(dHeader, data=velements)


class DataFrame(DataSeries, DefaultDataColumnAccess):
    """

        behave like a dataseries, but have columns

    """

    def __init__(self, default=None):

        DefaultDataColumnAccess.__init__(self, {}, global_default=default)# for a col in row
        DataSeries.__init__(self, []) # for each row

        self.title = None
        self.filepath = None

    def setTitle(self, value):

        self.title = value

    def setFilepath(self, value):

        self.filepath = value

    def __getitem__(self, item):

        if type(item) == int:
            return DataRow(self.column2idx, DataSeries.__getitem__(self, item))
        elif type(item) == str or (item in self.column2idx):

            seriesData = [None] * len(self.data)

            idx = self.getColumnIndex(item)
            for i in range(0, len(self.data)):
                val = self.data[i][idx]
                seriesData[i] = val

            return DataSeries(seriesData)

        elif type(item) == tuple or type(item) == list:

            if len(item) == 1:
                return DataRow(self.column2idx, DataSeries.__getitem__(self, item[0]))

            elif len(item) == 2:
                rowData = DataSeries.__getitem__(self, item[0])
                row = DataRow(self.column2idx, rowData)

                return row[item[1]]

        elif type(item) == slice:

            elemIdx = range(item.start, item.stop, item.step)
            return [DataRow(self.column2idx, DataSeries.__getitem__(self, x)) for x in elemIdx]

    def __setitem__(self, key, value):

        if type(key) == int:
            DataSeries.__setitem__(self, key, value)
        elif type(key) == tuple or type(key) == list:

            if len(key) == 1:
                DataSeries.__setitem__(self, key, value)

            elif len(key) == 2:
                rowData = DataSeries.__getitem__(self, key[0])
                row = DataRow(self.column2idx, rowData)
                row[key[1]] = value
                row_tuple = row.to_tuple()

                DataSeries.__setitem__(self, key[0], row_tuple)

        elif type(key) == slice:

            elemIdx = range(key.start, key.stop, key.step)

            if len(elemIdx) == len(value):
                for x in elemIdx:
                    DataSeries.__setitem__(x, value[x])

    def get_abs_largest(self, idcol, targetCols):

        id2val = defaultdict(0)
        for row in self:
            idVal = row[idcol]
            for col in targetCols:

                colVal = row[col]

                if abs(colVal) > abs(id2val[idVal]):
                    id2val[idVal] = colVal

        return sorted([(x, id2val[x]) for x in id2val], key=lambda x: x[1], reverse=True)




    def merge(self, other):

        if not type(other)==DataFrame:
            raise DataFrameException("Can only merge two dataframes")

        owncols = set(self.column2idx)
        othercols = set(other.column2idx)

        for x in owncols:
            if x not in othercols:
                raise DataFrameException("Missing col in other: " + str(x))

        for row in other:
            self.addRow( row )

    def addColumn(self, name, default=None):
        newColIdx = DefaultDataColumnAccess.addColumn(self, name, default)

        for i in range(0, len(self)):
            oldData = list(self.data[i])
            oldData.append(default)

            self.data[i] = tuple(oldData)

        return newColIdx


    def __len__(self):
        return len(self.data)

    def toRowDict(self, oTuple):

        if len(oTuple) != len(self.column2idx):
            return None

        dReturn = {}

        for x in self.column2idx:
            dReturn[x] = oTuple[self.column2idx[x]]

        return dReturn



    def namedRows(self, idxcol, datacols):


        if idxcol == None:
            idxcol = [x for x in range(0, len(self.data))]
        else:
            idxcol = [x[idxcol] for x in self.data]

        res = {}

        for datacolName in datacols:
            dataDict = {}

            datacol = datacols[datacolName]

            for i in range(0, len(self.data)):
                key = idxcol[i]
                value = self.data[i][datacol]

                dataDict[key] = value

            res[datacolName] = DataRow.fromDict(dataDict)

        return res


    def toDataRow(self, idxcol = None, datacol = None):

        if datacol == None:
            datacol = len(self.column2idx)-1

        if idxcol == None:
            idxcol = [x for x in range(0, len(self.data))]
        else:
            idxcol = [x[idxcol] for x in self.data]

        dataDict = {}

        for i in range(0, len(self.data)):

            key = idxcol[i]
            value = self.data[i][datacol]

            dataDict[key] = value

        return DataRow.fromDict(dataDict)


    def applyByRow(self, oColumn, oFunc):

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.data)):
            vElemLine = list(self.data[i])

            self[i, iColumnIndex] = oFunc(vElemLine)

    def applyToRow(self, oFunc):
        """

        :param oFunc: must return tuple/list of length of header
        :return:
        """

        for i in range(0, len(self.data)):
            vElemLine = list(self.data[i])
            vReturnLine = oFunc(vElemLine)
            self[i] = vReturnLine

    def getColumn(self, oColumn):

        vReturn = []

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.data)):
            vReturn.append(self[i, iColumnIndex])

        return vReturn

    def filterRows(self, filter):

        ndata = []

        for row in self.data:

            rowResult = filter(row)

            if rowResult:
                ndata.append(row)

        self.data = ndata

    def selectColumns(self, colDict):

        for x in colDict:
            if not x in self.column2idx:
                raise DataFrameException("Unknown column " + x + "\nAvailable columns: " + str([x for x in self.getHeader()]))

        outDF = DataFrame()
        outDF.addColumns([colDict[x] for x in colDict])

        for row in self:

            drDict = {}
            for x in colDict:
                drDict[colDict[x]] = row[x]

            outDF.addRow( DataRow.fromDict(drDict))

        return outDF

    def setColumn(self, oColumn, vNewValues):

        if len(vNewValues) != len(self.data):
            raise ValueError("elements have different lengths")

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.data)):
            self.data[i, iColumnIndex] = vNewValues[i]

    def addRow(self, row):

        if not (type(row) == DataRow):
            raise DataFrameException( 'Trying to insert invalid row type: ' + str(type(row)))

        rowheaders = set([x for x in row.column2idx])
        selfheaders = set([x for x in self.column2idx])
        setDifference = rowheaders.difference(selfheaders)

        if len(setDifference) > 0:
            raise DataFrameException( 'Row does not contain all needed headers: ' + str(selfheaders))


        rowDict = row.to_dict()

        newrow = [None] * len(self.column2idx)
        for x in self.column2idx:

            if x in rowDict:
                newrow[ self.column2idx[x] ] = rowDict[x]
            else:
                newrow[ self.column2idx[x]] = self.idx2default[self.column2idx[x]]

        self.data.append(tuple(newrow))


    def updateRowIndexed(self, idColumn, drows, ignoreMissingCols=False, addIfNotFound=False):

        # create index
        idColIdx = self.column2idx[idColumn]

        colIndex2Rows = defaultdict(set)

        for i in range(0, len(self.data)):
            rowT = self.data[i]
            colIndex2Rows[rowT[idColIdx]].add(i)

        addedRows = 0
        updatedRows=0

        addedRowGenes = []

        for drow in drows:

            drowDict = drow

            if type(drow) == DataRow:
                drowDict = drow.to_dict()

            findElement = drowDict[idColumn]

            if findElement in colIndex2Rows:

                rowIdxs = colIndex2Rows[findElement]


                for rowIdx in rowIdxs:
                    rowT = self.data[rowIdx]

                    rowL = list(rowT)

                    for col in drowDict:
                        colidx = self.getColumnIndex(col)
                        rowL[colidx] = drowDict[col]

                    rowL = tuple(rowL)
                    self.data[rowIdx] = rowL

                updatedRows += 1


            else:
                if addIfNotFound:
                    addedRowGenes.append(drowDict[idColumn])
                    addedRows += 1

                    if ignoreMissingCols:
                        if type(drow) == dict:
                            drow = DataRow.fromDict(drow)

                        rowheaders = set([x for x in drow.column2idx])

                        selfheaders = set([x for x in self.column2idx])
                        setDifference = rowheaders.difference(selfheaders)

                        for col in setDifference:
                            drow.addColumn(col, self.get_default(self.getColumnIndex(col)))

                        self.addRow(drow)

                    else:
                        self.addRow(drow)

        #print("Updated Rows", updatedRows)
        #print("Added Rows", addedRows)
        #print("Index Elems", random.sample([x for x in colIndex2Rows], min(len(colIndex2Rows), 10)))
        #print("Index Elems Added", random.sample(addedRowGenes, min(len(addedRowGenes), 10)))

    def updateRow(self, findID, idColumn, drow, addIfNotFound=False):

        idColIdx = self.column2idx[idColumn]

        drowDict = drow

        if type(drow) == DataRow:
            drowDict = drow.to_dict()

        found = False

        #print("length data", len(self.data))

        for i in range(0, len(self.data)):

            rowT = self.data[i]

            if rowT[idColIdx] == findID:

                rowL = list(rowT)

                for col in drowDict:
                    colidx = self.getColumnIndex(col)

                    rowL[colidx] = drowDict[col]

                self.data[i] = tuple(rowL)
                found=True


        if not found and addIfNotFound:
            self.addRow(drow)



    def getRow(self, oColumn, oValue, oDefaultValue = None):

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.data)):

            if self.data[i][iColumnIndex] == oValue:
                return DataRow(self.data[i], self.column2idx)

        return oDefaultValue

    def getRows(self, oColumn, vValues):

        iColumnIndex = self.getColumnIndex(oColumn)

        vReturn = []
        for i in range(0, len(self.data)):

            if self.data[i][iColumnIndex] in vValues:
                vReturn.append( DataRow(self.data[i], self.column2idx) )

        return vReturn


    def findRow(self, oColumn, oValue, oDefaultValue=None):

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.data)):

            if self.data[i][iColumnIndex] == oValue:
                return self[i]

        return oDefaultValue

    def findRows(self, oColumn, vValues):

        iColumnIndex = self.getColumnIndex(oColumn)

        vReturn = []

        for i in range(0, len(self.data)):

            if self.data[i][iColumnIndex] in vValues:
                vReturn.append(self.data[i])

        return vReturn

    def __str__(self):

        sortedHeader = sorted(self.column2idx.items(), key=operator.itemgetter(1))

        vHeader = [str(x[0]) for x in sortedHeader]

        sStr = "\t".join(vHeader)

        for oLine in self.data:
            sStr += "\n"
            sStr += "\t".join([str(x) for x in oLine])

        return sStr

    def _makeStr(self, sep='\t', include_header=True):

        sortedHeader = sorted(self.column2idx.items(), key=operator.itemgetter(1))
        vHeader = [str(x[0]) for x in sortedHeader]
        
        allData = []
        if include_header:
            allData.append(sep.join(vHeader))

        for oLine in self.data:
            allData.append(sep.join([str(x) for x in oLine]))

        return "\n".join(allData)

    def _writeToFile(self, content, filename):

        with open(filename, 'w') as file:
            file.write(content)

    def _makeXLSX(self, outFile):
        wb = Workbook()
        # grab the active worksheet
        ws = wb.active

        # Data can be assigned directly to cells
        sortedHeader = sorted(self.column2idx.items(), key=operator.itemgetter(1))
        vHeader = [str(x[0]) for x in sortedHeader]

        ws.append(vHeader)

        for oLine in self.data:
            ws.append([str(x) for x in oLine])

        # Save the file
        wb.save( outFile )

    def _makeLatex(self):

        columns = self.getHeader()

        latexTableRows = []

        latexTableRows.append("\\begin{tabular}{" + "l"*len(columns) + "}")
        latexTableRows.append(" & ".join(columns) + "\\\\")
        latexTableRows.append("\\midrule")

        for row in self:

            rowVals = [str(row[self.getColumnIndex(x)]) for x in columns]

            latexTableRows.append(" & ".join(rowVals) + "\\\\ \hline")


        latexTableRows.append("\\bottomrule")
        latexTableRows.append("\\end{tabular}")
        latexTableRows.append("")
        latexTableRows.append("")
        latexTableRows.append("")

        return "\n".join(latexTableRows)


    def _makeHTMLString(self, html_element_id=None):

        headpart = """
                <link rel="stylesheet" href="https://cdn.datatables.net/1.10.15/css/jquery.dataTables.min.css">
                <script src="https://code.jquery.com/jquery-1.12.4.js"></script>
                <script src="https://cdn.datatables.net/1.10.15/js/jquery.dataTables.min.js"></script>
        """

        bodypart = """
        {% if title %}
        {{title}}
        {% endif %}       
        <table id="{{html_element_id}}" class="display" cellspacing="0" width="100%">
                <thead>
                <tr>
                {% for column in columns %}
                    <th>{{column}}</th>
                    {% endfor %}
                </tr>
                </thead>

                <tbody>
                {%- for row in rows %}
                <tr>
                    {% for idx in indices %}
                    <td>{{ row[idx] }}</td>
                    {%- endfor -%}
                </tr>
                {% endfor -%}
                </tbody>

                <tfoot>
                <tr>
                {% for column in columns %}
                    <th>{{column}}</th>
                    {% endfor %}
                </tr>
                </tfoot>

                </table>

                <script>
                $(document).ready(function() {
                    $('#{{html_element_id}} tfoot th').each( function () {
                        var title = $(this).text();
                        $(this).html( '<input type="text" placeholder="Search '+title+'" />' );
                    } );

                    // DataTable
                    var table = $('#{{html_element_id}}').DataTable({
                                                        "columnDefs": [
                                                            { "type": "numeric-comma" }
                                                        ]
                                                    } );

                    table.columns().every( function (i) {
                        var that = this;

                        
                        $( 'input', this.footer() ).on( 'keyup change', function () {

                            var searchStr = this.value;

                            if ((searchStr.length >= 1) && ((searchStr[0] == '<') || (searchStr[0] == '>')))
                            {
                                if (searchStr.length > 1)
                                {
                                    var searchNum = searchStr.substring(1,searchStr.length);
                                    searchNum = +searchNum
        
                                    if (searchStr[0] == '<')
                                    {
                                    
        
                                        $.fn.dataTable.ext.search.push(
                                            function( settings, data, dataIndex ) {
    
                                                var rownum = +data[i];
                                        
                                                if (rownum < searchNum)
                                                {
                                                    return true;
                                                }
                                                
                                                return false;
                                            }
                                        );
                                    }
                                    else {
                                        $.fn.dataTable.ext.search.push(
                                            function( settings, data, dataIndex ) {
    
                                                var rownum = +data[i];
                                        
                                                if (rownum > searchNum)
                                                {
                                                    return true;
                                                }
                                                
                                                return false;
                                            }
                                        );
                                    }

                                    table.draw();
                                    
                                    $.fn.dataTable.ext.search.pop()
                                }




                            } else {
                                if ( that.search() !== this.value ) {
                                that.search( this.value )
                                    .draw();
                                }
                            }

                            
                        } );
                    } );

                } );
                </script>
        """

        sortedHeader = sorted(self.column2idx.items(), key=operator.itemgetter(1))
        vHeader = [str(x[0]) for x in sortedHeader]
        vIndices = [x[1] for x in sortedHeader]

        if html_element_id == None:
            html_element_id = "dftable"

        jinjaTemplate = Template(bodypart)
        output = jinjaTemplate.render(rows=self.data, indices=vIndices, columns=vHeader, title=self.title, html_element_id=html_element_id)

        return (headpart, output)


    def _makeTableFilter(self, outFile):

        (headpart, bodypart) = self._makeHTMLStringFilterTable("dftable")

        if self.title != None:
            bodypart = "<h1>"+self.title+"</h1>" + bodypart

        htmlfile="""

        <html>
            <head>
        """ + headpart + """
            </head>
            <body>
        """ + bodypart + """
            </body>
        </html>
        """

        with open(outFile, 'w') as outHtml:
            outHtml.write(htmlfile)

        def copyFolders(root_src_dir, root_target_dir):

            for src_dir, dirs, files in os.walk(root_src_dir):
                dst_dir = src_dir.replace(root_src_dir, root_target_dir)
                if not os.path.exists(dst_dir):
                    os.mkdir(dst_dir)
                for file_ in files:
                    src_file = os.path.join(src_dir, file_)
                    dst_file = os.path.join(dst_dir, file_)
                    if os.path.exists(dst_file):
                        os.remove(dst_file)

                    shutil.copy(src_file, dst_dir)


        sourceDir = os.path.dirname(__file__) + "/../data/tablefilter"
        targetDir = os.path.dirname(outFile) + "/tablefilter"

        print("copy tablefilter files from", sourceDir, "to", targetDir)
        copyFolders(sourceDir, targetDir)

    def _makeHTMLStringFilterTable(self, html_element_id=None):

        headpart = """
        """

        bodypart = """
        {% if title %}
        {{title}}
        {% endif %}
        
        <button id="csvButton" type="button">Save current table!</button>
        
        <table id="{{html_element_id}}" class="display" cellspacing="0" width="100%">
                <thead>
                <tr>
                {% for column in columns %}
                    <th>{{column}}</th>
                    {% endfor %}
                </tr>
                </thead>

                <tbody>
                {%- for row in rows %}
                <tr>
                    {% for idx in indices %}
                    <td>{{ row[idx] }}</td>
                    {%- endfor -%}
                </tr>
                {% endfor -%}
                </tbody>

                <tfoot>
                <tr>
                {% for column in columns %}
                    <th>{{column}}</th>
                    {% endfor %}
                </tr>
                </tfoot>

                </table>

<script src="tablefilter/tablefilter.js"></script>

<script data-config>
    var filtersConfig = {
        base_path: 'tablefilter/',
        alternate_rows: true,
        rows_counter: true,
        btn_reset: true,
        loader: true,
        status_bar: true,
        mark_active_columns: true,
        highlight_keywords: true,
        sticky_headers: true,
        col_types: [{{coltypes}}],
        custom_options: {
            cols:[],
            texts: [],
            values: [],
            sorts: []
        },
        col_widths: [],
        extensions:[{ name: 'sort' }]
    };

    var tf = new TableFilter("{{html_element_id}}", filtersConfig);
    tf.init();

function download_csv(csv, filename) {
    var csvFile;
    var downloadLink;

    // CSV FILE
    csvFile = new Blob([csv], {type: "text/csv"});

    // Download link
    downloadLink = document.createElement("a");

    // File name
    downloadLink.download = filename;

    // We have to create a link to the file
    downloadLink.href = window.URL.createObjectURL(csvFile);

    // Make sure that the link is not displayed
    downloadLink.style.display = "none";

    // Add the link to your DOM
    document.body.appendChild(downloadLink);

    // Lanzamos
    downloadLink.click();
}

function isHidden(el) {
    var style = window.getComputedStyle(el);
    return ((style.display === 'none') || (style.visibility === 'hidden'))
}

function export_table_to_csv(html, filename) {
	var csv = [];
	var rows = document.querySelectorAll("table tr");
	
    for (var i = 0; i < rows.length; i++) {
		var row = [], cols = rows[i].querySelectorAll("td, th");

        if (!isHidden(rows[i]))
        {
            for (var j = 0; j < cols.length; j++) 
            {
                colText = ""+cols[j].innerText;
                colText = colText.replace(/(\\r\\n|\\n|\\r)/gm, ';')
                row.push(colText);

            }

            if (row.length > 0)
            {
                csv.push(row.join("\\t"));
            }		

        }
		    
	}

    // Download CSV
    download_csv(csv.join("\\n"), filename);
}

document.addEventListener('readystatechange', event => {

    if (event.target.readyState === "interactive") {      //same as:  document.addEventListener("DOMContentLoaded"...   // same as  jQuery.ready
            console.log("Ready state");

        document.getElementById("csvButton").addEventListener("click", function () {
            var html = document.getElementById("{{html_element_id}}").outerHTML;
            export_table_to_csv(html, "table.tsv");
        });

    }

    if (event.target.readyState === "complete") {
        console.log("Now external resources are loaded too, like css,src etc... ");
        
        document.getElementById("csvButton").addEventListener("click", function () {
            var html = document.getElementById("{{html_element_id}}").outerHTML;
            export_table_to_csv(html, "table.tsv");
        });
    }

});

                </script>

        """

        column2types = defaultdict(set)

        for ridx, row in enumerate(self):

            if ridx > 100:
                break

            for column in self.column2idx:

                ctype = type(row[column])
                column2types[column].add(ctype)

        ntype = type(None)
        for x in column2types:

            print(x, column2types[x])
            if ntype in column2types[x]:
                column2types[x].remove(ntype)

            print(x, column2types[x])


        jsCols = []
        for column in self.getHeader():

            nonNumberType = len(column2types[column]) == 0

            if not nonNumberType:
                for ctype in column2types[column]:
                    if not ctype in [int, float, complex]:
                        nonNumberType = True

            if nonNumberType:
                jsCols.append("\"string\"")
            else:
                jsCols.append("\"number\"")



        sortedHeader = sorted(self.column2idx.items(), key=operator.itemgetter(1))
        vHeader = [str(x[0]) for x in sortedHeader]
        vIndices = [x[1] for x in sortedHeader]

        if html_element_id == None:
            html_element_id = "dftable"

        jinjaTemplate = Template(bodypart)
        output = jinjaTemplate.render(rows=self.data, indices=vIndices, columns=vHeader, title=self.title,
                                      html_element_id=html_element_id, coltypes=", ".join(jsCols))

        return (headpart, output)



    def _makeHTML(self, outFile):

        (headpart, bodypart) = self._makeHTMLString()

        if self.title != None:
            bodypart = "<h1>"+self.title+"</h1>" + bodypart

        htmlfile="""

        <html>
            <head>
        """ + headpart + """
            </head>
            <body>
        """ + bodypart + """
            </body>
        </html>
        """

        with open(outFile, 'w') as outHtml:
            outHtml.write(htmlfile)







    def export(self, outFile, exType=ExportTYPE.TSV, html_element_id=None, include_header=True):

        outputText = None

        if exType == ExportTYPE.XLSX and outFile != None:
            self._makeXLSX(outFile)
            return

        if exType == ExportTYPE.HTML and outFile != None:
            self._makeHTML(outFile)
            return

        if exType == ExportTYPE.HTML_STRING:
            return self._makeHTMLString(html_element_id)

        if exType == ExportTYPE.TSV:
            outputText = self._makeStr('\t', include_header=include_header)
        elif exType == ExportTYPE.CSV:
            outputText = self._makeStr(';', include_header=include_header)

        elif exType == ExportTYPE.LATEX:
            outputText = self._makeLatex()

        elif exType == None:
            outputText = self._makeStr('\t')

        elif exType == ExportTYPE.TABLE_FILTER:
            self._makeTableFilter(outFile)
            return

        if outFile == None:
            print(outputText)
            return
        else:
            self._writeToFile( outputText, outFile )
            return

        raise DataFrameException('Invalid export type {n} with output file {m}'.format(n=str(exType), m=str(outFile)))


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
    def createLineTuple(cls, sLine, oHeader, cDelim, vNumberHeader=None, oEmptyValue=None, replacements=None):

        sLine = sLine.strip()

        aLine = sLine.split(cDelim)

        vLine = list(aLine)

        origVH = vNumberHeader

        if vNumberHeader is None:
            vNumberHeader = [False] * len(oHeader)

        for e in range(0, len(vNumberHeader)):
            if e >= len(vLine):
                print("incorrect e")
                print(origVH)
                print(len(vLine), vLine)
                print(len(vNumberHeader), vNumberHeader)

            if vNumberHeader[e]:
                oRes = toNumber(vLine[e], vLine[e])
                vLine[e] = oRes

                # if oRes == None:
                #    print("some error " + str(aLine))


        while len(vLine) < len(oHeader):
            vLine.append(oEmptyValue)

        if replacements != None:
            for i in range(0, len(vLine)):
                if vLine[i] in replacements:
                    vLine[i] = replacements[vLine[i]]

        return tuple(vLine)

    @classmethod
    def parseFromFile(cls, sFileName, oHeader=None, skipLines=0, cDelim='\t', bConvertTextToNumber=True, encoding="utf-8", skipChar=None, replacements=None,header=True):


        if type(sFileName) == str:

            if not fileExists(sFileName):
                print("error loading file: " + sFileName)
                print("file does not exist")

                return None

            sFileName = open(sFileName, 'r')

        oNewDataFrame = DataFrame()

        vLines = readLines(sFileName=sFileName, encoding=encoding, iSkipLines=skipLines, skipChar=skipChar)

        if len(vLines) == 0:
            return None

        iStartLine = 0

        if oHeader is None:

            if header==True:
                # retrieve header from first line
                iStartLine += 1
                oHeaderRes = cls.createHeader(vLines[0], cDelim)
                oHeader = oHeaderRes[0]
            else:
                numCols = len(vLines[0].split(cDelim))
                oHeader = {}
                for i in range(0,numCols):
                    oHeader["var{}".format(i)] = i

        else:

            oHeaderRes = cls.createHeader(oHeader, cDelim)
            oHeader = oHeaderRes[0]

        oNewDataFrame.column2idx = oHeader

        vNumberHeader = [False] * len(oHeader)

        for i in range(iStartLine, len(vLines)):

            sLine = vLines[i]

            if i == iStartLine and bConvertTextToNumber:

                sLine = vLines[i].strip()
                aLine = sLine.split(cDelim)
                vLine = list(aLine)

                for e in range(0, len(vLine)):

                    if replacements != None:
                        if vLine[e] in replacements:
                            vLine[e] = replacements[vLine[e]]

                    if isNumber(vLine[e]):
                        vNumberHeader[e] = True

            assert(len(vNumberHeader) == len(oHeader))
            assert(len(vLine) == len(oHeader))

            lineTuple = cls.createLineTuple(sLine=sLine, oHeader=oHeader, cDelim=cDelim, vNumberHeader=vNumberHeader, oEmptyValue=None, replacements=replacements)

            oNewDataFrame.data.append(lineTuple)

        return oNewDataFrame