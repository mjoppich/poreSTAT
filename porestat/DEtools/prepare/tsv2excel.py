import csv
import sys
from xlsxwriter.workbook import Workbook

# Add some command-line logic to read the file names.

for tsv_file in sys.argv[1:]:

    #tsv_file = sys.argv[1]
    print(tsv_file)

    ain = tsv_file.split(".")
    ain[-1] = "xlsx"
    xlsx_file = ".".join(ain)


    # Create an XlsxWriter workbook object and add a worksheet.
    workbook = Workbook(xlsx_file)
    worksheet = workbook.add_worksheet()

    # Create a TSV file reader.
    tsv_reader = csv.reader(open(tsv_file,'rt'),delimiter="\t")

    # Read the row data from the TSV file and write it to the XLSX file.
    for row, data in enumerate(tsv_reader):

        ndata = []
        for x in data:
            try:
                y = float(x)
                ndata.append(y)
            except:
                ndata.append(x)

        worksheet.write_row(row, 0, ndata)

    # Close the XLSX file.
    workbook.close()