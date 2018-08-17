import argparse
import os
import shutil
from enum import Enum
import matplotlib.pyplot as plt
import mpld3
from mpld3 import plugins

from porestat.utils.DataFrame import DataFrame, ExportTYPE
from ..utils.ArgParseExt import FileStubType


class PlotStyle(Enum):
    DEFAULT='default',
    XKCD='xkcd',
    FIVETHIRTYEIGHT='fivethirtyeight',
    CLASSIC='classic',
    BMH='bmh',
    DARK='dark_background'


class PlotStyleAction(argparse.Action):

    def __init__(self, option_strings, dest, nargs=None, const=None, default=None, type=None, choices=None, required=False, help=None, metavar=None):
        super(PlotStyleAction, self).__init__(option_strings, dest, nargs, const, default, type, choices, required, help, metavar)

        self.help = 'Sets type of file export and must be one of {m}'.format(m=', '.join([str(x.value) for x in PlotStyle]))

    def __call__(self, parser, args, values, option_string=None):

        try:
            eVal = PlotStyle[values]
            args.__dict__[ self.dest ] = eVal

        except:

            raise argparse.ArgumentError(None, 'PlotStyle can not be {n}, '
                                               'it must be one of {m}'.format(n=values,
                                                                              m=', '.join(
                                                                                  [str(x.name) for x in PlotStyle])))


class PlotSaveTypeAction(argparse.Action):

    def __init__(self, option_strings, dest, nargs=None, const=None, default=None, type=None, choices=None, required=False, help=None, metavar=None):
        super(PlotSaveTypeAction, self).__init__(option_strings, dest, nargs, const, default, type, choices, required, help, metavar)

        self.help = 'Sets type of file export and must be one of {m}'.format(m=', '.join([str(x.value) for x in PlotSaveTYPE]))

    def __call__(self, parser, args, values, option_string=None):

        if not len(values) == 1:
            raise argparse.ArgumentError(None, "ExportTYPE must be one value")


        value = values[0]

        if not type(value) == PlotSaveTYPE:

            try:
                eVal = PlotSaveTYPE[values.upper()]
                args.__dict__[ self.dest ] = eVal

            except:

                raise argparse.ArgumentError(None, 'ExportTYPE can not be {n}, '
                                                   'it must be one of {m}'.format(n=values,
                                                                                  m=', '.join([str(x.value) for x in PlotSaveTYPE])))

        else:
            args.__dict__[self.dest] = value


class PlotSaveTYPE(Enum):
    PNG='png'
    HTML='HTML'
    HTML_STRING='HTML_STRING'
    JSON='json'
    D3='D3'
    MPL='mpl'

    @classmethod
    def getFileExtension(cls, outType):

        if outType == PlotSaveTYPE.PNG:
            return 'png'

        if outType == PlotSaveTYPE.HTML:
            return 'html'

        if outType == PlotSaveTYPE.JSON:
            return 'json'


class PlotConfig:

    @classmethod
    def fromParserArgs(cls, simArgs):

        pltCfg = PlotConfig()

        if 'save_plot' in simArgs.__dict__ and simArgs.save_plot != None:
            pltCfg.saveToFile( simArgs.save_plot )

        if 'plot_style' in simArgs.__dict__ and simArgs.plot_style != None:
            pltCfg.setStyle( simArgs.plot_style )

        if 'save_plot_type' in simArgs.__dict__ and simArgs.save_plot_type != None:
            pltCfg.setOutputType(simArgs.save_plot_type)

        return pltCfg



    @classmethod
    def addParserArgs(cls, parser):

        parser.add_argument('-sp', '--save-plot', type=FileStubType('w'), nargs=1, help='path-prefix to file where plots are saved. Final file will be save-plot.xxx.png')
        parser.add_argument('-spt', '--save-plot-type', action=PlotSaveTypeAction, default=PlotSaveTYPE.PNG, type=PlotSaveTYPE, nargs=1)
        parser.add_argument('-ps', '--plot-style', action=PlotStyleAction, default=PlotStyle.BMH, type=PlotStyle, nargs=1)

        return parser

    def __init__(self):

        self.save_to_file = False
        self.save_file = None
        self.saved_plot = 0

        self.outputType = PlotSaveTYPE.PNG
        self.style = PlotStyle.DEFAULT
        self.transparent_bg = True

        self.d3js = None
        self.mpld3js=None

        self.d3js = mpld3.getD3js()
        self.mpld3js = mpld3.getmpld3js(True)

        self.createdPlots = []

        self.centeredStyle = {
                    'display': "block",
                    'margin': "auto",
                    'height': "750px",
                    'width': "50vw"
                }

    def __str__(self):

        allElements = [self.save_to_file, self.save_file, self.saved_plot, self.style, self.transparent_bg, self.createdPlots]
        allStrElements = [str(x) for x in allElements]

        return "[" + ", ".join(allStrElements) + "]"

    def saveToFile(self, filePath):

        self.save_to_file = True
        self.save_file = filePath

        if type(self.save_file) == list:
            assert(len(self.save_file) > 0)
            self.save_file = self.save_file[0]

        if self.save_file == None:
            print("pltcfg save file is empty", filePath, type(filePath))

        plt.ioff()

    def setOutputType(self, outType):
        self.outputType = outType

        if outType == PlotSaveTYPE.MPL:
            self.save_to_file = False
            plt.ion()

        if self.usesMPLD3():
            plt.ioff()


    def setStyle(self, style):
        self.style = style

    def startPlot(self):

        if self.style == PlotStyle.XKCD:
            plt.xkcd()
        else:
            plt.style.use( self.style.value )

    def getCreatedPlots(self):

        return self.createdPlots

    def resetCreatedPlots(self):
        self.createdPlots = []

    def usesMPLD3(self):

        if self.outputType in [ PlotSaveTYPE.HTML_STRING, PlotSaveTYPE.HTML, PlotSaveTYPE.JSON, PlotSaveTYPE.D3 ]:
            return True

        return False

    def addHTMLPlot(self, html):

        if self.outputType == PlotSaveTYPE.HTML_STRING:

            self.createdPlots.append(html)


    def makeTable(self, df):

        assert(type(df) == DataFrame)

        if self.outputType == PlotSaveTYPE.HTML_STRING:

            (header, body) = df.export(None, exType=ExportTYPE.HTML_STRING, html_element_id="dftable_" + str(len(self.createdPlots)))

            self.createdPlots.append(header + body)

            return

        else:
            df.export(None)

        if self.save_to_file != None and self.save_to_file and self.save_file != None:

            fileExt = PlotSaveTYPE.getFileExtension(self.outputType)

            if fileExt == None:
                print("Could not make table: invalid outputType", self.outputType)
                return

            exactFilename = self.save_file + ".%02d." + fileExt
            exactFilename = exactFilename % self.saved_plot
            self.saved_plot += 1

            if self.outputType == PlotSaveTYPE.HTML:
                df.export(exactFilename, exType=ExportTYPE.HTML)

    def makePlot(self, noTightLayout=False, figHeight='100%', figWidth='100%', centeredStyle=None):

        current_figure = plt.gcf()

        if self.usesMPLD3():
            plugins.connect(current_figure, plugins.MousePosition(fontsize=14))


        if self.outputType == PlotSaveTYPE.HTML_STRING:


            if centeredStyle == None:
                centeredStyle = self.centeredStyle



            outString = mpld3.fig_to_html(current_figure,
                                          template_type='notebook',
                                          d3_url=self.d3js,
                                          mpld3_url=self.mpld3js,
                                          figHeight=figHeight,
                                          figWidth=figWidth,
                                          styles=centeredStyle
                                          )
            self.createdPlots.append(outString)
            plt.close(current_figure)
            return

        if self.outputType == PlotSaveTYPE.D3:
            mpld3.show(current_figure)

        if self.save_to_file:

            exactFilename = self.save_file + ".%02d." + PlotSaveTYPE.getFileExtension( self.outputType )
            exactFilename = exactFilename % self.saved_plot
            self.saved_plot += 1

            if self.outputType == PlotSaveTYPE.HTML:
                mpld3.save_html(current_figure, exactFilename, template_type='simple', d3_url=self.d3js, mpld3_url=self.mpld3js)
            elif self.outputType == PlotSaveTYPE.JSON:
                mpld3.save_json(current_figure, exactFilename, d3_url=self.d3js, mpld3_url=self.mpld3js)
            elif self.outputType == PlotSaveTYPE.PNG:
                plt.savefig(exactFilename, transparent=self.transparent_bg, bbox_inches='tight')

            self.createdPlots.append(exactFilename)

        else: # if self.outputType == PlotSaveTYPE.MPL

            if not noTightLayout:

                legends = current_figure.legends

                makeTightLayout = True
                for lgd in legends + [ax.legend_ for ax in current_figure.axes]:

                    if lgd == None:
                        continue

                    lgd.draggable()

                    makeTightLayout = makeTightLayout and lgd._bbox_to_anchor == None


                if makeTightLayout:
                    plt.tight_layout()

            plt.show()

        plt.close(current_figure)


    def prepareHTMLOutput(self, folderPath, fileName, relativeImport=False):


        filenames = os.path.splitext(fileName)

        baseFileName = ".".join(filenames[0:len(filenames)-1])
        filesPath = os.path.join(folderPath, baseFileName)

        os.makedirs(folderPath, exist_ok=True)
        os.makedirs(filesPath, exist_ok=True)



        d3js = mpld3.getD3js()
        mpld3js = mpld3.getmpld3js(True)
        mathjaxJS = os.path.join(str(os.path.dirname(os.path.realpath(__file__))), "../data/MathJax.js")
        mathjaxSVGJS = os.path.join(str(os.path.dirname(os.path.realpath(__file__))), "../data/MathJaxSVG.js")

        outD3JS = os.path.join(filesPath, os.path.basename(d3js))
        outMPLD3JS = os.path.join(filesPath, os.path.basename(mpld3js))
        outMathJaxJS = os.path.join(filesPath, os.path.basename(mathjaxJS))
        outMathJaxSVGJS = os.path.join(filesPath, os.path.basename(mathjaxSVGJS))

        shutil.copyfile(d3js, outD3JS)
        shutil.copyfile(mpld3js, outMPLD3JS)
        shutil.copyfile(mathjaxJS, outMathJaxJS)
        shutil.copyfile(mathjaxSVGJS, outMathJaxSVGJS)

        if relativeImport:
            outD3JS = "./"+os.path.join(baseFileName, os.path.basename(d3js))
            outMPLD3JS = "./"+os.path.join(baseFileName, os.path.basename(mpld3js))
            outMathJaxJS = "./"+os.path.join(baseFileName, os.path.basename(mathjaxJS))
            outMathJaxSVGJS = "./"+os.path.join(baseFileName, os.path.basename(mathjaxSVGJS))


        with open(os.path.join(folderPath, fileName), 'w') as htmlFile:

            htmlFile.write(
                """
                <html>
                <head>
                    
                    <style>
                    * {
                        font-family: "Arial", Verdana, sans-serif;
                        }
                    </style>
                        
                    <script type="text/x-mathjax-config">
                          MathJax.Hub.Config({
                                extensions: ["tex2jax.js"],
                                jax: ["input/TeX", "output/HTML-CSS"],
                                tex2jax: {
                                    inlineMath: [ ['$','$'] ],
                                    displayMath: [ ['$$','$$'], ["\\[","\\]"] ],
                                    processEscapes: true
                                },
                                "HTML-CSS": { fonts: ["TeX"] },
                                TeX: {
                                Macros: {
                                    mathdefault: ["{#1}",1]
                                }
                          }
                          });
                    </script>   
                """)

            htmlFile.write("""     
                <script type="text/javascript" src="{d3jsURL}"></script>
                <script type="text/javascript" src="{mpld3jsURL}"></script>
                <script type="text/javascript" src="{mathjaxURL}"></script>
                
                </head>
                <body>
                """.format(mpld3jsURL=outMPLD3JS, d3jsURL=outD3JS, mathjaxURL=outMathJaxJS, mathjaxSVGURL=outMathJaxSVGJS))

            for x in self.createdPlots:
                htmlFile.write(x + "\n")

                htmlFile.flush()

            htmlFile.write( "<script type=\"text/javascript\" src=\"{mathjaxSVGURL}\"></script>".format(mathjaxSVGURL=outMathJaxSVGJS))
            htmlFile.write("</body></html>\n")

