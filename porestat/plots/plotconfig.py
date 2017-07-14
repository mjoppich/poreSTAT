import argparse
from enum import Enum
import matplotlib.pyplot as plt
import mpld3
from mpld3 import plugins


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

        try:
            eVal = PlotSaveTYPE[values.upper()]
            args.__dict__[ self.dest ] = eVal

        except:

            raise argparse.ArgumentError(None, 'ExportTYPE can not be {n}, '
                                               'it must be one of {m}'.format(n=values,
                                                                              m=', '.join([str(x.value) for x in PlotSaveTYPE])))
class PlotSaveTYPE(Enum):
    PNG='png'
    HTML='HTML'
    HTML_STRING='HTML_STRING'
    JSON='json'
    D3='D3'

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

        parser.add_argument('-sp', '--save-plot', type=str, help='path-prefix to file where plots are saved. Final file will be save-plot.xxx.png')
        parser.add_argument('-spt', '--save-plot-type', action=PlotSaveTypeAction, default=PlotSaveTYPE.PNG)
        parser.add_argument('-ps', '--plot-style', action=PlotStyleAction, default=PlotStyle.BMH)

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

    def __str__(self):

        allElements = [self.save_to_file, self.save_file, self.saved_plot, self.style, self.transparent_bg, self.createdPlots]
        allStrElements = [str(x) for x in allElements]

        return "[" + ", ".join(allStrElements) + "]"

    def saveToFile(self, filePath):

        self.save_to_file = True
        self.save_file = filePath

        plt.ioff()

    def setOutputType(self, outType):
        self.outputType = outType

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

    def makePlot(self):

        current_figure = plt.gcf()
        mpld3.fig_to_html()

        if self.usesMPLD3():
            plugins.connect(current_figure, plugins.MousePosition(fontsize=14))


        if self.outputType == PlotSaveTYPE.HTML_STRING:
            outString = mpld3.fig_to_html(current_figure, template_type='notebook', d3_url=self.d3js, mpld3_url=self.mpld3js)
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
                mpld3.save_html(current_figure, exactFilename, template_type='notebook', d3_url=self.d3js, mpld3_url=self.mpld3js)
            elif self.outputType == PlotSaveTYPE.JSON:
                mpld3.save_json(current_figure, exactFilename, d3_url=self.d3js, mpld3_url=self.mpld3js)
            elif self.outputType == PlotSaveTYPE.PNG:
                plt.savefig(exactFilename, transparent=self.transparent_bg, bbox_inches='tight')

            self.createdPlots.append(exactFilename)

        else:

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
