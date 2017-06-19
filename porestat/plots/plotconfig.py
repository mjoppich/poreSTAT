import argparse
from enum import Enum
import matplotlib.pyplot as plt

class PlotStyle(Enum):
    DEFAULT='default',
    XKCD='xkcd',
    FIVETHIRTYEIGHT='fivethirtyeight',
    CLASSIC='classic',
    BMH='bmh',
    DARK='dark_background'


class PlotStyleAction(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):

        try:
            eVal = PlotStyle[values]

            self.style = eVal
            args.__dict__[ self.dest ] = self.style

        except:

            raise argparse.ArgumentError(None, 'PlotStyle can not be {n}, '
                                               'it must be one of {m}'.format(n=values,
                                                                              m=', '.join([str(x.name) for x in PlotStyle])))

    def __repr__(self):
        return self.style

class PlotConfig:

    @classmethod
    def fromParserArgs(cls, simArgs):

        pltCfg = PlotConfig()

        if simArgs.save_plot != None:
            pltCfg.saveToFile( simArgs.save_plot )

        if simArgs.plot_style != None:
            pltCfg.setStyle( simArgs.plot_style )

        return pltCfg



    @classmethod
    def addParserArgs(cls, parser):

        def toPlotStyle(input):
            return PlotStyle[input]

        parser.add_argument('-sp', '--save-plot', type=str, help='path-prefix to file where plots are saved. Final file will be save-plot.xxx.png')
        parser.add_argument('-ps', '--plot-style', action=PlotStyleAction)

        return parser

    def __init__(self):

        self.save_to_file = False
        self.save_file = None
        self.saved_plot = 0

        self.style = PlotStyle.DEFAULT
        self.transparent_bg = True

        self.createdPlots = []

    def __str__(self):

        allElements = [self.save_to_file, self.save_file, self.saved_plot, self.style, self.transparent_bg, self.createdPlots]
        allStrElements = [str(x) for x in allElements]

        return "[" + ", ".join(allStrElements) + "]"

    def saveToFile(self, filePath):

        self.save_to_file = True
        self.save_file = filePath + ".%02d.png"

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

    def makePlot(self):


        if self.save_to_file:

            exactFilename = self.save_file % self.saved_plot

            self.saved_plot += 1

            plt.savefig( exactFilename, transparent=self.transparent_bg, bbox_inches='tight')

            self.createdPlots.append(exactFilename)
        else:

            current_figure = plt.gcf()
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
