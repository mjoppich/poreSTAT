from enum import Enum
import matplotlib.pyplot as plt


class PlotStyle(Enum):



class PlotStyle(Enum):
    DEFAULT=0,
    XKCD=1,
    FIVETHIRTYEIGHT=2,
    CLASSIC=3

    def __call__(self, string):
        if not name in VALID_BAR_NAMES:
            raise argparse.ArgumentError(None,'Bar can not be {n}, '
                               'it must be one of {m}'.format(
                                    n=name, m=', '.join(VALID_BAR_NAMES)))
        self.name = string

    def __init__(self):
        pass

    def __repr__(self):
        return name

class PlotConfig:

    @classmethod
    def fromParserArgs(cls, simArgs):

        pltCfg = PlotConfig()

        if simArgs.save_plot != None:
            pltCfg.saveToFile( simArgs.save_plot )

        if simArgs.plot_style != None:
            pltCfg.setStyle( simArgs.plot_style )



    @classmethod
    def addParserArgs(cls, parser):

        def toPlotStyle(input):
            return PlotStyle[input]

        parser.add_argument('-sp', '--save-plot', type=str, help='path to file where plots are saved')
        parser.add_argument('-ps', '--plot-style', type=toPlotStyle, choices=[x for x in PlotStyle], help='read types ('+ ",".join([str(x) for x in PlotStyle]) +')')

        return parser

    def __init__(self):

        self.save_to_file = False
        self.save_file = None
        self.saved_plot = 0

        self.style = None
        self.transparent_bg = True

    def saveToFile(self, filePath):

        self.save_to_file = True
        self.save_file = filePath

    def setStyle(self, style):

        self.style = style

    def makePlot(self):

        if self.style == PlotStyle.XKCD:
            plt.xkcd()
        elif self.style == PlotStyle.FIVETHIRTYEIGHT:
            plt.style.use('fivethirtyeight')
        elif self.style== PlotStyle.CLASSIC:
            plt.style.use('classic')

        plt.tight_layout()

        if self.save_to_file:

            exactFilename = self.save_file % self.saved_plot

            plt.savefig( exactFilename, transparent=self.transparent_bg)
        else:
            plt.show()
