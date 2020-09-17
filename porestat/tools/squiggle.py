from porestat.plots.plotconfig import PlotConfig

from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory, PSToolException, PSToolInterface

from ..hdf5tool.Fast5File import Fast5Directory, Fast5TYPE, Fast5TYPEAction, Fast5File
import argparse
import sys, os
from ..utils.Files import eprint
import matplotlib.pyplot as plt


class SquigglePlotFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):
        super(SquigglePlotFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which + ' help')
        parser.add_argument('-r', '--read', type=argparse.FileType('r'), help='minion read folder', required=False)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):
        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return SquigglePlotter(simArgs)


class Environment(object):
    pass


class SquigglePlotter(PSToolInterface):

    def __init__(self, args):
        super(SquigglePlotter, self).__init__(args)

    def exec(self):

        print(self.args.read.name)

        read = Fast5File(self.args.read.name)
        read.type = Fast5TYPE.BASECALL_1D

        evs = read._read_raw_signal()

        for x in evs:

            plt.figure()
            plt.plot([x for x in range(0, len(x))], x)
            plt.show()



        evs = read._read_events(read.type)




        data = []

        colnames = evs.dtype.names
        for row in evs:

            rowData = {}
            for idx, col in enumerate(colnames):
                rowData[col] = row[idx]

            data.append(rowData)

        print(len(data))

        print(data[-2])
        print(data[-1])

        for i in range(0, len(data)-1):

            if not data[i]['start'] + data[i]['length'] <= data[i+1]['start']:
                continue
                print(data[i]['start'] + data[i]['length'], data[i+1]['start'])
                print("ERROR", data[i])
                print("ERROR", data[i+1])
                print()


        xvals = []
        yvals = []
        for dp in data:

            start = dp.get("start", None)
            meanVal = dp.get('mean', None)

            if start == None or meanVal == None:
                continue

            xvals.append(start)
            yvals.append(meanVal)


        minx = min(xvals)
        xvals = [x-minx for x in xvals]

        duration = max(xvals)-min(xvals)
        print(duration)

        def slidingWindow(sequence, winSize, step=1):
            """Returns a generator that will iterate through
            the defined chunks of input sequence.  Input sequence
            must be iterable."""

            # Verify the inputs
            try:
                it = iter(sequence)
            except TypeError:
                raise Exception("**ERROR** sequence must be iterable.")

            winSize -= 1

            lbefore = winSize // 2
            lafter = winSize // 2

            if winSize % 2 == 1:
                lbefore += 1

            # Do the work
            for i in range(0, len(sequence)):

                istart = max(0, i-lbefore)
                iend = min(i+lafter, len(sequence))

                yield sequence[istart:iend]


        winXVals = []
        winYVals = []

        lastEnd = data[idx]

        steps = 10
        stepDur = duration / steps
        winFacets = []
        limits = []
        for i in range(1,steps+1):
            limits.append(i * stepDur)

        curFacet = 0

        for idx,x in enumerate(slidingWindow(yvals, 5)):
            level=sum(x) / len(x)

            xstart = data[idx]['start']

            if idx+1 < len(data):
                xend = data[idx+1]['start']
            else:
                xend = xstart + data[idx]['length']

            xstart -= minx
            xend -= minx

            winFacets.append(curFacet)
            winFacets.append(curFacet)

            winXVals.append(xstart)
            winXVals.append(xend)

            winYVals.append(level)
            winYVals.append(level)

            if xend > limits[curFacet]:
                curFacet+=1


        import seaborn as sns
        import pandas as pd

        df = pd.DataFrame.from_dict({
            'time': winXVals,
            'Current': winYVals,
            'facet': winFacets
        })

        self.args.pltcfg.startPlot()



        g = sns.FacetGrid(df, row="facet", size=8, sharex=False)
        fig = plt.gcf()
        fig.set_tight_layout(True)

        g = g.map(plt.step, "time", "Current")

        for i, nax in enumerate([x for x in g.axes]):
            ax = nax[0]
            ax.set_title("Squiggle Plot Read {} (Duration {:.5}s, {:.5}s-{:.5}s)".format(os.path.basename(self.args.read.name), duration, i*stepDur, (i+1)*stepDur))
            ax.set_xlim(i*stepDur, (i+1) * stepDur)

        self.args.pltcfg.makePlot(noTightLayout=True)

