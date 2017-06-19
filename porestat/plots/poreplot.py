import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import genfromtxt
import os

from collections import Counter
import datetime as dt
from matplotlib.ticker import Formatter
from .plotconfig import PlotConfig
from matplotlib import gridspec


class TimestampDateFormatter(Formatter):
    def __init__(self, fmt='%Y-%m-%d %H:%M:%S'):
        self.fmt = fmt

    def __call__(self, x, pos=0):
        'Return the label for time x at position pos'

        date = dt.datetime.fromtimestamp(x)

        return date.strftime(self.fmt)

class TimestampTimeFormatter(Formatter):
    def __init__(self, fmt='%H:%M:%S'):
        self.fmt = fmt

    def __call__(self, x, pos=0):
        'Return the label for time x at position pos'

        m, s = divmod(x, 60)
        h, m = divmod(m, 60)

        return "%d:%02d:%02d" % (h, m, s)

class PlotDirectionTYPE(Enum):
    VERTICAL='vertical'
    HORIZONTAL='horizontal'

class PorePlot:

    @classmethod
    def getColorVector(cls, elems, colormap="Viridis"):

        colors = []
        for i in range(0, elems):
            color = cls.getColorLin(0, elems, i)
            colors.append(color)
        return colors

    @classmethod
    def getColorMap(cls, colormap="Viridis"):
        cmap = plt.cm.get_cmap(colormap)
        return cmap

    @classmethod
    def getColor(cls, colormap="Viridis", value=0.5):
        cmap = cls.getColorMap(colormap=colormap)

        return cmap(value)

    @classmethod
    def getColorLin(cls, min, max, val, colormap="viridis"):

        value = val / (max-min)

        return cls.getColor(colormap=colormap, value=value)

    @classmethod
    def plotLoadOut(cls, pore2length, title='Title', xlabel='channels', ylabel='flowcell inlet', pltcfg = PlotConfig()):
        """

        :param pore2length: dictionary channelID -> read-length
        :param title: title of the plot
        :param xlabel: xlabel of the plot
        :param ylabel: ylabel of the plot
        :return:
        """

        p2c = Counter()
        p2l = {}

        minCount = 0
        maxCount = 0

        this_dir, this_filename = os.path.split(__file__)
        chipLayout = genfromtxt(this_dir + '/../data/chip_layout.csv', delimiter=',', dtype=int)
        pores = chipLayout.shape
        minPore = np.amin(chipLayout)
        maxPore = np.amax(chipLayout)

        minAvgLength = 1000000
        maxAvgLength = 0

        scolormap = "plasma"

        for i in range(1, maxPore+1):

            p2c[i] = 0

            if i in pore2length:
                count = len(pore2length[i])
                p2c[i] = count

                if count > 0:
                    average =np.average(pore2length[i])
                    p2l[i] = average
                    minAvgLength = min(minAvgLength, average)
                    maxAvgLength = max(maxAvgLength, average)
                else:
                    p2l[i] = -1

                minCount = min(count, minCount)
                maxCount = max(count, maxCount)

            else:
                p2c[i] = 0
                p2l[i] = -1

        minPoreRad = 50
        maxPoreRad = 200

        p2area = {}

        for x in p2c:

            if p2c[x] == 0:
                p2area[x] = minPoreRad / 2
                continue

            count = p2c[x] - minCount
            frac = count / (maxCount - minCount)
            area = minPoreRad + frac * (maxPoreRad-minPoreRad)
            p2area[x] = area

        p2coord = {}

        offset = -1
        if minPore == 0:
            offset = 0

        for x in p2c:
            row, col = divmod(x+offset, pores[0])
            p2coord[x] = (row, col)


        xvec = []
        yvec = []
        area = []
        color = []

        for i in range(0, chipLayout.shape[0]):
            for j in range(0,chipLayout.shape[1]):

                channelID = int(chipLayout[i,j])

                xvec.append(i)
                yvec.append(j)
                area.append(p2area[channelID])

                if p2l[channelID] == -1:
                    color.append((1, 1, 1, 1))

                else:
                    frac = p2l[channelID] / (maxAvgLength - minAvgLength)
                    color.append(cls.getColor(colormap=scolormap, value=frac))

        pltcfg.startPlot()

        fig = plt.figure(figsize=(10,10))
        fig.set_dpi(100)

        gs = gridspec.GridSpec(2,1, height_ratios=[10, 1])

        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])

        axarr = [ ax0, ax1 ]

        axarr[0].set_title( title )
        axarr[0].set_xlabel( xlabel , fontsize=18)
        axarr[0].set_ylabel( ylabel , fontsize=18)

        axarr[0].axes.get_xaxis().set_ticks( [x for x in range(0, pores[0])] )
        axarr[0].axes.get_yaxis().set_ticks( [x for x in range(0, pores[1])])

        axarr[0].axes.get_xaxis().set_ticklabels( [str(x) for x in range(1, pores[0]+1)] )
        axarr[0].axes.get_yaxis().set_ticklabels( [str(x) for x in range(1, pores[1]+1)])

        axarr[0].tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='on', labeltop='off',
                        labelright='off', labelbottom='on')

        #axarr[0].Axes(fig, [1, pores[0], 1, pores[1]])
        axarr[0].scatter(xvec, yvec, s=area, c=color, alpha=0.5)

        sizeSteps = 5
        sizes = [ x for x in range(int(min(area)), int(max(area)), int((max(area)-min(area))/sizeSteps)) ]

        legendPlts = []
        for size in sizes:
            obj =axarr[0].scatter([], [], s=size, marker='o', color='#555555')
            legendPlts.append(obj)

        sizeLabels = [str(x) for x in sizes]
        sizeLabels[0] += " reads"
        plt.legend(tuple(legendPlts), tuple(sizeLabels), scatterpoints=1, ncol=sizeSteps, loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True)

        #ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
        cb1 = mpl.colorbar.ColorbarBase(axarr[1], cmap=cls.getColorMap(scolormap),
                                        norm=mpl.colors.Normalize(vmin=minAvgLength, vmax=maxAvgLength),
                                        orientation='horizontal')

        axarr[1].set_title('Legend (Avg Read Length)')

        pltcfg.makePlot()



    @classmethod
    def plotTimeLine(cls, readsPerTime, labels, title, colors = None, bins = 100, pltcfg = PlotConfig()):

        """

        :param readsPerTime: vector of counters of timestamp -> #reads
        :param labels: vector of string with label for each element in readsPerTime
        :return:
        """

        histInput = []

        for x in readsPerTime:

            readsTime = x

            timePoints = sorted(list(readsTime.keys()))

            histInTime = []
            for x in readsTime:
                for i in range(0, readsTime[x]):
                    histInTime.append(x)

            histInput.append( histInTime )

        if colors == None or len(colors) != len(histInput):

            colors = []
            for i in range(0, len(histInput)):

                color = cls.getColorLin(0, len(histInput)-1, i)
                colors.append(color)

        pltcfg.startPlot()
        fig, ax = plt.subplots()

        formatter = TimestampDateFormatter()
        ax.xaxis.set_major_formatter(formatter)

        linebc, bins, patches  = ax.hist(histInput, bins, histtype='bar', stacked=True, ls='dotted', color=colors, label=labels)
        ax.set_title( title )

        plt.legend()

        fig.autofmt_xdate()
        pltcfg.makePlot()

    @classmethod
    def plotHistogram(cls, someData, labels, title, bins = 100, xlabel=None, ylabel=None, pltcfg = PlotConfig()):
        pltcfg.startPlot()
        fig, ax = plt.subplots()
        linebc, bins, patches = ax.hist( someData , bins, histtype='bar', stacked=False, label=labels)
        ax.set_title( title )

        if xlabel != None:
            ax.set_xlabel( xlabel )

        if ylabel != None:
            ax.set_ylabel( ylabel )

#        ax.axes.get_xaxis().set_ticks( [i for i in range(1, len(stepLabels)+1)] )
#        ax.axes.get_xaxis().set_ticklabels( stepLabels, rotation=90 )

        if len(someData)>1:
            ax.legend()

        pltcfg.makePlot()

    @classmethod
    def plotCumHistogram(cls, someData, labels, title, bins=100, xlabel=None, ylabel=None, pltcfg=PlotConfig()):

        pltcfg.startPlot()
        fig, ax = plt.subplots()

        dataVec = []
        for x in labels:
            dataVec.append( someData[x] )

        linebc, bins, patches = ax.hist(dataVec, bins, histtype='step', cumulative=1, normed=1, stacked=False, label=labels)
        ax.set_title(title)

        if xlabel != None:
            ax.set_xlabel(xlabel)

        if ylabel != None:
            ax.set_ylabel(ylabel)

            #        ax.axes.get_xaxis().set_ticks( [i for i in range(1, len(stepLabels)+1)] )
            #        ax.axes.get_xaxis().set_ticklabels( stepLabels, rotation=90 )

        if len(someData) > 1:
            ax.legend()

        plt.xscale('log')

        pltcfg.makePlot()

    @classmethod
    def plotSingleViolin(cls, data, title, ax):

        ax.violinplot(data, showmeans=True, showextrema=True, showmedians=True, points=500)
        ax.set_title(title)

    @classmethod
    def plotViolin(cls, someData, labels, title, pltcfg = PlotConfig(), axisManipulation = None, plotDirection=PlotDirectionTYPE.VERTICAL):

        if type(someData) == list:
            shape = (1, 1)
            someData = { title: someData }
        else:
            elems = len(someData)
            
            if plotDirection == PlotDirectionTYPE.VERTICAL:
                shape = (1, elems)
            else:
                shape = (elems, 1)

        pltcfg.startPlot()
        fig, ax = plt.subplots(nrows=shape[0], ncols=shape[1])

        if labels == None:
            labels = [x for x in someData]

        for i in range(0, len(labels)):

            x = labels[i]
            pos = divmod(i, shape[1])

            axisManipulator = axisManipulation[labels[i]] if not axisManipulation is None and labels[i] in axisManipulation else None

            cls.plotSingleViolin( someData[x], labels[i], ax[i] )
            
            if axisManipulator != None:
                axisManipulator(ax[i])

        plt.tight_layout()
        pltcfg.makePlot()

    @classmethod
    def plotSingleBoxplot(cls, data, title, ax):
        ax.boxplot(data, notch=True, patch_artist=True)
        ax.set_title(title)

    @classmethod
    def plotBoxplot(cls, someData, labels, title, pltcfg = PlotConfig(),  axisManipulation = None, plotDirection=PlotDirectionTYPE.VERTICAL):

        if type(someData) == list:
            shape = (1, 1)
            someData = { title: someData }
        else:
            elems = len(someData)

            if plotDirection == PlotDirectionTYPE.VERTICAL:
                shape = (1, elems)
            else:
                shape = (elems, 1)

        pltcfg.startPlot()
        fig, ax = plt.subplots(nrows=shape[0], ncols=shape[1])

        if labels == None:
            labels = [x for x in someData]

        for i in range(0, len(labels)):

            x = labels[i]
            pos = divmod(i, shape[1])

            axisManipulator = axisManipulation[labels[i]] if not axisManipulation is None and labels[i] in axisManipulation else None

            cls.plotSingleBoxplot( someData[x], labels[i], ax[i] )
            
            if axisManipulator != None:
                axisManipulator(ax[i])

        plt.tight_layout()
        pltcfg.makePlot()

    @classmethod
    def plotSingleBarplot(cls, barx, bary, title, ax):
        ax.boxplot( barx, bary )
        ax.set_title(title)

    @classmethod
    def plotBarplot(cls, someData, labels, title, pltcfg = PlotConfig(),  axisManipulation = None, plotDirection=PlotDirectionTYPE.VERTICAL):

        if type(someData) == list:
            shape = (1, 1)
            someData = { title: someData }
        else:
            elems = len(someData)

            if plotDirection == PlotDirectionTYPE.VERTICAL:
                shape = (1, elems)
            else:
                shape = (elems, 1)

        pltcfg.startPlot()
        fig, ax = plt.subplots(nrows=shape[0], ncols=shape[1])

        if labels == None:
            labels = [x for x in someData]

        for i in range(0, len(labels)):

            x = labels[i]
            pos = divmod(i, shape[1])

            axisManipulator = axisManipulation[labels[i]] if not axisManipulation is None and labels[i] in axisManipulation else None

            xvals = [key for key in someData[x]]
            yvals = [someData[x][key] for key in someData[x]]

            cls.plotSingleBarplot( xvals, yvals, labels[i], ax[i] )
            
            if axisManipulator != None:
                axisManipulator(ax[i])

        plt.tight_layout()
        pltcfg.makePlot()

    @classmethod
    def yieldPlot(cls, dataDict, title, xlabel, ylabel, pltcfg = PlotConfig()):

        pltcfg.startPlot()
        fig, ax = plt.subplots()

        formatter = TimestampTimeFormatter()
        ax.xaxis.set_major_formatter(formatter)

        labels = sorted([x for x in dataDict])
        vDataSeries = []

        for series in labels:
            vDataSeries.append( sorted(dataDict[series], key=lambda tup: tup[0]) )

        colors = PorePlot.getColorVector(len(labels))

        for i in range(0, len(labels)):
            xes = [x[0] for x in vDataSeries[i]]
            yes = [x[1] for x in vDataSeries[i]]

            ax.plot(xes, yes, color=colors[i])

        ax.set_xlabel( xlabel )
        ax.set_ylabel( ylabel )
        ax.set_title(title)

        # Put a legend below current axis
        ax.legend(labels, loc='upper center', bbox_to_anchor=(0.5, -0.2),
                  fancybox=True, shadow=True, ncol=1)

        fig.autofmt_xdate()

        pltcfg.makePlot()