import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import genfromtxt
import os

from collections import Counter
import datetime as dt
from matplotlib.ticker import Formatter


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

        date = dt.datetime.fromtimestamp(x)

        return date.strftime(self.fmt)

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
    def plotLoadOut(cls, pore2length, title='Title', xlabel='channels', ylabel='flowcell inlet'):
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

        fig, axarr = plt.subplots(2, figsize=(10, 10), gridspec_kw = {'height_ratios':[10, 1]})
        fig.set_dpi(100)

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

        #ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
        cb1 = mpl.colorbar.ColorbarBase(axarr[1], cmap=cls.getColorMap(scolormap),
                                        norm=mpl.colors.Normalize(vmin=minAvgLength, vmax=maxAvgLength),
                                        orientation='horizontal')

        axarr[1].set_title('Legend (Avg Read Length)')

        fig.tight_layout()

        plt.show()



    @classmethod
    def plotTimeLine(cls, readsPerTime, labels, title, colors = None, bins = 100):

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


        fig, ax = plt.subplots()

        formatter = TimestampDateFormatter()
        ax.xaxis.set_major_formatter(formatter)

        linebc, bins, patches  = ax.hist(histInput, bins, histtype='bar', stacked=True, ls='dotted', color=colors, label=labels)
        ax.set_title( title )

        plt.legend()

        fig.autofmt_xdate()
        plt.show()

    @classmethod
    def plotHistogram(cls, someData, labels, title, bins = 100, xlabel=None, ylabel=None):

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

        plt.tight_layout()

        plt.show()

    @classmethod
    def plotViolin(cls, someData, labels, title, bins = 100, xlabel=None, ylabel=None):

        fig, ax = plt.subplots()
        ax.violinplot(someData, showmeans=True, showextrema=True, showmedians=True)

        ax.axes.get_xaxis().set_ticks( [i for i in range(1, len(labels)+1)] )
        ax.axes.get_xaxis().set_ticklabels( labels, rotation=90 )

        ax.axes.get_yaxis().set_ticks([i for i in range(minQual, maxQual+1)])
        ax.axes.get_yaxis().set_ticklabels([str(chr(i)) for i in range(minQual, maxQual+1)])

        plt.tight_layout()

        plt.show()





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

        plt.tight_layout()

        plt.show()

    @classmethod
    def yieldPlot(cls, dataDict, title, xlabel, ylabel):

        fig, ax = plt.subplots()

        formatter = TimestampTimeFormatter()
        ax.xaxis.set_major_formatter(formatter)


        colors = PorePlot.getColorVector(len(timeAndLength))

        for i in range(0, len(timeAndLength)):
            ax.plot(timeAndLength[i], color=colors[i])

        plt.legend( labels, loc='upper left')

        fig.autofmt_xdate()
        plt.tight_layout()
        plt.show()