from enum import Enum

import mpld3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import re

import pandas
from mpld3 import plugins
from mpld3.plugins import PointHTMLTooltip
from numpy import genfromtxt
import os
import math

from collections import Counter, OrderedDict
import datetime as dt
from matplotlib.ticker import Formatter

from porestat.utils.OrderedSet import OrderedSet
from .plotconfig import PlotConfig
from matplotlib import gridspec, colors
import scipy.cluster.hierarchy as sch
import seaborn as sns


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

        return self.format(x)

    def format(self, x):
        m, s = divmod(x, 60)
        h, m = divmod(m, 60)

        return "%d:%02d:%02d" % (h, m, s)

class MultiAxesPointHTMLTooltip(mpld3.plugins.PluginBase):
    """A Plugin to enable an HTML tooltip:
    formated text which hovers over points.

    Parameters
    ----------
    points : matplotlib Collection or Line2D object
        The figure element to apply the tooltip to
    labels : list
        The labels for each point in points, as strings of unescaped HTML.
    hoffset, voffset : integer, optional
        The number of pixels to offset the tooltip text.  Default is
        hoffset = 0, voffset = 10
    css : str, optional
        css to be included, for styling the label html if desired
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from mpld3 import fig_to_html, plugins
    >>> fig, ax = plt.subplots()
    >>> points = ax.plot(range(10), 'o')
    >>> labels = ['<h1>{title}</h1>'.format(title=i) for i in range(10)]
    >>> plugins.connect(fig, MultiAxesPointHTMLTooltip(points[0], labels))
    >>> fig_to_html(fig)
    """

    JAVASCRIPT = """
    mpld3.register_plugin("htmltooltip", HtmlTooltipPlugin);
    HtmlTooltipPlugin.prototype = Object.create(mpld3.Plugin.prototype);
    HtmlTooltipPlugin.prototype.constructor = HtmlTooltipPlugin;
    HtmlTooltipPlugin.prototype.requiredProps = ["ids", "labels"];
    HtmlTooltipPlugin.prototype.defaultProps = {hoffset:0,
                                                voffset:10};
    function HtmlTooltipPlugin(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };

    HtmlTooltipPlugin.prototype.draw = function(){

           var tooltip = d3.select("body").append("div")
                    .attr("class", "mpld3-tooltip")
                    .style("position", "absolute")
                    .style("z-index", "10")
                    .style("visibility", "hidden");

       combinedElements = null;
       for (var i = 0; i < this.props.ids.length; i++)
       {

        console.log("adding info for element: " + this.props.ids[i] );

       var obj = mpld3.get_element(this.props.ids[i]);
       var labels = this.props.labels[i];
       var self = this;
       let pltelem = i;

       obj.elements().on("mouseover", function(d, i){
       console.log(pltelem);
                              tooltip.html(self.props.labels[pltelem][i])
                                     .style("visibility", "visible");})
           .on("mousemove", function(d, i){
                  tooltip
                    .style("top", d3.event.pageY + this.props.voffset + "px")
                    .style("left",d3.event.pageX + this.props.hoffset + "px");
                 }.bind(this))
           .on("mouseout",  function(d, i){
                           tooltip.style("visibility", "hidden");});

       }


    };
    """

    def __init__(self, points, labels=None,
                 hoffset=0, voffset=10, css=None):
        self.points = points
        self.labels = labels
        self.voffset = voffset
        self.hoffset = hoffset
        self.css_ = css or ""

        from matplotlib.lines import Line2D

        ids = []
        for point in points:

            if isinstance(point, Line2D):
                suffix = "pts"
            else:
                suffix = None

            points_id = self.get_id(point, suffix)
            ids.append(points_id)

        print(ids)

        self.dict_ = {"type": "htmltooltip",
                      # "id": get_id(points, suffix),
                      "ids": ids,
                      "labels": labels,
                      "hoffset": hoffset,
                      "voffset": voffset}

    def get_id(self, obj, suffix="", prefix="el", warn_on_invalid=True):
        """Get a unique id for the object"""
        if not suffix:
            suffix = ""
        if not prefix:
            prefix = ""

        objid = prefix + str(os.getpid()) + str(id(obj)) + suffix

        if warn_on_invalid and not self.html_id_ok(objid):
            print('"{0}" is not a valid html ID. This may cause problems')

        return objid

    def html_id_ok(self, objid, html5=False):
        """Check whether objid is valid as an HTML id attribute.

        If html5 == True, then use the more liberal html5 rules.
        """
        if html5:
            return not re.search('\s', objid)
        else:
            return bool(re.match("^[a-zA-Z][a-zA-Z0-9\-\.\:\_]*$", objid))



class PlotDirectionTYPE(Enum):
    VERTICAL='vertical'
    HORIZONTAL='horizontal'

class PorePlot:

    @classmethod
    def getToolTipCSS(cls):

        return """
        table.tooltip
        {
          border-collapse: collapse;
        }
        .tooltip th
        {
          color: #ffffff;
          background-color: #000000;
        }
        .tooltip td
        {
          background-color: #cccccc;
        }
        table.tooltip, .tooltip th, .tooltip td
        {
          font-family:Arial, Helvetica, sans-serif;
          border: 1px solid black;
          text-align: right;
        }
        """

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



        def makeHTMLTable( poreNum, poreReads, poreAvgLength):

            poreNumLine = "<tr><th>Pore #</th><th>{pn}</th></tr>".format(pn=poreNum)
            poreReadsLine = "<tr><td>Reads</td><td>{rc}</td></tr>".format(rc=poreReads)
            poreLengthLine = "<tr><td>Avg Length</td><td>{:.2f}</td></tr>".format(poreAvgLength)

            htmlStr = "<table class='tooltip'>"+poreNumLine+poreReadsLine+poreLengthLine+"</table>"

            return htmlStr

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


        htmlDescr = []
        for i in range(0, chipLayout.shape[0]):
            for j in range(0,chipLayout.shape[1]):

                channelID = int(chipLayout[i,j])
                htmlDescr.append( makeHTMLTable(channelID, p2c[channelID], p2l[channelID]) )


        #axarr[0].Axes(fig, [1, pores[0], 1, pores[1]])
        elems = axarr[0].scatter(xvec, yvec, s=area, c=color, alpha=0.5)


        if pltcfg.usesMPLD3():
            tooltip = MultiAxesPointHTMLTooltip([elems], [htmlDescr], voffset=10, hoffset=10, css=cls.getToolTipCSS())
            plugins.connect(fig, tooltip)
        else:

            sizeSteps = 5

            frac = count / (maxCount - minCount)
            area = minPoreRad + frac * (maxPoreRad-minPoreRad)

            readCounts = [ x for x in range(int(minCount), int(maxCount), (int(maxCount)-int(minCount))/ sizeSteps)]
            sizes = [x for x in range(int(min(area)), int(max(area)), int((max(area) - min(area)) / sizeSteps))]

            legendPlts = []
            for size in sizes:
                obj = axarr[0].scatter([], [], s=size, marker='o', color='#555555')
                legendPlts.append(obj)

            sizeLabels = [str(x) for x in readCounts]
            sizeLabels[0] += " reads"

            cls.makeLegend(fig, axarr[0], tuple(legendPlts), tuple(sizeLabels), pltcfg, bbox_to_anchor=(0.5, -0.2))


        #ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
        cb1 = mpl.colorbar.ColorbarBase(axarr[1], cmap=cls.getColorMap(scolormap),
                                        norm=mpl.colors.Normalize(vmin=minAvgLength, vmax=maxAvgLength),
                                        orientation='horizontal')

        axarr[1].set_title('Legend (Avg Read Length)')

        pltcfg.makePlot()


    @classmethod
    def plot_scatter_densities(cls, xdata, ydata, title, xlabel, ylabel, addInfos=None, iRoundTo = 3, figSize=(30,30), pltcfg=PlotConfig(), textIsNumber=True):

        fig, ax = plt.subplots()

        x = np.asarray(xdata)
        y = np.asarray(ydata)

        print("Going to display elements", len(xdata), len(ydata), min(xdata), max(xdata), min(ydata), max(ydata))

        # Set bin numbers in both axes
        nxbins = 50
        nybins = 50

        # Set the cutoff for resolving the individual points
        minperbin = 1

        # Make the density histrogram
        H, xedges, yedges = np.histogram2d(x,y, bins=(nxbins, nybins))
        # Reorient the axes
        #H = H[::-1]

        def transformX(xnew):
            val = ((xnew - extent[0]) / (extent[1] - extent[0])) * nxbins
            return val

        def transformY(ynew):
            val = ((ynew - extent[2]) / (extent[3] - extent[2])) * nybins
            return val

        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        print("extent", extent)

        # Figure out which bin each x,y point is in
        xbinsize = xedges[1] - xedges[0]
        ybinsize = yedges[1] - yedges[0]


        print("binsize", xbinsize, ybinsize)


        xi = ((x - xedges[0]) / xbinsize).astype(np.integer)
        yi = nybins - 1 - ((y - yedges[0]) / ybinsize).astype(np.integer)

        # Subtract one from any points exactly on the right and upper edges of the region
        xim1 = xi - 1
        yim1 = yi - 1
        xi = np.where(xi < nxbins, xi, xim1)
        yi = np.where(yi < nybins, yi, yim1)

        myCounts = np.zeros((nxbins, nybins))

        for i in range(0, len(x)):

            binx = xi[i]
            biny = yi[i]

            myCounts[binx, biny] += 1

        lowdensityx = []
        lowdensityy = []
        htmlDescr = []


        def makeHTML(xlabel, ylabel, xvalue, yvalue):

            allLines = []

            for elem, val in [(xlabel, xvalue), (ylabel, yvalue)]:

                if isinstance(val, (float, int)):
                    allLines.append("<tr><td>{}</td><td>{:.6f}</td></tr>".format(elem, val))
                else:
                    allLines.append("<tr><td>{}</td><td>{}</td></tr>".format(elem, val))

            htmlStr = "<table class='tooltip'>"+ " ".join(allLines) +"</table>"

            return htmlStr


        for i in range(0, len(x)):

            elemx = x[i]
            elemy = y[i]

            binx = xi[i]
            biny = yi[i]

            if myCounts[binx, biny] <= 2:
                lowdensityx.append(elemx)
                lowdensityy.append(elemy)

                htmlDescr.append(makeHTML(xlabel, ylabel, elemx, elemy))

        for i in range(0, nxbins):
            for j in range(0, nybins):
                if myCounts[i,j] <= 2:
                    myCounts[i,j] = 0


        cmap = PorePlot.getColorMap("Blues")

        lowdensityx = [transformX(xe) for xe in lowdensityx]
        lowdensityy = [transformY(ye) for ye in lowdensityy]


        myCounts = myCounts.transpose()
        myCounts = myCounts[::-1]

        cp1 = ax.imshow(myCounts, interpolation='nearest', origin="lower", norm=colors.LogNorm(), cmap=cmap)
        fig.colorbar(cp1)


        ax.set_xlim([-1, nxbins])
        ax.set_ylim([-1, nybins])


        if len(lowdensityx) < 10000:

            elems = ax.plot(lowdensityx, lowdensityy, linestyle='None', marker='o', mfc='k', mec='k', ms=3)

            if pltcfg.usesMPLD3() and len(htmlDescr) > 0:
                tooltip = MultiAxesPointHTMLTooltip([elems[0]], [htmlDescr], voffset=10, hoffset=10, css=cls.getToolTipCSS())
                plugins.connect(fig, tooltip)

        ax.set_title(title, size=20)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        tickCount = 10
        xTickSize = nxbins/tickCount
        yTickSize = nybins/tickCount

        xticks = [x * xTickSize for x in range(0, tickCount)]
        yticks = [x * yTickSize for x in range(0, tickCount)]

        ax.set_xticks(xticks)
        ax.set_yticks(yticks)

        xvel = extent[1]-extent[0]
        yvel = extent[3]-extent[2]

        axfmt = "{:6.3f}"

        ax.set_xticklabels( [axfmt.format(  extent[0] + (x/nxbins) *xvel ) for x in xticks] )
        ax.set_yticklabels( [axfmt.format(  extent[2] + (x/nxbins) *yvel ) for x in yticks] )

        plugins.connect(fig, plugins.MousePosition(fontsize=14))

        pltcfg.makePlot(noTightLayout=True)

    @classmethod
    def plotscatter(cls, xdata, ydata, title, xlabel, ylabel, addInfos=None, pltcfg = PlotConfig()):

        def makeHTML(addinfo):

            allLines = []

            for elem in addinfo:

                val = addinfo[elem]

                if isinstance(val, (float, int)):
                    allLines.append("<tr><td>{}</td><td>{:.6f}</td></tr>".format(elem, val))
                else:
                    allLines.append("<tr><td>{}</td><td>{}</td></tr>".format(elem, val))

            htmlStr = "<table class='tooltip'>"+ " ".join(allLines) +"</table>"

            return htmlStr

        htmlDescr = []

        for i in range(0, len(xdata)):
            if addInfos != None:
                htmlDescr.append( makeHTML( addInfos[i] ))

        pltcfg.startPlot()
        fig, ax = plt.subplots()

        elems = ax.plot(xdata, ydata, c='b', alpha=0.5, marker='o', linestyle='')

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        if pltcfg.usesMPLD3() and len(htmlDescr) > 0:
            tooltip = MultiAxesPointHTMLTooltip([elems[0]], [htmlDescr], voffset=10, hoffset=10, css=cls.getToolTipCSS())
            plugins.connect(fig, tooltip)

        pltcfg.makePlot()


    @classmethod
    def volcanoPlot(cls, genenames, foldchanges, pvalues, title, xlabel, ylabel, pltcfg = PlotConfig()):

        if len(genenames) != len(foldchanges) or len(genenames) != len(pvalues):
            raise ValueError("genenames, foldchanges and pvalues must have same length, but is " + str((len(genenames), len(foldchanges), len(pvalues))))


        def logpv(pv):

            if pv == None or pv == 0:
                return 0
            return math.log2(pv)

        def makeHTML(gene, fc, pv):
            geneLine = "<tr><th>Gene</th><th>{gene}</th></tr>".format(gene=gene)
            fcLine = "<tr><td>log2 FC</td><td>{:.6f}</td></tr>".format(fc)
            logpvLine = "<tr><td>log2 PV</td><td>{:.6f}</td></tr>".format(-logpv(pv))

            pvalLine = "<tr><td>p-Value</td><td>{:.6f}</td></tr>".format(pv)

            htmlStr = "<table class='tooltip'>"+geneLine+fcLine+logpvLine+pvalLine+"</table>"

            return htmlStr



        htmlDescr = []
        fcData = []
        pvData = []
        colors = []

        signfcData = []
        signpvData = []
        signHTMLDescr = []
        signColors = []

        for i in range(0, len(genenames)):
            fc = foldchanges[i]
            pv = pvalues[i]

            if fc == None or pv == None:
                continue



            if pv < 0.05:
                signfcData.append(fc)
                signpvData.append(-logpv(pv))
                signColors.append('red')
                signHTMLDescr.append( makeHTML( genenames[i], fc, pv ))

            else:
                fcData.append(fc)
                pvData.append(-logpv(pv))
                colors.append('blue')
                htmlDescr.append( makeHTML( genenames[i], fc, pv ))

        pltcfg.startPlot()
        fig, ax = plt.subplots()

        elems = ax.plot(fcData, pvData, c='b', alpha=0.5, marker='o', linestyle='')
        elemsSign = ax.plot(signfcData, signpvData, c='r', alpha=0.5, marker='o', linestyle='')

        print(len(signHTMLDescr))
        print(len(htmlDescr))

        ax.set_xlabel(xlabel)
        ax.set_ylabel("neg " + ylabel)
        ax.set_title(title)

        if pltcfg.usesMPLD3():
            tooltip = MultiAxesPointHTMLTooltip([elems[0], elemsSign[0]], [htmlDescr, signHTMLDescr], voffset=10, hoffset=10, css=cls.getToolTipCSS())
            plugins.connect(fig, tooltip)

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
    def plotHistogram(cls, someData, labels, title, bins = 100, xlabel=None, ylabel=None, plotDirection=PlotDirectionTYPE.VERTICAL, pltcfg = PlotConfig()):

        if type(someData) == dict:

            if labels == None:
                labels = [x for x in someData]

            someData = [someData[x] for x in labels]


        plotDir = 'horizontal'
        if plotDirection ==  PlotDirectionTYPE.VERTICAL:
            plotDir = 'vertical'


        pltcfg.startPlot()
        fig, ax = plt.subplots()
        linebc, bins, patches = ax.hist( someData , bins, histtype='bar', stacked=False, label=labels, orientation=plotDir)
        ax.set_title( title )

        if xlabel != None:
            ax.set_xlabel( xlabel )

        if ylabel != None:
            ax.set_ylabel( ylabel )

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_rotation(65)

#        ax.axes.get_xaxis().set_ticks( [i for i in range(1, len(stepLabels)+1)] )
#        ax.axes.get_xaxis().set_ticklabels( stepLabels, rotation=90 )

        if len(someData)>1:
            ax.legend()

        pltcfg.makePlot()

    @classmethod
    def plotCumHistogram(cls, someData, labels, title, bins = 100, xLogAxis=True,xlabel=None, normed=True, ylabel=None, plotDirection=PlotDirectionTYPE.VERTICAL, pltcfg = PlotConfig()):

        pltcfg.startPlot()
        fig, ax = plt.subplots()

        if type(someData) == dict:

            if labels == None:
                labels = [x for x in someData]

            someData = [someData[x] for x in labels]


        plotDir = 'horizontal'
        if plotDirection ==  PlotDirectionTYPE.VERTICAL:
            plotDir = 'vertical'

        if bins == -1:
            bins = max([len(x) for x in someData])


        linebc, bins, patches = ax.hist(someData, bins, histtype='step', cumulative=1, normed=normed, stacked=False, label=labels, orientation=plotDir)
        ax.set_title(title)

        if xlabel != None:
            ax.set_xlabel(xlabel)

        if ylabel != None:
            ax.set_ylabel(ylabel)

            #        ax.axes.get_xaxis().set_ticks( [i for i in range(1, len(stepLabels)+1)] )
            #        ax.axes.get_xaxis().set_ticklabels( stepLabels, rotation=90 )

        if len(someData) > 1:
            ax.legend()

        if xLogAxis:
            plt.xscale('log')

        pltcfg.makePlot()

    @classmethod
    def plotSingleViolin(cls, data, title, ax, vert=True):

        if len(data) == 0:
            plotData = np.array([float('nan'), float('nan')], dtype=float)
            plotPos = [1]
        else:
            if type(data[0]) == list:

                plotData = data
                plotPos = [i for i in range(1, len(plotData)+1)]

                for i in range(0, len(plotData)):

                    if type(plotData[i]) == list:
                        if len(plotData[i]) == 0:
                            plotData[i] = np.array([float('nan'), float('nan')], dtype=float)
                        else:
                            plotData[i] = np.array(plotData[i], dtype=float)

            else:
                plotPos = [1]
                plotData = np.array(data, dtype=float)

        ax.violinplot(plotData, positions=plotPos, showmeans=True, showextrema=True, showmedians=True, points=100, vert=vert, bw_method=0.1)
        ax.set_title(title)

    @classmethod
    def plotViolinSNS(cls, someData, labels, title, pltcfg=PlotConfig(), xTitle=None, yTitle=None,axisManipulation=None,
                   plotDirection=PlotDirectionTYPE.VERTICAL, shareX=None, shareY=None):

        if type(someData) == list:
            someData = {title: someData}
        else:

            if plotDirection == PlotDirectionTYPE.VERTICAL:
                orientation='v'
            else:
                orientation = 'h'

        pltcfg.startPlot()


        fig, ax = plt.subplots()

        if labels == None:
            labels = [x for x in someData]

        pdData = OrderedDict()

        for i in range(0, len(labels)):

            x = labels[i]
            strX = str(x)
            pdData[strX] = someData[x]


        pandasDF = pandas.DataFrame.from_dict(pdData,orient='index').T.dropna()

        sns.violinplot(data=pandasDF, orient=orientation, ax=ax)

        if xTitle != None:
            ax.set_xlabel(xTitle)

        if yTitle != None:
            ax.set_ylabel(yTitle)

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_rotation(65)

        fig.suptitle(title)
        pltcfg.makePlot()


    @classmethod
    def plotViolin(cls, someData, labels, title, pltcfg = PlotConfig(), axisManipulation = None, plotDirection=PlotDirectionTYPE.VERTICAL, shareX=None, shareY=None):

        if type(someData) == list:
            shape = (1, 1)
            someData = { title: someData }
        else:
            elems = len(someData)
            
            if plotDirection == PlotDirectionTYPE.VERTICAL:
                shape = (1, elems)

                if shareX == None:
                    shareX = False

                if shareY == None:
                    shareY = True


                vert=True

            else:
                shape = (elems, 1)

                if shareX == None:
                    shareX = True

                if shareY == None:
                    shareY = False


                vert=False

        pltcfg.startPlot()

        fig, ax = plt.subplots(nrows=shape[0], ncols=shape[1], sharex=shareX, sharey=shareY)

        if labels == None:
            labels = [x for x in someData]

        for i in range(0, len(labels)):

            x = labels[i]
            axisManipulator = axisManipulation[labels[i]] if not axisManipulation is None and labels[i] in axisManipulation else None

            if shape == (1,1):
                cls.plotSingleViolin( someData[x], labels[i], ax, vert=vert)
            else:
                cls.plotSingleViolin( someData[x], labels[i], ax[i], vert=vert)
            
            if axisManipulator != None:
                axisManipulator(ax[i])

        fig.suptitle(title)
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

        if shape == (1,1):
            ax = [ax]

        if labels == None:
            labels = [x for x in someData]

        for i in range(0, len(labels)):

            x = labels[i]
            pos = divmod(i, shape[1])

            axisManipulator = axisManipulation[labels[i]] if not axisManipulation is None and labels[i] in axisManipulation else None

            cls.plotSingleBoxplot( someData[x], labels[i], ax[i] )
            ax[i].set_title(labels[i])


            if axisManipulator != None:
                axisManipulator(ax[i])

        fig.suptitle(title)

        plt.tight_layout()
        pltcfg.makePlot()

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

        fig.suptitle(title)

        plt.tight_layout()
        pltcfg.makePlot()

    @classmethod
    def yieldPlot(cls, dataDict, title, xlabel, ylabel, pltcfg = PlotConfig()):

        pltcfg.startPlot()
        fig, ax = plt.subplots()

        formatter = TimestampTimeFormatter()

        labels = sorted([x for x in dataDict])
        vDataSeries = []

        for series in labels:
            vDataSeries.append( sorted(dataDict[series], key=lambda tup: tup[0]) )

        colors = PorePlot.getColorVector(len(labels))

        allHandles = []
        allLabels = []
        for i in range(0, len(labels)):
            xes = [x[0] for x in vDataSeries[i]]
            yes = [x[1] for x in vDataSeries[i]]

            handle = ax.plot(xes, yes, color=colors[i])

            allHandles.append(handle)
            allLabels.append(labels[i])

        ax.set_xlabel( xlabel )
        ax.set_ylabel( ylabel )
        ax.set_title(title)

        xTicks = ax.get_xticks().tolist()
        newLabels = [formatter.format(x) for x in xTicks]

        plt.xticks( xTicks, newLabels )

        if not pltcfg.usesMPLD3():
            allHandles = None # fix some mpl error

        cls.makeLegend(fig, ax, allHandles, allLabels, pltcfg, bbox_to_anchor=(0.5, -0.1))

        pltcfg.makePlot()


    @classmethod
    def plotBarsNoHierarchy(cls, plotData, title, xlabel, ylabel, xlabelrotation='horizontal', pltcfg=PlotConfig(), noBarLabels=False, gridLines=True):

        def autolabel(rects):
            """
            Attach a text label above each bar displaying its height
            """

            maxHeight = 0
            for rect in rects:
                maxHeight = max( [maxHeight, rect.get_height()])

            for rect in rects:
                height = rect.get_height()

                heightOffset = min(1.05*height - height, 10)

                if maxHeight <= 1.0:
                    ax.text(rect.get_x() + rect.get_width() / 2., height + heightOffset, '%0.2f' % height,
                            ha='center', va='bottom')
                else:
                    ax.text(rect.get_x() + rect.get_width() / 2., height+heightOffset, '%d' % int(height),ha='center', va='bottom')

        pltcfg.startPlot()
        fig, ax = plt.subplots()

        createdAxes = {}
        colorVector = PorePlot.getColorVector( len(plotData) )

        runIDs = sorted([x for x in plotData])

        runCount = 0

        width = 0.9
        offset = -width/2.0

        mpld3Rects = []
        for run in runIDs:

            data = plotData[run]

            rects = ax.bar(runCount+1, data, width, color=colorVector[runCount], label=run)
            createdAxes[run] = rects
            mpld3Rects.append(rects[0])

            runCount += 1

        # add some text for labels, title and axes ticks
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlabel(xlabel)

        ax.set_xlim([0, runCount+1])

        ax.set_xticks([i for i in range(1, len(runIDs)+1)])
        ax.set_xticklabels(runIDs, rotation=xlabelrotation)

        ax.grid(gridLines)

        allRects = []
        allLabels = []

        for (run, rect) in createdAxes.items():
            allRects.append(rect)
            allLabels.append(run)

        if not noBarLabels:
            if pltcfg.usesMPLD3():

                def makeHTML( label ):
                    runID = label
                    runCount = plotData[label]

                    runIDLine = "<tr><th>RunID</th><th>{rid}</th></tr>".format(rid=runID)
                    countLine = "<tr><td>Read Count</td><td>{rc}</td></tr>".format(rc=runCount)

                    htmlStr = "<table class='tooltip'>" + runIDLine + countLine + "</table>"

                    return htmlStr

                tooltip = MultiAxesPointHTMLTooltip(mpld3Rects, [[makeHTML(x)] for x in allLabels], voffset=10, hoffset=10, css=cls.getToolTipCSS())
                plugins.connect(fig, tooltip)
            else:
                for rect in allRects:
                    autolabel(rect)

        cls.makeLegend(fig, ax, allRects, allLabels, pltcfg)

        pltcfg.makePlot()

    @classmethod
    def plotBars(cls, plotData, title, xlabel, ylabel, xlabelrotation='horizontal', pltcfg=PlotConfig(), noBarLabels=False, gridLines=True):

        def autolabel(rects):
            """
            Attach a text label above each bar displaying its height
            """

            maxHeight = 0
            for rect in rects:
                maxHeight = max( [maxHeight, rect.get_height()])

            for rect in rects:
                height = rect.get_height()

                heightOffset = min(1.05*height - height, 10)

                if maxHeight <= 1.0:
                    ax.text(rect.get_x() + rect.get_width() / 2., height + heightOffset, '%0.2f' % height,
                            ha='center', va='bottom')
                else:
                    ax.text(rect.get_x() + rect.get_width() / 2., height+heightOffset, '%d' % int(height),ha='center', va='bottom')

        allRuns = []
        allGroups = OrderedSet()
        for x in sorted([run for run in plotData]):

            allRuns.append(x)

            for y in plotData[x]:
                allGroups.add(y)

        allGroups = list(allGroups)

        barWidth = 0.9

        width = barWidth / len(allRuns)  # the width of the bars

        ind = np.arange( len(allGroups) )  # the x locations for the groups

        pltcfg.startPlot()
        fig, ax = plt.subplots()

        createdAxes = {}
        colorVector = PorePlot.getColorVector( len(allRuns) )

        runCount = 0

        mpld3Rects = []
        mpld3HTMLs = []

        def makeHTML(label):
            runID = label
            runCounts = plotData[label]

            runIDLine = "<tr><th>RunID</th><th>{rid}</th></tr>".format(rid=runID)
            groupLines = ""
            totalCounts = 0

            for group in allGroups:
                countLine = "<tr><td>" + group + "</td><td>{rc}</td></tr>".format(rc=runCounts[group])
                groupLines += countLine
                totalCounts += runCounts[group]

            groupLines += "<tr><td>Total</td><td>{rc}</td></tr>".format(rc=totalCounts)

            htmlStr = "<table class='tooltip'>" + runIDLine + groupLines + "</table>"

            return htmlStr

        for runs in plotData:

            data = []

            for i in ind:
                data.append(plotData[runs][allGroups[i]])

            lpos = ind + 0.5 - (barWidth/2.0) + 0.5 * width + runCount * width

            rects = ax.bar(lpos, data, width, color=colorVector[runCount], label=runs)
            createdAxes[runs] = rects

            for i in range(0, len(rects)):
                mpld3Rects.append(rects[i])
                mpld3HTMLs.append( [makeHTML(runs)] )

            runCount += 1

        # add some text for labels, title and axes ticks
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlabel(xlabel)

        ax.set_xticks(ind + 0.5)
        ax.set_xticklabels(allGroups, rotation=xlabelrotation)

        ax.grid(gridLines)

        allRects = []
        allLabels = []

        for (run, rect) in createdAxes.items():
            allRects.append(rect)
            allLabels.append(run)

        if not noBarLabels:
            if pltcfg.usesMPLD3():
                tooltip = MultiAxesPointHTMLTooltip(mpld3Rects, mpld3HTMLs, voffset=10,
                                                    hoffset=10, css=cls.getToolTipCSS())
                plugins.connect(fig, tooltip)
            else:
                for rect in allRects:
                    autolabel(rect)

        cls.makeLegend(fig, ax, allRects, allLabels, pltcfg)

        pltcfg.makePlot()

    @classmethod
    def makeLegend(cls, fig, ax, handles, labels, pltcfg, bbox_to_anchor=(0.5, -0.1)):

        if pltcfg.usesMPLD3():
            interactive_legend = plugins.InteractiveLegendPlugin(handles, labels, alpha_unsel=0.5, alpha_over=1.5, start_visible=True)
            plugins.connect(fig, interactive_legend)
        else:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

            if handles == None:
                ax.legend(labels, loc='upper center', bbox_to_anchor=bbox_to_anchor, fancybox=True, shadow=True, ncol=1)
            else:
                ax.legend(handles, labels, loc='upper center', bbox_to_anchor=bbox_to_anchor, fancybox=True,
                          shadow=True, ncol=1)


    @classmethod
    def heat_map_cluster(cls, oHeatMap, vXLabels, vYLabels, sTitle, sLegendText = None,
                         oRange=None, bCluster = True, oHeatMapValues = None, bIsNumber = True, iRoundTo = 3, figSize=(30,30),
                         sXDescr = "", sYDescr = "", iXAxisFontSize = 12, iYAxisFontSize=12, iTitleFontSize=12, iLegendSize=25,
                         pltcfg=PlotConfig()):

        leftbottom = (0.075,0.225)
        widthheight = (0.75,0.75)

        # Generate random features and distance matrix.
        D = oHeatMap
        Dt = np.array(oHeatMap)
        Dt = Dt.transpose()

        pltcfg.startPlot()
        fig = plt.figure(figsize=figSize)

        if len(vYLabels) > 1 and bCluster:
            # Compute and plot first dendrogram.
            ax1 = fig.add_axes([0.85,leftbottom[1],0.15,widthheight[1]])
            Y = sch.linkage(D, method='centroid')
            Z1 = sch.dendrogram(Y, orientation='right', labels = vYLabels)
            ax1.set_xticks([])
            ax1.set_yticks([])

            idx1 = Z1['leaves']
        else:
            idx1 = [x for x in range(0, len(vYLabels))]

        if len(vXLabels) > 1 and bCluster:
            # Compute and plot second dendrogram.
            ax2 = fig.add_axes([leftbottom[0],0.05,widthheight[0],0.15])
            Y = sch.linkage(Dt, method='centroid')
            Z2 = sch.dendrogram(Y, orientation='bottom', labels = vXLabels, leaf_rotation=90)
            ax2.set_xticks([])
            ax2.set_yticks([])

            idx2 = Z2['leaves']

        else:
            idx2 = [x for x in range(0, len(vXLabels))]

        # Plot distance matrix.
        axmatrix = fig.add_axes([leftbottom[0], leftbottom[1],widthheight[0], widthheight[1]])


        D = D[idx1,:]
        D = D[:,idx2]

        cmap = cls.getColorMap('viridis')
        cmap.set_bad(alpha=0.0)
        cmap.set_over(alpha=0.0)
        cmap.set_under(alpha=0.0)

        if oRange is None:
            oRange = (np.min(D), np.max(D))

        im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=cmap, vmin=oRange[0], vmax=oRange[1])

        if not oHeatMapValues is None:

            (xRange, yRange) = D.shape

            oHeatMapValues = oHeatMapValues[idx1,:]
            oHeatMapValues = oHeatMapValues[:,idx2]

            for x in range(0, xRange):
                for y in range(0, yRange):
                    if bIsNumber:

                        if oHeatMapValues[x,y] == int(oHeatMapValues[x,y]):
                            axmatrix.text(y, x, str(int(oHeatMapValues[x,y])), va='center', ha='center')
                        else:
                            axmatrix.text(y, x, str(round(oHeatMapValues[x,y], iRoundTo)), va='center', ha='center')
                    else:
                        axmatrix.text(y, x, str(oHeatMapValues[x,y]).replace(' ','\n'), va='center', ha='center')

        axmatrix.yaxis.set_label_position('left')
        axmatrix.xaxis.set_label_position('top')

        axmatrix.yaxis.set_ticks_position('right')
        axmatrix.xaxis.set_ticks_position('bottom')

        axmatrix.xaxis.set_label_text(sXDescr, fontdict={'size': iTitleFontSize})
        axmatrix.yaxis.set_label_text(sYDescr, fontdict={'size': iTitleFontSize})

        axmatrix.set_xticks(range(len(idx2)))
        axmatrix.set_yticks(range(len(idx1)))

        axmatrix.set_xticklabels( [vXLabels[i] for i in idx2] )
        axmatrix.set_yticklabels( [vYLabels[i] for i in idx1] )

        ax = im.axes
        for item in ([ax.xaxis.label] + ax.get_xticklabels()):
            item.set_fontsize( iXAxisFontSize )

        for item in ([ax.yaxis.label] + ax.get_yticklabels()):
            item.set_fontsize( iYAxisFontSize )

        for item in [ax.title]:
            item.set_fontsize( iTitleFontSize )

        for label in im.axes.xaxis.get_ticklabels():
            label.set_rotation(90)

        # Plot colorbar.
        axcolor = fig.add_axes([0.85,0.05,0.02,0.15])
        minValue = np.amin(Dt)
        maxValue = np.amax(Dt)

        #cbar = mpl.colorbar.ColorbarBase(axcolor, cmap=cmap,
        #                                norm=mpl.colors.Normalize(vmin=minValue, vmax=maxValue),
        #                                orientation='horizontal')

        #pylab.colorbar(im, cax=axcolor)
        cbar = fig.colorbar(im, cax=axcolor)
        #cbar.set_label(sTitle,size=iLegendSize)

        iStep = float(cbar.vmax - cbar.vmin) / 8.0

        vTicks = []
        for i in np.arange(cbar.vmin, cbar.vmax, iStep):
            vTicks.append(i)

        if int(cbar.vmax) - vTicks[len(vTicks)-1] < vTicks[1]-vTicks[0]:
            vTicks[len(vTicks)-1] = cbar.vmax
        else:
            vTicks.append(cbar.vmax)

        cbar.set_ticks(vTicks)
        cbar.ax.tick_params(labelsize=iLegendSize)

        if sLegendText != None:
            cbar.set_label(sLegendText, size=iLegendSize)

        plt.title( sTitle )

        pltcfg.makePlot(noTightLayout=True)