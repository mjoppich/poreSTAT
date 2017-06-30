import argparse
import HTSeq
import matplotlib

from porestat.plots.poreplot import PorePlot

from porestat.plots.plotconfig import PlotConfig
from ..utils.DataFrame import DataFrame, DataRow
from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory,PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Files import fileExists
from scipy import stats
import sys
import numpy as np
from matplotlib import pyplot as plt
import mpld3
import random
import os


class SimilarityAnalysisFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(SimilarityAnalysisFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser = subparsers.add_parser('similarity', help='expls help')
        parser.add_argument('-c', '--counts', nargs='+', type=str, help='counts summary file', required=False)

        def fileOpener( filename ):
            open(filename, 'w').close()
            return filename

        parser.add_argument('-o', '--output', type=fileOpener, help='output location, default: std out', default=sys.stdout)

        parser = PlotConfig.addParserArgs(parser)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return SimilarityAnalysis(simArgs)



class SimilarityAnalysis(ParallelPSTInterface):

    def __init__(self, args):

        super(SimilarityAnalysis, self).__init__(args)
        self.counts = None

    def _makePropDict(self):

        return None

    def readCounts(self, args):

        for countsFile in args.counts:

            if not fileExists(countsFile):
                raise PSToolException("Read info file does not exist: " + str(countsFile))

        self.counts = {}

        for countsFile in args.counts:
            self.counts[countsFile] = DataFrame.parseFromFile(countsFile)

    def prepareInputs(self, args):
        self.readCounts(args)

        self.writeLinesToOutput(args.output, "\t".join(['gene', 'coverage', 'rank']) + "\n", mode='w')

        return []

    def execParallel(self, data, environment):

        return None


    def joinParallel(self, existResult, newResult, oEnvironment):

        return None


    def makeResults(self, parallelResult, oEnvironment, args):

        genes = {}
        unionGenes = set()

        for x in self.counts:
            df = self.counts[x]
            allgenes = df['gene']

            dfGenes = allgenes.to_set()

            unionGenes = unionGenes.union(dfGenes)

            genes[x] = dfGenes

        N = len(self.counts)
        unionGenes = list(unionGenes)

        dataRows = [x for x in self.counts]

        cntData = np.zeros((N,len(unionGenes)))

        for i in range(0, N):

            expName = dataRows[i]
            df = self.counts[expName]

            for j in range(0,len(unionGenes)):

                genename = unionGenes[j]
                geneExpr = df.findRow('gene', genename)

                if geneExpr == None:
                    cntData[i, j] = 0.0
                else:
                    cntData[i, j] = geneExpr['coverage']


        tooltipCSS = """
        table
        {
          border-collapse: collapse;
        }
        th
        {
          color: #ffffff;
          background-color: #000000;
        }
        td
        {
          background-color: #cccccc;
        }
        table, th, td
        {
          font-family:Arial, Helvetica, sans-serif;
          border: 1px solid black;
          text-align: right;
        }
        """

        def makeHTMLTable( gene, cond1, cond1val, cond2, cond2val):

            cond1 = cond1.split( os.path.sep )[-2]
            cond2 = cond2.split(os.path.sep)[-2]

            lines = []
            lines.append("<tr><th></th><th>{gn}</th></tr>".format(gn=gene))
            lines.append( "<tr><td>{c}</td><td>{cv:.2f}</td></tr>".format(c=cond1, cv=cond1val) )
            lines.append( "<tr><td>{c}</td><td>{cv:.2f}</td></tr>".format(c=cond2, cv=cond2val))



            htmlStr = "<table>"+"".join(lines)+"</table>"

            return htmlStr

        dataShape = cntData.shape

        fig, ax = plt.subplots(dataShape[0], dataShape[0], sharex="col", sharey="row", figsize=(8, 8))
        fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95,
                            hspace=0.1, wspace=0.1)

        allPltAnnot = []
        for i in range(dataShape[0]):
            for j in range(dataShape[0]):
                xdata = []
                ydata = []
                colors = []

                X = cntData[j, :]
                Y = cntData[i, :]

                htmlData = []
                for k in range(0, len(unionGenes)):
                    #xdata.append(x[0])
                    #ydata.append(x[1])

                    if i != j and X[k] == Y[k] and X[k] > 0:
                        print(unionGenes[k])

                    htmlData.append( makeHTMLTable(unionGenes[k], dataRows[j], X[k], dataRows[i], Y[k]) )
                    colors.append( 'blue' )

                print(X)
                print(Y)

                points = ax[i,j].scatter(X, Y, c=colors, s=40, alpha=0.6)
                allPltAnnot.append((points, htmlData))

        combinedPoints = []
        combinedHTML = []

        for tooltipdata in allPltAnnot:
            combinedPoints.append(tooltipdata[0])
            combinedHTML.append(tooltipdata[1])

        # remove tick labels
        #for axi in ax.flat:
        #    for axis in [axi.xaxis, axi.yaxis]:
        #        axis.set_major_formatter(plt.NullFormatter())

        # Here we connect the linked brush plugin
        tooltip = MultiAxesPointHTMLTooltip( combinedPoints, combinedHTML, voffset=10, hoffset=10, css=tooltipCSS)
        mpld3.plugins.connect(fig, tooltip)
        mpld3.plugins.connect(fig, mpld3.plugins.MousePosition(fontsize=14))

        mpld3.plugins.connect(fig, mpld3.plugins.LinkedBrush(points))
        mpld3.show()




        return

    def get_id(self, obj, suffix="", prefix="el", warn_on_invalid=True):
        """Get a unique id for the object"""
        if not suffix:
            suffix = ""
        if not prefix:
            prefix = ""
        import os

        objid = prefix + str(os.getpid()) + str(id(obj)) + suffix

        if warn_on_invalid and not self.html_id_ok(objid):
            print('"{0}" is not a valid html ID. This may cause problems')

        return objid

    def html_id_ok(self, objid, html5=False):
        """Check whether objid is valid as an HTML id attribute.

        If html5 == True, then use the more liberal html5 rules.
        """
        import re

        if html5:
            return not re.search('\s', objid)
        else:
            return bool(re.match("^[a-zA-Z][a-zA-Z0-9\-\.\:\_]*$", objid))


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

       console.log(labels);


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

        ids = []
        for point in points:

            if isinstance(points, matplotlib.lines.Line2D):
                suffix = "pts"
            else:
                suffix = None

            points_id = mpld3.utils.get_id(point, suffix)
            ids.append(points_id)

        print(ids)

        self.dict_ = {"type": "htmltooltip",
                      # "id": get_id(points, suffix),
                      "ids": ids,
                      "labels": labels,
                      "hoffset": hoffset,
                      "voffset": voffset}







