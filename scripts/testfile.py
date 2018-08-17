from collections import Counter

import matplotlib.pyplot as plt
import numpy as np

import mpld3
from mpld3 import plugins
from porestat.plots.poreplot import PorePlot

from porestat.plots.plotconfig import PlotConfig, PlotSaveTYPE

import matplotlib.colors as colors

pltcfg = PlotConfig()
pltcfg.saveToFile("/mnt/c/Users/mjopp/Desktop/mpltest.html")
pltcfg.setOutputType(PlotSaveTYPE.HTML_STRING)

plotData = {'NOT ALIGNED': Counter({'UNKNOWN': 8637}), 'ALIGNED': Counter({'BASECALL_1D': 15882, 'BASECALL_2D': 5876, 'BASECALL_1D_COMPL': 331}), 'UNALIGNED': Counter({'BASECALL_1D': 1981, 'BASECALL_2D': 1598, 'BASECALL_1D_COMPL': 477})}


PorePlot.plotBars(plotData, title="blubb", xlabel="x", ylabel="y", pltcfg = pltcfg)

pltcfg.mpld3js = "mpld3.v0.3.1.dev1.js"
pltcfg.d3js = "d3.v3.min.js"

pltcfg.prepareHTMLOutput("/mnt/c/Users/mjopp/Desktop/", "test.html", relativeImport=True)
