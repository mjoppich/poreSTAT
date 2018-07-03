import matplotlib.pyplot as plt
import numpy as np

import mpld3
from mpld3 import plugins
from porestat.plots.plotconfig import PlotConfig, PlotSaveTYPE

pltcfg = PlotConfig()
pltcfg.saveToFile("/mnt/c/Users/mjopp/Desktop/mpltest.html")
pltcfg.setOutputType(PlotSaveTYPE.HTML)

fig, ax = plt.subplots()

#Create x,y arrays of normally distributed points
npts = 100000
x = np.random.standard_normal(npts)
y = np.random.standard_normal(npts)

#Set bin numbers in both axes
nxbins = 100
nybins = 100

#Set the cutoff for resolving the individual points
minperbin = 1

#Make the density histrogram
H, yedges, xedges = np.histogram2d(y,x,bins=(nybins,nxbins))
#Reorient the axes
H =  H[::-1]

extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]

#Figure out which bin each x,y point is in
xbinsize = xedges[1]-xedges[0]
ybinsize = yedges[1]-yedges[0]
xi = ((x-xedges[0])/xbinsize).astype(np.integer)
yi = nybins-1-((y-yedges[0])/ybinsize).astype(np.integer)

#Subtract one from any points exactly on the right and upper edges of the region
xim1 = xi-1
yim1 = yi-1
xi = np.where(xi < nxbins,xi,xim1)
yi = np.where(yi < nybins,yi,yim1)

#Get all points with density below the threshold
lowdensityx = x[H[yi,xi] <= minperbin]
lowdensityy = y[H[yi,xi] <= minperbin]

ax.plot(lowdensityx,lowdensityy,linestyle='None',marker='o',mfc='k',mec='k',ms=3)
cp1 = ax.imshow(H,interpolation='nearest',extent=extent,vmin=minperbin)
fig.colorbar(cp1)

ax.set_title('An Image', size=20)

plugins.connect(fig, plugins.MousePosition(fontsize=14))

pltcfg.makePlot(noTightLayout=True)