import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from collections import Counter

class PorePlot:

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
    def plotLoadOut(cls, pore2length, pores=(16,8)):

        p2c = Counter()
        p2l = {}

        minCount = 0
        maxCount = 0
        minPore = pores[0]*pores[1]

        minAvgLength = 1000000
        maxAvgLength = 0

        scolormap = "plasma"

        for i in pore2length:
            p2c[i] = 0

            minPore = min(i, minPore)

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
                #minCount = 0

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

        for x in p2c:

            coords = p2coord[x]
            xvec.append( coords[0] )
            yvec.append( coords[1] )
            area.append( p2area[x] )

            if p2l[x] == -1:
                color.append( (1,1,1,1) )

            else:
                frac = p2l[x] / (maxAvgLength - minAvgLength)
                color.append(cls.getColor(colormap=scolormap, value=frac))

        fig, axarr = plt.subplots(2, figsize=(10, 10), gridspec_kw = {'height_ratios':[10, 1]})
        fig.set_dpi(100)

        axarr[0].set_title('title')
        axarr[0].set_xlabel('xlabel', fontsize=18)
        axarr[0].set_ylabel('ylabel', fontsize=18)

        axarr[0].axes.get_xaxis().set_ticks( [x for x in range(0, pores[1])] )
        axarr[0].axes.get_yaxis().set_ticks( [x for x in range(0, pores[0])])

        axarr[0].axes.get_xaxis().set_ticklabels( [str(x) for x in range(1, pores[1]+1)] )
        axarr[0].axes.get_yaxis().set_ticklabels( [str(x) for x in range(1, pores[0]+1)])

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