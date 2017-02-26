import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

class PorePlot:

    @classmethod
    def getColor(cls, colormap="Viridis", value=0.5):
        cmap = plt.cm.get_cmap(colormap)

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

        for i in pore2length:
            p2c[i] = 0

            minPore = min(i, minPore)

            if i in pore2length:
                count = len(pore2length[i])
                p2c[i] = count

                average =np.average(pore2length[i])
                p2l[i] = average
                minAvgLength = min(minAvgLength, average)
                maxAvgLength = max(maxAvgLength, average)

                minCount = min(count, minCount)
                maxCount = max(count, maxCount)

            else:
                p2c[i] = 0
                minCount = 0

        minPoreRad = 50
        maxPoreRad = 200

        p2area = {}

        for x in p2c:

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

            frac = p2l[x] / (maxAvgLength - minAvgLength)
            color.append(cls.getColor(colormap='plasma', value=frac))




        fig=plt.figure(figsize=(10, 10))
        fig.set_dpi(100)

        plt.Axes(fig, [1, pores[0], 1, pores[1]])

        print(area)

        plt.scatter(xvec, yvec, s=area, c=color, alpha=0.5)
        plt.axis('off')

        plt.show()