import matplotlib.pyplot as plt
import pandas as pd
import argparse
import numpy as np

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")




if __name__ == '__main__':

    #--defile {defile} --enrichment {resfile} --elems {maxelems}
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-d', '--defile', type=argparse.FileType("r"), required=False)
    parser.add_argument('-e', '--enrichment', type=argparse.FileType("r"), required=True)
    parser.add_argument('-c', '--elems', type=int, default=20, required=False)

    args = parser.parse_args()


    enrichmentDF = pd.read_csv(args.enrichment, sep="\t")



    def custom_div_cmap(numcolors=11, name='custom_div_cmap',
                        mincol='blue', midcol='white', maxcol='red'):
        """ Create a custom diverging colormap with three colors
        
        Default is blue to white to red with 11 colors.  Colors can be specified
        in any way understandable by matplotlib.colors.ColorConverter.to_rgb()
        """

        from matplotlib.colors import LinearSegmentedColormap 
        
        cmap = LinearSegmentedColormap.from_list(name=name, 
                                                colors = [x for x in [mincol, midcol, maxcol] if not x is None],
                                                N=numcolors)
        return cmap
    color1 = "#883656"
    color2 = "#4d6841"
    color3 = "#afafaf"
    rvb = custom_div_cmap(150, mincol=color1,maxcol=color2, midcol=color3)

    def plotORAresult( dfin, title, numResults=10, filename=None, less=None, more=None):
        #https://www.programmersought.com/article/8628812000/
        
        def makeTitle(colDescr, colID, colSize, setSize):
            out = []
            for x,y,z, s in zip(colDescr, colID, colSize, setSize):
                out.append("{} ({}, r={}/{})".format(x, y, z, s))

            return out


        df_raw = dfin.sort_values("qvalue").copy()

        # Prepare Data
        #df = df_raw[['cty', 'manufacturer']].groupby('manufacturer').apply(lambda x: x.mean())
        df = pd.DataFrame(df_raw)

        #determine plot type

        termIDColumn = df.columns[0]
        termNameColumn = df.columns[1]
        qvalueColumn = "qvalue"

        df = df_raw[[termNameColumn, qvalueColumn]]

        print(df_raw.columns)
        print("GeneRatio" in df_raw.columns)
        print("BgRatio" in df_raw.columns)

        rvb = None

        if "GeneRatio" in df_raw.columns and "BgRatio" in df_raw.columns:
            #ORA
            df["sig_geneset_size"] = df_raw.GeneRatio.str.split("/").apply(lambda x: x[0]).astype(int)
            df["geneset_size"] = df_raw.BgRatio.str.split("/").apply(lambda x: x[0]).astype(int)

            rvb = custom_div_cmap(150, mincol=color2, maxcol=color3, midcol=None)

            colorValues = [rvb(x/df.geneset_size.max()) for x in df.geneset_size]



        elif "NES" in df_raw.columns:
            #GSEA
            df["sig_geneset_size"] = df_raw.core_enrichment.str.split("/").apply(len)
            df["geneset_size"] = df_raw.setSize

            df["dot_color"] = df_raw["NES"]
            rvb = custom_div_cmap(150, mincol=color1,maxcol=color2, midcol=color3)

            maxNES = max(df_raw.NES.max(), -df_raw.NES.min())
            print(maxNES)
            nesValues = (1.0+(df_raw.NES / maxNES)) / 2.0
            colorValues = [rvb(x) for x in nesValues]


        else:
            raise ValueError()

        df['termtitle'] = makeTitle(df_raw[termNameColumn], df_raw[termIDColumn], df["sig_geneset_size"],df["geneset_size"])


        df.sort_values('qvalue', inplace=True, ascending=True)
        df.reset_index()
        df = df[:numResults]
        colorValues = colorValues[:numResults]
        
        df = df.iloc[::-1]
        colorValues = colorValues[::-1]
        
        maxNLog = max(-np.log(df.qvalue))
        maxLine = ((maxNLog// 10)+1)*10
        #print(maxNLog, maxLine)
        
        
        # Draw plot
        fig, ax = plt.subplots(figsize=(10,10), dpi= 80)
        ax.hlines(y=df.termtitle, xmin=0, xmax=maxLine, color='gray', alpha=0.7, linewidth=1, linestyles='dashdot')
        ax.vlines(x=-np.log(0.05), ymin=0, ymax=numResults, color='red', alpha=0.7, linewidth=1, linestyles='dashdot')
        
        sizeFactor = 10    
        scatter = ax.scatter(y=df.termtitle, x=-np.log(df.qvalue), s=df.geneset_size*sizeFactor, c=colorValues, alpha=0.7, )

        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6, func=lambda x: x/sizeFactor)
        labels = [x for x in labels]

        # Title, Label, Ticks and Ylim
        ax.set_title(title, fontdict={'size':12})
        ax.set_xlabel('Neg. Log. Adj. p-Value')
        ax.set_yticks(df.termtitle)
        ax.set_yticklabels(df.termtitle, fontdict={'horizontalalignment': 'right'})
        plt.grid(b=None)
        plt.tight_layout()
        plt.yticks(fontsize=16)
        
        if filename != None:
            for x in filename:
                plt.savefig(x, bbox_inches='tight')
                print(x)

            plt.close()
        else:
            plt.show()

    df = pd.read_csv(args.enrichment.name, sep="\t")
    print(df.columns)
    plotORAresult(df, "Enrichments", numResults=args.elems, filename=["{}.png".format(args.enrichment.name)])