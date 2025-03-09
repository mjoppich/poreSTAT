from collections import defaultdict

import matplotlib.pyplot as plt

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")

import numpy as np
import argparse

def transformPVal(pval):
    return -np.log10(pval)


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-p', '--pvalue', type=int, default=7)
    parser.add_argument('-n', '--set-name', type=int, default=1)
    parser.add_argument('-d', '--set-description', type=int, nargs='+', default=None, required=False)
    parser.add_argument('-c', '--num-sets', type=int, default=20)
    parser.add_argument('-t', '--tsv', type=argparse.FileType('r'), nargs="+", required=True, help='de files')
    args = parser.parse_args()



    for inputfile in args.tsv:

        """

GO.ID   Description     GeneRatio       BgRatio pvalue  p.adjust        qvalue  geneID  Count
GO:0060326      cell chemotaxis 36/283  271/16876       3.74377091841125e-22    1.36834827067931e-18    9.75350844533457e-19    Pdgfb/Ccl3/Il16/Hspb1/Ccl8/Ccl4/Edn1/Ednrb/Il1b/Vcam1/Pde4b/Stap1/Cxcl3/Cxcl1/F7/Ccl22/Ccl17/Cxcl10/Ccl7/Ccl2/Serpine1/Ccr7/Lpar1/Itga9/Thbs1/Trem3/Trem1/Ccrl2/Cxcr4/Ccr2/Gpr18/Ch25h/Fpr2/Cx3cr1/Cxcl2/Cxcl11     36
GO:0030595      leukocyte chemotaxis    31/283  200/16876       3.21738285191329e-21    5.87976716187153e-18    4.19106450446599e-18    Pdgfb/Ccl3/Il16/Ccl8/Ccl4/Edn1/Ednrb/Il1b/Pde4b/Stap1/Cxcl3/Cxcl1/F7/Ccl22/Ccl17/Cxcl10/Ccl7/Ccl2/Serpine1/Ccr7/Itga9/Thbs1/Trem3/Trem1/Ccr2/Gpr18/Ch25h/Fpr2/Cx3cr1/Cxcl2/Cxcl11   31
GO:0050900      leukocyte migration     37/283  331/16876       3.98295987362438e-20    4.8525727793657e-17     3.45888620604222e-17    Pdgfb/Ccl3/Il16/Ccl8/Ccl4/Edn1/Ednrb/Tnf/Il1b/Il1a/Vcam1/Pde4b/Stap1/Cxcl3/Cxcl1/Itgal/F7/Ccl22/Ccl17/Cxcl10/Ccl7/Ccl2/Serpine1/Ccr7/F11r/Itga9/Thbs1/Trem3/Trem1/Ccr2/Gpr18/Ch25h/Spn/Fpr2/Cx3cr1/Cxcl2/Cxcl11     37
GO:0071222      cellular response to lipopolysaccharide 27/283  170/16876       6.81501284123456e-19    6.22721798367808e-16    4.43872546896198e-16    Il12b/Havcr2/Cmpk2/Nos2/Acod1/Tnf/Il6/Il1b/Il1a/Pde4b/Stap1/Cxcl3/Cxcl1/Plscr4/Malt1/Nlrp3/Cxcl10/Ccl2/Serpine1/Cd14/Cx3cr1/Raet1c/Cxcl2/Cxcl11/Irak2/Cd80/Gbp10    27
GO:0032496      response to lipopolysaccharide  31/283  246/16876       1.61735828301614e-18    1.1822889048848e-15     8.42728789571567e-16    Il12b/Irak3/Havcr2/Cmpk2/Nos2/Ednrb/Acod1/Tnf/Il6/Il1b/Il1a/Pde4b/Stap1/Cxcl3/Cxcl1/Plscr4/Ptgs2/Malt1/Nlrp3/Cxcl10/Ccl2/Serpine1/Ccr7/Cd14/Cx3cr1/Raet1c/Cxcl2/Cxcl11/Irak2/Cd80/Gbp10     31
GO:0071219      cellular response to molecule of bacterial origin       27/283  179/16876       2.69537104194482e-18    1.64193019305139e-15    1.1703584787392e-15     Il12b/Havcr2/Cmpk2/Nos2/Acod1/Tnf/Il6/Il1b/Il1a/Pde4b/Stap1/Cxcl3/Cxcl1/Plscr4/Malt1/Nlrp3/Cxcl10/Ccl2/Serpine1/Cd14/Cx3cr1/Raet1c/Cxcl2/Cxcl11/Irak2/Cd80/Gbp10    27
GO:0097529      myeloid leukocyte migration     28/283  199/16876       4.1439555434795e-18     2.16373678734537e-15    1.54229924362583e-15    Pdgfb/Ccl3/Ccl8/Ccl4/Edn1/Ednrb/Il1b/Il1a/Pde4b/Stap1/Cxcl3/Cxcl1/Ccl22/Ccl17/Cxcl10/Ccl7/Ccl2/Serpine1/Ccr7/Itga9/Thbs1/Trem3/Trem1/Ccr2/Fpr2/Cx3cr1/Cxcl2/Cxcl11  28
GO:0002237      response to molecule of bacterial origin        31/283  266/16876       1.5858735283255e-17     7.24545968253713e-15    5.1645223455337e-15     Il12b/Irak3/Havcr2/Cmpk2/Nos2/Ednrb/Acod1/Tnf/Il6/Il1b/Il1a/Pde4b/Stap1/Cxcl3/Cxcl1/Plscr4/Ptgs2/Malt1/Nlrp3/Cxcl10/Ccl2/Serpine1/Ccr7/Cd14/Cx3cr1/Raet1c/Cxcl2/Cxcl11/Irak2/Cd80/Gbp10     31
GO:0071216      cellular response to biotic stimulus    27/283  203/16876       7.28357648766237e-17    2.95794134026733e-14    2.10840372011279e-14    Il12b/Havcr2/Cmpk2/Nos2/Acod1/Tnf/Il6/Il1b/Il1a/Pde4b/Stap1/Cxcl3/Cxcl1/Plscr4/Malt1/Nlrp3/Cxcl10/Ccl2/Serpine1/Cd14/Cx3cr1/Raet1c/Cxcl2/Cxcl11/Irak2/Cd80/Gbp10    27

        """

        colnames = None
        
        allResults = list()
        for lidx, line in enumerate(inputfile):

            if lidx == 0:
                colnames = line.strip().split("\t")
                continue

            line = line.strip().split("\t")

            elemName = line[args.set_name]

            if not args.set_description is None:
                argStrs = []
                for argIdx in args.set_description:
                    argStr = "{}={}".format(colnames[argIdx], line[argIdx])
                    argStrs.append(argStr)

                elemName = "{}\n({})".format(elemName, ", ".join(argStrs))

            pval = transformPVal(float(line[ args.pvalue ]))

            allResults.append((elemName, pval))

        allResults = sorted(allResults, key=lambda x: x[1], reverse=True)

        if len(allResults) > args.num_sets:
            allResults = allResults[:args.num_sets]

        elemLabels = [x[0] for x in allResults]
        elemValues = [x[1] for x in allResults]

        y_pos = np.arange(len(elemLabels))
        fig, ax = plt.subplots(figsize=(8, 5+len(elemLabels)*0.5))

        ax.barh(y_pos, elemValues, align='center')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(elemLabels)
        ax.invert_yaxis()  # labels read top-to-bottom
        ax.set_xlabel('-log10(adj_pval)')
        ax.set_title('Top enriched Gene-Sets')

        ax.axvline(x=transformPVal(0.05), ymin=y_pos[0], ymax=y_pos[-1], c="r")

        outname = inputfile.name + ".png"
        plt.savefig(outname, bbox_inches="tight")
        print(outname)