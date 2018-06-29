import random, os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

import argparse

from porestat.analysis.alignmentReport import ReportFactory
from porestat.analysis.foldchange_distribution import FoldChangeDistributionFactory
from porestat.analysis.foldchange_similarity import FoldChangeSimilarityFactory

from porestat.analysis.alignmentStatistics import AlignmentStatisticAnalysisFactory
from porestat.analysis.foldchange_top_regulated import FoldChangeTopRegulatedFactory
from porestat.analysis.read_counts import ReadCountAnalysisFactory
from porestat.analysis.similarity_analysis import SimilarityAnalysisFactory

from porestat.tools.PTToolInterface import PSToolException
from porestat.utils import eprint

import random


def argParseID(elem):

    return elem

if __name__ == '__main__':

    sys.stderr.write("Enrichment analysis requires the following packages to be installed: {}\n".format(" ".join(['libcurl4-openssl-dev', 'libssl-dev', 'libmariadbclient-dev', 'libcurl4-openssl-dev', 'libxml2-dev'])))


    parser = argparse.ArgumentParser(prog='poreAnalysis')

    parser.register('type', None, argParseID)

    #parser.add_argument('command', type=str, help='foo help')
    subparsers = parser.add_subparsers(help='sub-command help')

    cmd2tool = {}

    allTools = [
        ReadCountAnalysisFactory(parser, subparsers, 'read_counts'),
        AlignmentStatisticAnalysisFactory(parser, subparsers, 'alignment_stat'),
        SimilarityAnalysisFactory(parser, subparsers, 'similarity'),
        ReportFactory(parser, subparsers, 'summary'),
        FoldChangeDistributionFactory(parser, subparsers, 'foldchange'),
        ReportFactory(parser, subparsers, 'timeline'),
        FoldChangeSimilarityFactory(parser, subparsers, 'fcdist'),
        FoldChangeTopRegulatedFactory(parser, subparsers, 'topreg')
    ]

    for tool in allTools:
        cmd2tool[tool.which] = tool

    args = parser.parse_args()

    if not 'func' in vars(args):
        parser.print_help()
    else:

        try:
            calcObj = args.func(args)
            calcObj.exec()

        except PSToolException as e:
            eprint(e)

            if args.which:
                cmd2tool[args.which].print_usage()
