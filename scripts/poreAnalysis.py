import sys

#import porestat.utils.PickleArgparse as pickArg
import argparse

from porestat.analysis.read_counts import ReadCountAnalysisFactory

from porestat.tools.PTToolInterface import PSToolException
from porestat.utils import eprint

import random


def argParseID(elem):

    return elem

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='poreAnalysis')

    parser.register('type', None, argParseID)

    #parser.add_argument('command', type=str, help='foo help')
    subparsers = parser.add_subparsers(help='sub-command help')

    cmd2tool = {}
    cmd2tool['read_counts'] = ReadCountAnalysisFactory(parser, subparsers)
    cmd2tool['alignment_stat'] = AlignmentStatisticsFactory(parser, subparsers)


    args = parser.parse_args()

    if not 'func' in vars(args):
        parser.print_help()
    else:

        try:
            calcObj = args.func(args)

            calcObj.exec()

        except PSToolException as e:
            eprint(e)
