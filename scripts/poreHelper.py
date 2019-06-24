import random, os, sys


sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")
#import porestat.utils.PickleArgparse as pickArg
import argparse

from porestat.tools.PTToolInterface import PSToolException
from porestat.utils import eprint
from porestat.helper.filter_length import FilterLengthFactory
from porestat.helper.extract_entries_fasta import FAFQRecordExtractFactory
from porestat.helper.sequence_merge import FARecordSequenceMergeFactory
from porestat.helper.revcompl import FaRevComplFactory

from porestat.assemblies.GenomeSimilarity import GenomeSimilarityPlotFactory




def argParseID(elem):

    return elem

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='porestat')

    parser.register('type', None, argParseID)

    #parser.add_argument('command', type=str, help='foo help')
    subparsers = parser.add_subparsers(help='sub-command help')

    cmd2tool = {}

    allTools = [
        FilterLengthFactory(parser, subparsers, 'filter_length'),
        FAFQRecordExtractFactory(parser, subparsers, 'seq_extract'),
        GenomeSimilarityPlotFactory(parser, subparsers, 'simplot'),
        FARecordSequenceMergeFactory(parser, subparsers, 'seq_merge'),
        FaRevComplFactory(parser, subparsers, 'rev_compl')
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

