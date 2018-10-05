import random, os, sys

from porestat.tools.squiggle import SquigglePlotFactory

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

#import porestat.utils.PickleArgparse as pickArg
import argparse

from porestat.tools.kmer_coverage import KmerHistogramFactory


from porestat.tools.experiment_ls import ExperimentLsFactory
from porestat.tools.channel_occupancy import ChannelOccupancyFactory
from porestat.tools.nucleotide_distribution import NucleotideDistributionFactory
from porestat.tools.quality_distribution import QualityDistributionFactory
from porestat.tools.read_info import ReadInfoFactory
from porestat.tools.report import ReportFactory
from porestat.tools.extract_sequences import ExtractSequencesFactory
from porestat.tools.timelineReads import TimelineReadsFactory
from porestat.tools.quality_position import QualityPositionFactory
from porestat.tools.length_histogram import LengthHistogramFactory
from porestat.tools.yield_plot import YieldPlotFactory

from porestat.tools.PTToolInterface import PSToolException
from porestat.utils import eprint
from porestat.tools.demangle_files import DemangleFilesFactory





def argParseID(elem):

    return elem

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='porestat')

    parser.register('type', None, argParseID)

    #parser.add_argument('command', type=str, help='foo help')
    subparsers = parser.add_subparsers(help='sub-command help')

    cmd2tool = {}

    allTools = [
        ReadInfoFactory(parser, subparsers, 'info'),
        ExperimentLsFactory(parser, subparsers, 'expls'),
        ChannelOccupancyFactory(parser, subparsers, 'occupancy'),
        ReportFactory(parser, subparsers, 'summary'),
        ExtractSequencesFactory(parser, subparsers, 'seq'),
        TimelineReadsFactory(parser, subparsers, 'timeline'),
        NucleotideDistributionFactory(parser, subparsers, 'nucleotides'),
        QualityDistributionFactory(parser, subparsers, 'qual_dist'),
        QualityPositionFactory(parser, subparsers, 'qual'),
        LengthHistogramFactory(parser, subparsers, 'histo'),
        YieldPlotFactory(parser, subparsers, 'yield'),
        DemangleFilesFactory(parser, subparsers, 'demangle'),
        KmerHistogramFactory(parser, subparsers, 'kmer'),
        SquigglePlotFactory(parser, subparsers, 'squiggle')

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

