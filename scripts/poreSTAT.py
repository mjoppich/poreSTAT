import sys

#import porestat.utils.PickleArgparse as pickArg
import argparse

from porestat.tools.demangle_files import DemangleFilesFactory
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

import random, os

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")



def argParseID(elem):

    return elem

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='porestat')

    parser.register('type', None, argParseID)

    #parser.add_argument('command', type=str, help='foo help')
    subparsers = parser.add_subparsers(help='sub-command help')

    cmd2tool = {}
    cmd2tool['EXPLS'] = ExperimentLsFactory(parser, subparsers)
    cmd2tool['OCC'] = ChannelOccupancyFactory(parser, subparsers)
    cmd2tool['REPORT'] = ReportFactory(parser, subparsers)
    cmd2tool['FASTQ'] = ExtractSequencesFactory(parser, subparsers)
    cmd2tool['TIME'] = TimelineReadsFactory(parser, subparsers)
    cmd2tool['NUC_DIST'] = NucleotideDistributionFactory(parser, subparsers)
    cmd2tool['QUAL_DIST'] = QualityDistributionFactory(parser, subparsers)
    cmd2tool['QUAL_POS'] = QualityPositionFactory(parser, subparsers)
    cmd2tool['HIST'] = LengthHistogramFactory(parser, subparsers)
    cmd2tool['YIELD'] = YieldPlotFactory(parser, subparsers)
    cmd2tool['INFO'] = ReadInfoFactory(parser, subparsers)
    cmd2tool['DEMANGLE'] = DemangleFilesFactory(parser, subparsers)

    args = parser.parse_args()


    if not 'func' in vars(args):
        parser.print_help()
    else:

        try:
            calcObj = args.func(args)

            calcObj.exec()

        except PSToolException as e:
            eprint(e)

