import sys

#import porestat.utils.PickleArgparse as pickArg
import argparse

from porestat.tools.experiment_ls import ExperimentLsFactory
from porestat.tools.channel_occupancy import ChannelOccupancyFactory
from porestat.tools.nucleotide_distribution import NucleotideDistributionFactory
from porestat.tools.quality_distribution import QualityDistributionFactory
from porestat.tools.stats_summary import StatsSummaryFactory
from porestat.tools.extract_sequences import ExtractSequencesFactory
from porestat.tools.timelineReads import TimelineReadsFactory
from porestat.tools.quality_position import QualityPositionFactory

from porestat.tools.PTToolInterface import PSToolException
from porestat.utils import eprint

import random


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
    cmd2tool['STAT'] = StatsSummaryFactory(parser, subparsers)
    cmd2tool['FASTQ'] = ExtractSequencesFactory(parser, subparsers)
    cmd2tool['TIME'] = TimelineReadsFactory(parser, subparsers)
    cmd2tool['NUC_DIST'] = NucleotideDistributionFactory(parser, subparsers)
    cmd2tool['QUAL_DIST'] = QualityDistributionFactory(parser, subparsers)
    cmd2tool['QUAL_POS'] = QualityPositionFactory(parser, subparsers)

    args = parser.parse_args()

    if not 'func' in vars(args):
        parser.print_help()
    else:

        try:
            calcObj = args.func(args)

            calcObj.exec()

        except PSToolException as e:
            eprint(e)

