import sys
import argparse

from porestat.tools.experiment_ls import Experiment_ls
from porestat.tools.channel_occupancy import Channel_occupancy
from porestat.tools.stats_summary import StatsSummary
from porestat.tools.extract_sequences import ExtractSequences

import random

if __name__ == '__main__':
    from porestat.plots.poreplot import PorePlot

    counts = {}
    for x in range(1,513):

        counts[x] = []

        for j in range(0, random.randint(5,10)):
            counts[x].append(  random.randint(1000,60000) )

        if random.randint(0, 100) % 10 == 0:
            counts[x] = []

    #PorePlot.plotLoadOut(counts, pores=(32,16))
    #exit(0)

    parser = argparse.ArgumentParser(prog='porestat')
    #parser.add_argument('command', type=str, help='foo help')
    subparsers = parser.add_subparsers(help='sub-command help')

    cmd2tool = {}
    cmd2tool['EXPLS'] = Experiment_ls(parser, subparsers)
    cmd2tool['OCC'] = Channel_occupancy(parser, subparsers)
    cmd2tool['STAT'] = StatsSummary(parser, subparsers)
    cmd2tool['FASTQ'] = ExtractSequences(parser, subparsers)

    cmd2tool['TIME'] = TimelineReads(parser, subparsers)

    args = parser.parse_args()

    if not 'func' in vars(args):
        parser.print_help()
    else:

        args.func(args)
