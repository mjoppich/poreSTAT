import sys
import argparse

from porestat.tools.experiment_ls import Experiment_ls
from porestat.tools.channel_occupancy import Channel_occupancy
from porestat.plots.poreplot import PorePlot
from porestat.tools.stats_summary import StatsSummary

import random

if __name__ == '__main__':

    random.random()

    pore2count = {}
    for x in range(1, 16*8+1):

        pore2count[x] = []

        for i in range(0, random.randint(2,20)):
            pore2count[x].append( 1000 + random.randint(1000, 15000) )

    PorePlot.plotLoadOut(pore2count)


    parser = argparse.ArgumentParser(prog='porestat')
    #parser.add_argument('command', type=str, help='foo help')
    subparsers = parser.add_subparsers(help='sub-command help')

    cmd2tool = {}
    cmd2tool['EXPLS'] = Experiment_ls(parser, subparsers)
    cmd2tool['OCC'] = Channel_occupancy(parser, subparsers)
    cmd2tool['STAT'] = StatsSummary(parser, subparsers)

    args = parser.parse_args()

    args.func(args)
