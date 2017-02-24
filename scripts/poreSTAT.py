import sys
import argparse

from porestat.tools.experiment_ls import Experiment_ls
from porestat.tools.channel_occupancy import Channel_occupancy

if __name__ == '__main__':


    parser = argparse.ArgumentParser(prog='porestat')
    #parser.add_argument('command', type=str, help='foo help')
    subparsers = parser.add_subparsers(help='sub-command help')

    cmd2tool = {}
    cmd2tool['EXPLS'] = Experiment_ls(parser, subparsers)
    cmd2tool['OCC'] = Channel_occupancy(parser, subparsers)

    args = parser.parse_args()

    args.func(args)
