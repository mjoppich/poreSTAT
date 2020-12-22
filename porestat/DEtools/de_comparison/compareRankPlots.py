# -*- coding: utf-8 -*-
# python3
# Copyright (c) 2017 by Dr. Justin Klotz
import argparse
from collections import defaultdict

import numpy as np
from typing import List
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors as mpl_colors, cm as mpl_cm, colorbar as mpl_colorbar, figure as mpl_figure, \
    ticker as mpl_ticker
import random, sys, matplotlib.colors
# Make sure that we are using QT5
import pandas as pd
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE


# TODO: Verify set_visible works
# TODO: Verify re-plot works
# TODO: Determine which part of code takes longest (pproject from pypi)


class ParCoord:
    def __init__(self,
                 data_sets: list,
                 figsize=(8, 6)):

        # Verify user input
        num_dims = len(data_sets[0])
        # dimensions
        if num_dims < 2:
            raise ValueError('Must supply data with more than one dimension.')

        # Setup figure and axes
        x = range(num_dims)
        fig = plt.figure(figsize=figsize)
        axes = list()
        for idx in range(num_dims - 1):
            axes.append(fig.add_subplot(1, num_dims - 1, idx + 1))

        # Store, initialize values
        self.fig = fig
        self._axes = axes
        self._ax_color_bar = None
        self._x = x
        self._data_sets = None
        self._data_sets_info = None
        self._num_dims = num_dims
        self._colors = None
        self._scores = None
        self._scores_norm_min = None
        self._scores_norm_max = None
        self._color_style = None
        self._color_map_norm = None
        self._use_variable_line_width = False
        self._sorted_scores_indices = None

        self._set_data(data_sets)

    def _set_data(self,
                  data_sets: list):
        # make sure there are no nan's
        for data_set in data_sets:
            if np.nan in data_set:
                raise ValueError('Argument "data_sets" must not contain nan\'s')

        self._data_sets = data_sets

    def reset_data(self,
                   data_sets: list):
        self._scores = None
        self._sorted_scores_indices = None
        self._colors = None
        self._set_data(data_sets)

    def set_visible(self,
                    visible: List[bool]):
        if len(visible) != len(self._data_sets):
            raise ValueError('List "visible" passed to set_visible must be same length as number of data sets')
        if self._scores is not None:
            visible = [visible[idx] for idx in self._sorted_scores_indices]
        for ax in self._axes:
            for idx, line in enumerate(ax.lines):
                if not isinstance(visible[idx], bool):
                    raise ValueError('List "visible" passed to set_visible must only contain booleans.')
                line.set_visible(visible[idx])

    def save_fig(self,
                 path: str):
        try:
            self.fig.savefig(path)
        except:
            raise ValueError('Could not save file to path ' + path)

    def plot(self,
             num_ticks: int = 7,
             line_width=1.1,
             y_min: list or None = None,
             y_max: list or None = None,
             preRanges=None,
             labelFunc=lambda colIdx, val: '{:4.2f}'.format(val)
             ):

        # Calculate data set limits
        data_sets_info = list()
        for idx, m in enumerate(zip(*self._data_sets)):
            if y_min is not None:
                mn = y_min[idx]
            else:
                mn = min(m)
            if y_max is not None:
                mx = y_max[idx]
            else:
                mx = max(m)
            if mn == mx:
                mn -= 0.5
                mx = mn + 1.
            r = float(mx - mn)

            if preRanges == None or not idx in preRanges:
                data_sets_info.append({'min': mn, 'max': mx, 'range': r})
            else:
                data_sets_info.append(preRanges[idx])

        # Normalize the data sets
        norm_data_sets = list()
        for ds in self._data_sets:
            nds = [(value - data_sets_info[dimension]['min']) / data_sets_info[dimension]['range']
                   for dimension, value in enumerate(ds)]
            norm_data_sets.append(nds)

        # Clear axes in case they were plotted on previously
        for ax in self._axes:
            ax.clear()

        # Plot the data sets on all the subplots
        for i, ax in enumerate(self._axes[:self._num_dims]):
            for dsi in range(len(norm_data_sets)):
                if self._colors is not None:
                    if self._scores is not None and self._use_variable_line_width and self._scores[0] != self._scores[
                        -1]:
                        min_width = line_width
                        max_width = min_width + 3
                        line_width_iter = min_width + (max_width - min_width) * (self._scores[-1] - self._scores[dsi]) \
                                          / (self._scores[-1] - self._scores[0])
                    else:
                        line_width_iter = line_width
                    ax.plot(self._x, norm_data_sets[dsi], linewidth=line_width_iter, color=self._colors[dsi])
                else:
                    ax.plot(self._x, norm_data_sets[dsi], linewidth=line_width)
            ax.set_xlim([self._x[i], self._x[i + 1]])

        # Set the axis ticks
        for dimension, (ax, xx) in enumerate(zip(self._axes[:self._num_dims], self._x[:self._num_dims])):
            # x-axis
            ax.xaxis.set_major_locator(mpl_ticker.FixedLocator([xx]))
            # y-axis
            ax.yaxis.set_major_locator(mpl_ticker.FixedLocator(list(np.linspace(0, 1, num_ticks))))
            min_ = data_sets_info[dimension]['min']
            max_ = data_sets_info[dimension]['max']
            labels = [labelFunc(dimension, val) for val in np.linspace(min_, max_, num_ticks)]
            ax.set_yticklabels(labels, color='k', weight='semibold')  # backgroundcolor='0.75'

        # Move the final axis' ticks to the right-hand side
        if len(self._axes) < self._num_dims:
            ax_last = self._axes[-1].twinx()  # create last axis if this is not a re-plot
        else:
            ax_last = self._axes[-1]
        dimension_last = self._num_dims - 1
        # x-axis
        ax_last.xaxis.set_major_locator(mpl_ticker.FixedLocator([self._x[-2], self._x[-1]]))
        # y-axis
        ax_last.yaxis.set_major_locator(mpl_ticker.FixedLocator(list(np.linspace(0, 1, num_ticks))))
        num_ticks = len(ax_last.get_yticklabels())
        min_ = data_sets_info[dimension_last]['min']
        max_ = data_sets_info[dimension_last]['max']
        labels = [labelFunc(dimension_last, val) for val in np.linspace(min_, max_, num_ticks)]
        ax_last.set_yticklabels(labels, color='k', weight='semibold')

        # Stack the subplots
        self.fig.subplots_adjust(wspace=0)

        # Adjust plot borders, labels
        self._axes.append(ax_last)
        for ax in self._axes:
            ax.tick_params(direction='inout', length=10, width=1)
            ax.tick_params(axis='x', pad=20, labelsize=12)
            ax.set_ylim(0, 1)
            ax.spines['bottom'].set_visible(False)
            ax.spines['top'].set_visible(False)

    def set_colors(self,
                   colors: str or list):
        if isinstance(colors, str):
            colors = [colors] * len(self._data_sets)
        self._colors = colors
        self._scores = None

    def set_scores(self,
                   scores: list,
                   color_map_norm_type=mpl_colors.Normalize,
                   color_style: str = 'cool',
                   scores_norm_min: int or float = None,
                   scores_norm_max: int or float = None,
                   plot_high_scores_on_top: bool = True,
                   use_variable_line_width: bool = False):
        # Create color map
        if scores_norm_min is None:
            scores_norm_min = min(scores)
        if scores_norm_max is None:
            scores_norm_max = max(scores)
        color_map_norm = color_map_norm_type(vmin=scores_norm_min, vmax=scores_norm_max, clip=True)
        color_map = mpl_cm.ScalarMappable(cmap=color_style, norm=color_map_norm)
        # Sort data sets, scores
        sorted_scores_indices = list(np.argsort(scores))
        if not plot_high_scores_on_top:
            sorted_scores_indices = sorted_scores_indices[::-1]
        self._data_sets = [self._data_sets[idx] for idx in sorted_scores_indices]
        scores = [scores[idx] for idx in sorted_scores_indices]
        # Set plot colors
        colors = list()
        for score in scores:
            colors.append(color_map.to_rgba(score))
        # Store values
        self._color_style = color_style
        self._colors = colors
        self._scores = scores
        self._scores_norm_min = scores_norm_min
        self._scores_norm_max = scores_norm_max
        self._color_map_norm = color_map_norm
        self._use_variable_line_width = use_variable_line_width
        self._sorted_scores_indices = sorted_scores_indices

    def set_labels(self,
                   labels: list):
        for idx in range(self._num_dims - 1):
            self._axes[idx].set_xticklabels([labels[idx]], color='k')
        self._axes[self._num_dims - 1].set_xticklabels([labels[self._num_dims - 2], labels[self._num_dims - 1]],
                                                       color='k')

    def add_color_bar(self,
                      label: str or None = None):
        if self._scores is not None:
            if self._ax_color_bar is None:  # if axis for color bar not yet created
                # create axis for color bar
                self._ax_color_bar = mpl_colorbar.make_axes(self._axes, Location='right')[0]
            self._ax_color_bar.set_ylim([self._scores_norm_min, self._scores_norm_max])
            cb = mpl_colorbar.ColorbarBase(self._ax_color_bar,
                                           cmap=self._color_style,
                                           norm=self._color_map_norm,
                                           orientation='vertical',
                                           )
            if label:
                cb.set_label(label)
        else:
            raise ValueError('Set scores before adding a color bar.')


# Make colors and coordinate values for example
# red


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-p', '--prefixes', nargs='+', type=str, required=False, help="output base")
    parser.add_argument('-o', '--output', type=str, required=False, help="output base")
    args = parser.parse_args()


    assert(len(args.counts) == len(args.prefixes))


    upPref2geneFCP = defaultdict(list)
    downPref2geneFCP = defaultdict(list)

    for fidx, defile in enumerate(args.counts):

        print("Loading file", defile)

        indf = DataFrame.parseFromFile(defile.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
        })

        inHeaders = indf.getHeader()

        assert("ROB_log2FC" in inHeaders)
        assert("ROB_ADJ.PVAL" in inHeaders)

        for row in indf:

            geneSym = row.get("gene_symbol", row.get("id", None))

            if geneSym == None:
                continue

            geneFC = row.get("ROB_log2FC", None)
            genePVal = row.get("ROB_ADJ.PVAL", None)

            if geneFC == None or genePVal == None:
                continue

            if genePVal > 0.05:
                continue

            geneFC = float(geneFC)
            genePVal = float(genePVal)

            if genePVal > 0.05:
                continue

            if geneFC > -0.5:
                upPref2geneFCP[args.prefixes[fidx]].append((geneSym, geneFC, genePVal))

            if geneFC < 0.5:
                downPref2geneFCP[args.prefixes[fidx]].append((geneSym, geneFC, genePVal))


    print("Sorting")
    print([x for x in upPref2geneFCP], [len(upPref2geneFCP[x]) for x in upPref2geneFCP])

    for pref in upPref2geneFCP:
        upPref2geneFCP[pref] = sorted(upPref2geneFCP[pref], key=lambda x: x[2])

    print([x for x in downPref2geneFCP], [len(downPref2geneFCP[x]) for x in downPref2geneFCP])
    for pref in downPref2geneFCP:
        downPref2geneFCP[pref] = sorted(downPref2geneFCP[pref], key=lambda x: x[2])

    print("Intersecting")

    intersectGenesUp = None
    intersectGenesDown = None

    for pref in upPref2geneFCP:
        if intersectGenesUp == None:
            intersectGenesUp = set([x[0] for x in upPref2geneFCP[pref]])
        else:
            intersectGenesUp = intersectGenesUp.intersection([x[0] for x in upPref2geneFCP[pref]])

    for pref in upPref2geneFCP:
        if intersectGenesDown == None:
            intersectGenesDown = set([x[0] for x in downPref2geneFCP[pref]])
        else:
            intersectGenesDown = intersectGenesDown.intersection([x[0] for x in downPref2geneFCP[pref]])

    labels = [args.prefixes[0] + " log2FC", args.prefixes[0] + " rank",
              args.prefixes[1] + " rank", args.prefixes[1] + " log2FC"]

    print(labels, len(intersectGenesUp), len(intersectGenesDown))

    def makeDataDF(indata, ininter):


        newDF = []

        pref0data = {}
        for gsym, gfc, gpv in indata[args.prefixes[0]]:
            pref0data[gsym] = (gfc, gpv)

        pref1data = {}
        for gsym, gfc, gpv in indata[args.prefixes[1]]:
            pref1data[gsym] = (gfc, gpv)

        gene2idx0 = {}
        for gidx, geneT in enumerate(indata[args.prefixes[0]]):
            gene2idx0[geneT[0]] = gidx

        gene2idx1 = {}
        for gidx, geneT in enumerate(indata[args.prefixes[1]]):
            gene2idx1[geneT[0]] = gidx

        for gene in ininter:

            pref0 = pref0data[gene]
            pref1 = pref1data[gene]

            pref0Idx = gene2idx0[gene]
            pref1Idx = gene2idx1[gene]

            dfrow = {
                args.prefixes[0] + " log2FC": pref0[0],
                args.prefixes[0] + " rank": pref0Idx,
                args.prefixes[1] + " log2FC": pref1[0],
                args.prefixes[1] + " rank": pref1Idx,
            }

            newDF.append(dfrow)

        print("in DF", len(newDF))

        return newDF

    print("Making DF")

    upDF = makeDataDF(upPref2geneFCP, intersectGenesUp)
    downDF = makeDataDF(downPref2geneFCP, intersectGenesDown)


    def makePlots(compData, outsuffix):

        upregulatedPlot = outsuffix == "up"

        alldata = []

        for row in compData:
            newdata = []

            for l in labels:

                if l.endswith("log2FC") and outsuffix == "up":
                    newdata.append(-row[l])
                else:
                    newdata.append(row[l])

            alldata.append(tuple(newdata))

        fcMin = min([x[0] for x in alldata] + [x[3] for x in alldata])
        fcMax = max([x[0] for x in alldata] + [x[3] for x in alldata])
        fcRange = fcMax - fcMin

        raMin = min([x[1] for x in alldata] + [x[2] for x in alldata])
        raMax = max([x[1] for x in alldata] + [x[2] for x in alldata])
        raRange = raMax - raMin

        preRanges = {}
        preRanges[0] = {'min': fcMin, 'max': fcMax, 'range': fcRange}
        preRanges[3] = {'min': fcMin, 'max': fcMax, 'range': fcRange}

        preRanges[1] = {'min': raMin, 'max': raMax, 'range': raRange}
        preRanges[2] = {'min': raMin, 'max': raMax, 'range': raRange}


        def globalFormatFunc(col, val):
            print(col, val)

            if col in [0, 3] and upregulatedPlot:
                return '{:4.2f}'.format(-val)

            return '{:4.2f}'.format(val)


        par_co = ParCoord(alldata, figsize=(12, 48))
        par_co.set_colors('#1f77b4')
        par_co.plot(num_ticks=24, preRanges=preRanges, labelFunc=lambda colIdx, val: globalFormatFunc(colIdx, val))
        par_co.set_labels(labels)

        myfig = par_co.fig

        outfilename = args.output + "." + outsuffix + ".png"
        print(outfilename)
        myfig.savefig(outfilename, bbox_inches="tight")

    if len(upDF) > 0:
        makePlots(upDF, "up")
    else:
        fig = plt.figure()
        plt.title("No Data to show")
        outfilename = args.output + ".up.png"
        plt.savefig(outfilename, bbox_inches="tight")
        plt.close()



    if len(downDF) > 0:
        makePlots(downDF, "down")
    else:
        fig = plt.figure()
        plt.title("No Data to show")
        outfilename = args.output + ".down.png"
        plt.savefig(outfilename, bbox_inches="tight")
        plt.close()