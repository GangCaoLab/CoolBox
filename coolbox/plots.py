import sys
import math
import collections

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mpl_toolkits.axisartist as axisartist
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from matplotlib.patches import Rectangle

import numpy as np
from scipy import ndimage


from .utilities import (
    GenomeRange, cm2inch, file_to_intervaltree,
    change_chrom_names, rgb2hex,
    Interval, IntervalTree,
    opener, ReadBed, to_string
)


import logging

FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s"
logging.basicConfig(format=FORMAT)

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = [
    "PlotFrame", "PlotSpacer", "PlotXAxis",
    "PlotBigWig", "PlotABCompartment", "PlotBedGraph", "PlotBed",
    "PlotBoundaries", "PlotTADs", "PlotArcs",
    "PlotCool",
    "PlotVlines", "PlotHighLightRegions",
    "PlotHiCPeaks"
]


class PlotFrame(object):

    DEFAULT_WIDTH = 40
    DEFAULT_WIDTH_RATIOS = (0.93, 0.07)
    DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0, 'top': 1}

    def __init__(self, properties_dict, tracks, *args, **kwargs):

        self.properties = properties_dict
        self.tracks = tracks

        if 'width' not in self.properties:
            self.properties['width'] = PlotFrame.DEFAULT_WIDTH

        if 'width_ratios' not in self.properties:
            self.properties['width_ratios'] = PlotFrame.DEFAULT_WIDTH_RATIOS

        if 'margins' not in self.properties:
            self.properties['margins'] = PlotFrame.DEFAULT_MARGINS

        super().__init__()

    def get_tracks_height(self, default_height=3):
        """
        Get heights of all tracks.

        Return:
            heights (:obj:`list` of `float`): heights of all tracks.
        """
        heights = []
        for track in self.tracks.values():
            if 'height' in track.properties:

                # auto specify height for Cool Track
                if track.properties['height'] == 'cool_auto':
                    cool_height = track.get_tracks_height(
                        self.properties['width'] * self.properties['width_ratios'][0])
                    heights.append(cool_height)
                else:
                    heights.append(track.properties['height'])

            else:
                heights.append(default_height)
        return heights

    def plot(self, chrom, start, end):
        """
        Plot all tracks.
        """

        tracks_height = self.get_tracks_height()
        self.properties['height'] = sum(tracks_height)

        fig = plt.figure(figsize=cm2inch(self.properties['width'],
                                         self.properties['height']))
        if 'title' in self.properties:
            fig.suptitle(self.properties['title'])

        grids = matplotlib.gridspec.GridSpec(
            len(tracks_height), 2,
            height_ratios=tracks_height,
            width_ratios=self.properties['width_ratios'])
        axis_list = []
        for idx, track in enumerate(self.tracks.values()):
            axis = axisartist.Subplot(fig, grids[idx, 0])
            fig.add_subplot(axis)
            axis.axis[:].set_visible(False)
            axis.patch.set_visible(False)
            label_axis = plt.subplot(grids[idx, 1])
            label_axis.set_axis_off()
            try:
                track.plot(axis, label_axis, chrom, start, end)

            except Exception as e:
                log.error("Error occured when plot track 'name: {}, type:{}', {}:{}".format(
                    track.name, type(track),
                    type(e), str(e)))
            # plot coverages
            if hasattr(track, 'coverages'):
                for cov_idx, cov in enumerate(track.coverages):
                    cov.track = track
                    try:
                        cov.plot(axis, chrom, start, end)
                    except Exception as e:
                        log.error("Error occured when plot track's coverage "
                                  "'track name: {}, track type:{}, cov idx: {}, cov type: {}', "
                                  "{}:{}".format(
                            track.name, type(track), cov_idx, type(cov),
                            type(e), str(e)))

            axis_list.append(axis)

        margins = self.properties['margins']
        fig.subplots_adjust(wspace=0, hspace=0.0,
                            left=margins['left'],
                            right=margins['right'],
                            bottom=margins['bottom'],
                            top=margins['top'])

        plt.close()

        return fig


class TrackPlot(object):
    """
    The TrackPlot object is a holder for all tracks that are to be plotted.
    For example, to plot a bedgraph file a new class that extends TrackPlot
    should be created.

    It is expected that all TrackPlot objects have a plot method.

    """

    def __init__(self, *args, **kwargs):
        if not hasattr(self, 'properties'):
            self.properties = args[0]
        super().__init__()


class PlotSpacer(TrackPlot):

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        ax.set_xlim(start_region, end_region)
        pass


class PlotBedGraph(TrackPlot):

    DEFAULT_COLOR = '#a6cee3'

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)
        if 'color' not in self.properties:
            self.properties['color'] = PlotBedGraph.DEFAULT_COLOR
        self.interval_tree, ymin, ymax = file_to_intervaltree(self.properties['file'])

        if 'max_value' not in self.properties or self.properties['max_value'] == 'auto':
            self.properties['max_value'] = ymax

        if 'min_value' not in self.properties or self.properties['min_value'] == 'auto':
            self.properties['min_value'] = ymin

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.ax = ax
        self.label_ax = label_ax
        score_list = []
        pos_list = []

        if chrom_region not in list(self.interval_tree):
            chrom_region = change_chrom_names(chrom_region)

        for region in sorted(self.interval_tree[chrom_region][start_region - 10000:end_region + 10000]):
            score_list.append(float(region.data[0]))
            pos_list.append(region.begin + (region.end - region.begin) / 2)

        if 'color' not in self.properties:
            self.properties['color'] = PlotBedGraph.DEFAULT_COLOR

        if 'extra' in self.properties and self.properties['extra'][0] == '4C':
            # draw a vertical line for each fragment region center
            self.ax.fill_between(pos_list, score_list,
                                 facecolor=self.properties['color'],
                                 edgecolor='none')
            self.ax.vlines(pos_list, [0], score_list, color='olive', linewidth=0.5)
            self.ax.plot(pos_list, score_list, '-', color='slateblue', linewidth=0.7)
        else:
            try:
                self.ax.fill_between(pos_list, score_list, facecolor=self.properties['color'])
            except ValueError:
                log.warning("Invalid color {} for {}. "
                            "Using gray instead.".format(self.properties['color'], self.properties['file']))
                self.ax.fill_between(pos_list, score_list, facecolor='gray')

        self.ax.set_frame_on(False)
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)
        self.ax.set_xlim(start_region, end_region)

        ymax = self.properties['max_value']
        ymin = self.properties['min_value']

        if float(ymax) % 1 == 0:
            ymax_print = int(ymax)
        else:
            ymax_print = "{:.1f}".format(ymax)
        self.ax.set_ylim(ymin, ymax)
        ydelta = ymax - ymin
        small_x = 0.01 * (end_region - start_region)

        if 'show_data_range' in self.properties and \
                self.properties['show_data_range'] == 'no':
            pass
        else:
            # by default show the data range
            self.ax.text(start_region - small_x, ymax - ydelta * 0.2,
                         "[{}-{}]".format(ymin, ymax_print),
                         horizontalalignment='left', size='small',
                         verticalalignment='bottom')

        self.label_ax.text(0.15, 0.5, self.properties['title'],
                           horizontalalignment='left', size='large',
                           verticalalignment='center', transform=self.label_ax.transAxes)


class PlotBigWig(TrackPlot):

    DEFAULT_COLOR = "#33a02c"

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)
        import pyBigWig
        self.bw = pyBigWig.open(self.properties['file'])
        if 'color' not in self.properties:
            self.properties['color'] = PlotBigWig.DEFAULT_COLOR

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.ax = ax
        self.label_ax = label_ax

        genome_range = GenomeRange(chrom_region, start_region, end_region)
        log.debug("plotting {}".format(self.properties['file']))

        num_bins = self.__get_bins_num()
        self.__check_chrom_name(genome_range)
        scores_per_bin = self.__get_scores_per_bin(genome_range, num_bins)

        x_values = np.linspace(genome_range.start, genome_range.end, num_bins)

        if 'type' in self.properties and self.properties != 'fill':
            self.__plot_line_or_points(scores_per_bin, x_values)
        else:
            self.__plot_fill(scores_per_bin, x_values)

        ymin, ymax = self.__adjust_plot(genome_range)

        if "show_data_range" in self.properties and self.properties["show_data_range"] == 'no':
            pass
        else:
            self.__plot_data_range(ymin, ymax, genome_range)

        self.label_ax.text(0.15, 0.5, self.properties['title'],
                           horizontalalignment='left', size='large',
                           verticalalignment='center')

        return self.ax

    def __get_bins_num(self, default_num=700):
        num_bins = default_num
        if 'number_of_bins' in self.properties:
            try:
                num_bins = int(self.properties['number_of_bins'])
            except TypeError:
                num_bins = default_num
                log.warning("'number_of_bins' value: {} for bigwig file {} "
                            "is not valid. Using default value (700)".format(self.properties['number_of_bins'],
                                                                             self.properties['file']))
        return num_bins

    def __get_scores_per_bin(self, genome_range, num_bins, max_try_nums=5):
        # on rare occasions pyBigWig may throw an error, apparently caused by a corruption
        # of the memory. This only occurs when calling trackPlot from different
        # processors. Reloading the file solves the problem.
        num_tries = 0
        while num_tries < max_try_nums:
            num_tries += 1
            try:
                scores_per_bin = np.array(self.bw.stats(genome_range.chrom, genome_range.start,
                                                        genome_range.end, nBins=num_bins)).astype(float)
            except Exception as e:
                import pyBigWig
                self.bw = pyBigWig.open(self.properties['file'])

                log.warning("error found while reading bigwig scores ({}).\nTrying again. Iter num: {}".
                            format(e, num_tries))
                pass
            else:
                if num_tries > 1:
                    log.warning("After {} the scores could be computed".format(num_tries))
                break
        return scores_per_bin

    def __check_chrom_name(self, genome_range):
        if genome_range.chrom not in self.bw.chroms().keys():
            genome_range.change_chrom_names()

        if genome_range.chrom not in self.bw.chroms().keys():
            log.warning("Can not read region {} from bigwig file:\n\n"
                        "{}\n\nPlease check that the chromosome name is part of the bigwig file "
                        "and that the region is valid".format(str(genome_range), self.properties['file']))

    def __plot_line_or_points(self, scores_per_bin, x_values):
        if self.properties['type'].find(":") > 0:
            plot_type, size = self.properties['type'].split(":")
            try:
                size = float(size)
            except ValueError:
                raise ValueError("Invalid value: 'type = {}' in Track: {}\n"
                                 "A number was expected and found '{}'".format(
                                                                   self.properties['type'],
                                                                   self.properties['name'],
                                                                   size))
        else:
            plot_type = self.properties['type']
            size = None

        if plot_type == 'line':
            self.ax.plot(x_values, scores_per_bin, '-', linewidth=size, color=self.properties['color'])

        elif plot_type == 'points':
            self.ax.plot(x_values, scores_per_bin, '.', markersize=size, color=self.properties['color'])

        else:
            raise ValueError("Invalid: 'type = {}' in Track: {}\n".format(self.properties['type'],
                                                                          self.properties['name']))

    def __plot_fill(self, scores_per_bin, x_values):
        if 'positive_color' not in self.properties or 'negative_color' not in self.properties:
            self.ax.fill_between(x_values, scores_per_bin, linewidth=0.1,
                                 color=self.properties['color'],
                                 facecolor=self.properties['color'])
        else:
            self.ax.fill_between(x_values, 0, scores_per_bin, where=(scores_per_bin > 0),
                                 linewidth=0.1, color=self.properties['positive_color'],
                                 facecolor=self.properties['positive_color'])
            self.ax.fill_between(x_values, 0, scores_per_bin, where=(scores_per_bin < 0),
                                 linewidth=0.1, color=self.properties['negative_color'],
                                 facecolor=self.properties['negative_color'])

    def __adjust_plot(self, genome_range):
        self.ax.set_xlim(genome_range.start, genome_range.end)
        ymin, ymax = self.ax.get_ylim()
        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
            ymax = self.properties['max_value']
        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            ymin = self.properties['min_value']

        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            self.ax.set_ylim(ymax, ymin)
        else:
            self.ax.set_ylim(ymin, ymax)
        return ymin, ymax

    def __plot_data_range(self, ymin, ymax, genome_range):
        ydelta = ymax - ymin

        # set min max
        format_lim = lambda lim: int(lim) if float(lim) %1 == 0 else "{:.2f}".format(lim)
        ymax_print = format_lim(ymax)
        ymin_print = format_lim(ymin)
        small_x = 0.01 * genome_range.length
        # by default show the data range
        self.ax.text(genome_range.start - small_x, ymax - ydelta * 0.2,
                     "[ {} ~ {} ]".format(ymin_print, ymax_print),
                     horizontalalignment='left',
                     verticalalignment='top')


class PlotABCompartment(PlotBigWig):

    DEFAULT_POSITIVE_COLOR = "#ff9c9c"
    DEFAULT_NEGATIVE_COLOR = "#66ccff"

    def __init__(self, *args, **kwargs):
        PlotBigWig.__init__(self, *args, **kwargs)
        if 'positive_color' not in self.properties:
            self.properties['positive_color'] = PlotABCompartment.DEFAULT_POSITIVE_COLOR
        if 'negative_color' not in self.properties:
            self.properties['negative_color'] = PlotABCompartment.DEFAULT_NEGATIVE_COLOR


class PlotXAxis(TrackPlot):

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)
        if 'fontsize' not in self.properties:
            self.properties['fontsize'] = 15

    def plot(self, ax, label_axis, chrom_region, region_start, region_end):
        ax.set_xlim(region_start, region_end)
        ticks = ax.get_xticks()
        if ticks[-1] - ticks[1] <= 1e5:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"

        elif 1e5 < ticks[-1] - ticks[1] < 4e6:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"
        else:
            labels = ["{:,.1f} ".format((x / 1e6))
                      for x in ticks]
            labels[-2] += " Mbp"

        ax.axis["x"] = ax.new_floating_axis(0, 0.5)

        ax.axis["x"].axis.set_ticklabels(labels)
        ax.axis['x'].axis.set_tick_params(which='minor', bottom='on')

        ax.axis["x"].major_ticklabels.set(size=int(self.properties['fontsize']))

        if 'where' in self.properties and self.properties['where'] == 'top':
            ax.axis["x"].set_axis_direction("top")


class PlotBoundaries(TrackPlot):

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)

        line_number = 0
        interval_tree = {}
        intervals = []
        prev_chrom = None
        valid_intervals = 0

        with open(self.properties['file'], 'r') as file_h:
            for line in file_h.readlines():
                line_number += 1
                if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
                    continue
                try:
                    chrom, start, end = line.strip().split('\t')[0:3]
                except Exception as detail:
                    msg = 'Could not read line\n{}\n. {}'.format(line, detail)
                    raise IOError(msg)

                try:
                    start = int(start)
                    end = int(end)
                except ValueError as detail:
                    msg = "Error reading line: {}. One of the fields is not " \
                          "an integer.\nError message: {}".format(line_number, detail)
                    raise IOError(msg)

                assert start <= end, "Error in line #{}, end1 larger than start1 in {}".format(line_number, line)

                if prev_chrom and chrom != prev_chrom:
                    start_array, end_array = list(zip(*intervals))
                    start_array = np.array(start_array)
                    end_array = np.array(end_array)
                    # check if intervals are consecutive or 1bp positions demarcating the boundaries
                    if np.any(end_array - start_array == 1):
                        # The file contains only boundaries at 1bp position.
                        end_array = start_array[1:]
                        start_array = start_array[:-1]
                    interval_tree[prev_chrom] = IntervalTree()
                    for idx in range(len(start_array)):
                        interval_tree[prev_chrom].add(Interval(start_array[idx], end_array[idx]))
                        valid_intervals += 1
                    intervals = []

                intervals.append((start, end))

                # each interval spans from the smallest start to the largest end
                prev_chrom = chrom

        start, end = list(zip(*intervals))
        start = np.array(start)
        end = np.array(end)
        # check if intervals are consecutive or 1bp positions demarcating the boundaries
        if np.any(end - start == 1):
            # The file contains only boundaries at 1bp position.
            end = start[1:]
            start = start[:-1]
        interval_tree[chrom] = IntervalTree()
        for idx in range(len(start)):
            interval_tree[chrom].add(Interval(start[idx], end[idx]))
            valid_intervals += 1

        if valid_intervals == 0:
            log.warning("No valid intervals were found in file {}".format(self.properties['file']))

        file_h.close()
        self.interval_tree = interval_tree

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        """
        Plots the boundaries as triangles in the given ax.
        """
        x = []
        y = []
        if chrom_region not in self.interval_tree:
            chrom_region = change_chrom_names(chrom_region)
        for region in sorted(self.interval_tree[chrom_region][start_region:end_region]):
            """
                  /\
                 /  \
                /    \
            _____________________
               x1 x2 x3
            """
            x1 = region.begin
            x2 = x1 + float(region.end - region.begin) / 2
            x3 = region.end
            y1 = 0
            y2 = (region.end - region.begin)
            x.extend([x1, x2, x3])
            y.extend([y1, y2, y1])

        ax.plot(x, y, color='black')
        ax.set_xlim(start_region, end_region)


class PlotBed(TrackPlot):

    DEFAULT_COLOR = "#1f78b4"

    def __init__(self, *args, **kwarg):
        TrackPlot.__init__(self, *args, **kwarg)
        self.bed_type = None  # once the bed file is read, this is bed3, bed6 or bed12
        self.len_w = None  # this is the length of the letter 'w' given the font size
        self.interval_tree = {}  # interval tree of the bed regions

        from matplotlib import font_manager
        if 'fontsize' not in self.properties:
            self.properties['fontsize'] = 12
        else:
            self.properties['fontsize'] = float(self.properties['fontsize'])

        self.fp = font_manager.FontProperties(size=self.properties['fontsize'])

        if 'color' not in self.properties:
            self.properties['color'] = 'bed_rgb'
        if 'border_color' not in self.properties:
            self.properties['border_color'] = 'black'
        if 'labels' not in self.properties:
            self.properties['labels'] = 'auto'
        if 'style' not in self.properties:
            self.properties['style'] = 'flybase'
        if 'display' not in self.properties:
            self.properties['display'] = 'stacked'
        if 'interval_height' not in self.properties:
            self.properties['interval_height'] = 100

        if self.properties['labels'] != 'on':
            self.is_draw_labels = False
        else:
            self.is_draw_labels = True

        self.colormap = None
        # check if the color given is a color map
        if not matplotlib.colors.is_color_like(self.properties['color']) and self.properties['color'] != 'bed_rgb':
            # check if the color is a valid colormap name
            if self.properties['color'] not in matplotlib.cm.datad:
                log.warning("*WARNING* color: '{}' for Track {} is not valid. Color has "
                            "been set to {}".format(self.properties['color'], self.properties['name'],
                                                    PlotBed.DEFAULT_COLOR))
                self.properties['color'] = PlotBed.DEFAULT_COLOR
            else:
                self.colormap = self.properties['color']

        # to set the distance between rows
        self.row_scale = self.properties['interval_height'] * 2.3

        self.interval_tree, min_score, max_score = self.__process_bed()
        if self.colormap is not None:
            if 'min_value' in self.properties:
                min_score = self.properties['min_value']
            if 'max_value' in self.properties:
                max_score = self.properties['max_value']

            norm = matplotlib.colors.Normalize(vmin=min_score,
                                               vmax=max_score)

            cmap = matplotlib.cm.get_cmap(self.properties['color'])
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    def __get_length_w(self, fig_width, region_start, region_end):
        '''
        to improve the visualization of the genes it is good to have an estimation of the label
        length. In the following code I try to get the length of a 'W' in base pairs.
        '''
        if self.is_draw_labels:
            # from http://scipy-cookbook.readthedocs.org/items/Matplotlib_LaTeX_Examples.html
            inches_per_pt = 1.0 / 72.27
            font_in_inches = self.properties['fontsize'] * inches_per_pt
            region_len = region_end - region_start
            bp_per_inch = region_len / fig_width
            font_in_bp = font_in_inches * bp_per_inch
            self.len_w = font_in_bp
            log.debug("len of w set to: {} bp".format(self.len_w))
        else:
            self.len_w = 1

        return self.len_w

    def __process_bed(self):

        bed_file_h = ReadBed(opener(self.properties['file']))
        self.bed_type = bed_file_h.file_type

        if 'color' in self.properties and self.properties['color'] == 'bed_rgb' and \
           self.bed_type not in ['bed12', 'bed9']:
            log.warning("*WARNING* Color set to 'bed_rgb', but bed file does not have the rgb field. The color has "
                        "been set to {}".format(PlotBed.DEFAULT_COLOR))
            self.properties['color'] = PlotBed.DEFAULT_COLOR

        valid_intervals = 0
        interval_tree = {}

        max_score = float('-inf')
        min_score = float('inf')
        for bed in bed_file_h:
            if bed.score < min_score:
                min_score = bed.score
            if bed.score > max_score:
                max_score = bed.score

            if bed.chromosome not in interval_tree:
                interval_tree[bed.chromosome] = IntervalTree()

            interval_tree[bed.chromosome].add(Interval(bed.start, bed.end, bed))
            valid_intervals += 1

        if valid_intervals == 0:
            log.warning("No valid intervals were found in file {}".format(self.properties['file_name']))

        return interval_tree, min_score, max_score

    def __get_max_num_row(self, len_w, small_relative):
        ''' Process the whole bed regions at the given figure length and font size to
        determine the maximum number of rows required.
        :return:
        '''

        self.max_num_row = {}
        for chrom in self.interval_tree:
            row_last_position = []  # each entry in this list contains the end position
            self.max_num_row[chrom] = 0
            for region in sorted(self.interval_tree[chrom][0:500000000]):
                bed = region.data
                if self.is_draw_labels:
                    bed_extended_end = int(bed.end + (len(bed.name) * len_w))
                else:
                    bed_extended_end = (bed.end + 2 * small_relative)

                # get smallest free row
                if len(row_last_position) == 0:
                    free_row = 0
                    row_last_position.append(bed_extended_end)
                else:
                    # get list of rows that are less than bed.start, then take the min
                    idx_list = [idx for idx, value in enumerate(row_last_position) if value < bed.start]
                    if len(idx_list):
                        free_row = min(idx_list)
                        row_last_position[free_row] = bed_extended_end
                    else:
                        free_row = len(row_last_position)
                        row_last_position.append(bed_extended_end)

                if free_row > self.max_num_row[bed.chromosome]:
                    self.max_num_row[bed.chromosome] = free_row

        log.debug("max number of rows set to {}".format(self.max_num_row))
        return self.max_num_row

    def __get_y_pos(self, free_row):
        """
        The y_pos is set such that regions to be plotted do not overlap (stacked). To override this
        the properties['collapsed'] needs to be set.

        The algorithm uses a interval tree (self.region_interval) to check the overlaps
        and a sort of coverage vector 'rows used' to identify the row in which to plot
        :return: int y position
        """

        # if the domain directive is given, ypos simply oscilates between 0 and 100
        if self.properties['display'] == 'interlaced':
            ypos = self.properties['interval_height'] if self.counter % 2 == 0 else 1

        elif self.properties['display'] == 'collapsed':
            ypos = 0

        else:
            ypos = free_row * self.row_scale
        return ypos

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.counter = 0
        self.small_relative = 0.004 * (end_region - start_region)
        self.__get_length_w(ax.get_figure().get_figwidth(), start_region, end_region)
        if 'global_max_row' in self.properties and self.properties['global_max_row'] == 'yes':
            self.__get_max_num_row(self.len_w, self.small_relative)

        if chrom_region not in self.interval_tree.keys():
            chrom_region = change_chrom_names(chrom_region)

        genes_overlap = sorted(self.interval_tree[chrom_region][start_region:end_region])

        if self.properties['labels'] == 'auto':
            if len(genes_overlap) > 60:
                # turn labels off when too many intervals are visible.
                self.is_draw_labels = False
            else:
                self.is_draw_labels = True

        max_num_row_local = 1
        max_ypos = 0
        # check for the number of other intervals that overlap
        #    with the given interval
        #            1         2
        #  012345678901234567890123456
        #  1=========       4=========
        #       2=========
        #         3============
        #
        # for 1 row_last_position = [9]
        # for 2 row_last_position = [9, 14]
        # for 3 row_last_position = [9, 14, 19]
        # for 4 row_last_position = [26, 14, 19]

        row_last_position = []  # each entry in this list contains the end position
        # of genomic interval. The list index is the row
        # in which the genomic interval was plotted.
        # Any new genomic interval that wants to be plotted,
        # knows the row to use by finding the list index that
        # is larger than its start

        # check for overlapping genes including
        # label size (if plotted)

        for region in genes_overlap:
            """
            BED12 gene format with exon locations at the end
            chrX    20850   23076   CG17636-RA      0       -       20850   23017   0       3       946,765,64,     0,1031,2162,

            BED9
            bed with rgb at end
            chr2L   0       70000   ID_5    0.26864549832   .       0       70000   51,160,44

            BED6
            bed without rgb
            chr2L   0       70000   ID_5    0.26864549832   .
            """
            self.counter += 1
            bed = region.data

            if self.is_draw_labels:
                num_name_characters = len(bed.name) + 2  # +2 to account for an space before and after the name
                bed_extended_end = int(bed.end + (num_name_characters * self.len_w))
            else:
                bed_extended_end = (bed.end + 2 * self.small_relative)

            # get smallest free row
            if len(row_last_position) == 0:
                free_row = 0
                row_last_position.append(bed_extended_end)
            else:
                # get list of rows that are less than bed.start, then take the min
                idx_list = [idx for idx, value in enumerate(row_last_position) if value < bed.start]
                if len(idx_list):
                    free_row = min(idx_list)
                    row_last_position[free_row] = bed_extended_end
                else:
                    free_row = len(row_last_position)
                    row_last_position.append(bed_extended_end)

            rgb, edgecolor = self.__get_rgb_and_edge_color(bed)

            ypos = self.__get_y_pos(free_row)

            # do not plot if the maximum interval rows to plot is reached
            if 'gene_rows' in self.properties and free_row >= int(self.properties['gene_rows']):
                continue

            if free_row > max_num_row_local:
                max_num_row_local = free_row
            if ypos > max_ypos:
                max_ypos = ypos

            if self.bed_type == 'bed12':
                if self.properties['style'] == 'flybase':
                    self.__draw_gene_with_introns_flybase_style(ax, bed, ypos, rgb, edgecolor)
                else:
                    self.__draw_gene_with_introns(ax, bed, ypos, rgb, edgecolor)
            else:
                self.__draw_gene_simple(ax, bed, ypos, rgb, edgecolor)

            if not self.is_draw_labels:
                pass
            elif bed.start > start_region and bed.end < end_region:
                ax.text(bed.end + self.small_relative, ypos + (float(self.properties['interval_height']) / 2),
                        bed.name, horizontalalignment='left',
                        verticalalignment='center', fontproperties=self.fp)

        if self.counter == 0:
            log.warning("*Warning* No intervals were found for file {} "
                        "in Track '{}' for the interval plotted ({}:{}-{}).\n".
                        format(self.properties['file'], self.properties['name'], chrom_region, start_region, end_region))

        ymax = 0

        if 'global_max_row' in self.properties and self.properties['global_max_row'] == 'yes':
            ymin = self.max_num_row[chrom_region] * self.row_scale

        elif 'gene_rows' in self.properties:
            ymin = int(self.properties['gene_rows']) * self.row_scale
        else:
            ymin = max_ypos + self.properties['interval_height']

        log.debug("ylim {},{}".format(ymin, ymax))
        # the axis is inverted (thus, ymax < ymin)
        ax.set_ylim(ymin, ymax)

        if 'display' in self.properties:
            if self.properties['display'] == 'domain':
                ax.set_ylim(-5, 205)
            elif self.properties['display'] == 'collapsed':
                ax.set_ylim(-5, 105)

        ax.set_xlim(start_region, end_region)

        label_ax.text(0.15, 1, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='top', transform=label_ax.transAxes)

    def __get_rgb_and_edge_color(self, bed):
        rgb = self.properties['color']
        edgecolor = self.properties['border_color']

        if self.colormap:
            # translate value field (in the example above is 0 or 0.2686...) into a color
            rgb = self.colormap.to_rgba(bed.score)

        if self.properties['color'] == 'bed_rgb':
            # if rgb is set in the bed line, this overrides the previously
            # defined colormap
            if self.bed_type in ['bed9', 'bed12'] and len(bed.rgb) == 3:
                try:
                    rgb = [float(x) / 255 for x in bed.rgb]
                    if 'border_color' in self.properties:
                        edgecolor = self.properties['border_color']
                    else:
                        edgecolor = self.properties['color']
                except IndexError:
                    rgb = PlotBed.DEFAULT_COLOR
            else:
                rgb = PlotBed.DEFAULT_COLOR
        return rgb, edgecolor

    def __draw_gene_simple(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws an interval with direction (if given)
        """
        from matplotlib.patches import Polygon

        if bed.strand not in ['+', '-']:
            ax.add_patch(Rectangle((bed.start, ypos), bed.end - bed.start, self.properties['interval_height'],
                                   edgecolor=edgecolor, facecolor=rgb, linewidth=0.5))
        else:
            vertices = self.__draw_arrow(ax, bed.start, bed.end, bed.strand, ypos)
            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=rgb,
                                 linewidth=0.5))

    def __draw_gene_with_introns_flybase_style(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws a gene using different styles
        """
        from matplotlib.patches import Polygon
        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self.__draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
            return
        half_height = float(self.properties['interval_height']) / 2
        # draw 'backbone', a line from the start until the end of the gene
        ax.plot([bed.start, bed.end], [ypos + half_height, ypos + half_height], 'black', linewidth=0.5, zorder=-1)

        # get start, end of all the blocks
        positions = []
        for idx in range(0, bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if x0 < bed.thick_start < x1:
                positions.append((x0, bed.thick_start, 'UTR'))
                positions.append((bed.thick_start, x1, 'coding'))

            elif x0 < bed.thick_end < x1:
                positions.append((x0, bed.thick_end, 'coding'))
                positions.append((bed.thick_end, x1, 'UTR'))

            else:
                if x1 < bed.thick_start or x0 > bed.thick_end:
                    type = 'UTR'
                else:
                    type = 'coding'

                positions.append((x0, x1, type))

        # plot all blocks as rectangles except the last if the strand is + or
        # the first is the strand is -, which are drawn as arrows.
        if bed.strand == '-':
            positions = positions[::-1]

        first_pos = positions.pop()
        if first_pos[2] == 'UTR':
            _rgb = 'grey'
        else:
            _rgb = rgb

        vertices = self.__draw_arrow(ax, first_pos[0], first_pos[1], bed.strand, ypos)

        ax.add_patch(Polygon(vertices, closed=True, fill=True,
                             edgecolor=edgecolor,
                             facecolor=_rgb,
                             linewidth=0.5))

        for start_pos, end_pos, _type in positions:
            if _type == 'UTR':
                _rgb = 'grey'
            else:
                _rgb = rgb
            vertices = [(start_pos, ypos), (start_pos, ypos + self.properties['interval_height']),
                        (end_pos, ypos + self.properties['interval_height']), (end_pos, ypos)]

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=_rgb,
                                 linewidth=0.5))

    def __draw_arrow(self, ax, start, end, strand, ypos):
        """
        Draws a filled arrow
        :param ax:
        :param start:
        :param end:
        :param strand:
        :param ypos:
        :param rgb:
        :return: None
        """
        half_height = float(self.properties['interval_height']) / 2
        if strand == '+':
            x0 = start
            x1 = end  # - self.small_relative
            y0 = ypos
            y1 = ypos + self.properties['interval_height']

            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

            -----------------\
            ---------------- /

            """

            vertices = [(x0, y0), (x0, y1), (x1, y1), (x1 + self.small_relative, y0 + half_height), (x1, y0)]

        else:
            x0 = start  # + self.small_relative
            x1 = end
            y0 = ypos
            y1 = ypos + self.properties['interval_height']

            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

            /-----------------
            \-----------------
            """

            vertices = [(x0, y0), (x0 - self.small_relative, y0 + half_height), (x0, y1), (x1, y1), (x1, y0)]

        return vertices

    def __draw_gene_with_introns(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws a gene like in flybase gbrowse.
        """
        from matplotlib.patches import Polygon

        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self.__draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
            return
        half_height = float(self.properties['interval_height']) / 2
        quarter_height = float(self.properties['interval_height']) / 4
        three_quarter_height = quarter_height * 3

        # draw 'backbone', a line from the start until the end of the gene
        ax.plot([bed.start, bed.end], [ypos + half_height, ypos + half_height], 'black', linewidth=0.5, zorder=-1)

        for idx in range(0, bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if x1 < bed.thick_start or x0 > bed.thick_end:
                y0 = ypos + quarter_height
                y1 = ypos + three_quarter_height
            else:
                y0 = ypos
                y1 = ypos + self.properties['interval_height']

            if x0 < bed.thick_start < x1:
                vertices = ([(x0, ypos + quarter_height), (x0, ypos + three_quarter_height),
                             (bed.thick_start, ypos + three_quarter_height),
                             (bed.thick_start, ypos + self.properties['interval_height']),
                             (bed.thick_start, ypos + self.properties['interval_height']),
                             (x1, ypos + self.properties['interval_height']), (x1, ypos),
                             (bed.thick_start, ypos), (bed.thick_start, ypos + quarter_height)])

            elif x0 < bed.thick_end < x1:
                vertices = ([(x0, ypos),
                             (x0, ypos + self.properties['interval_height']),
                             (bed.thick_end, ypos + self.properties['interval_height']),
                             (bed.thick_end, ypos + three_quarter_height),
                             (x1, ypos + three_quarter_height),
                             (x1, ypos + quarter_height),
                             (bed.thick_end, ypos + quarter_height),
                             (bed.thick_end, ypos)])
            else:
                vertices = ([(x0, y0), (x0, y1), (x1, y1), (x1, y0)])

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 linewidth=0.1,
                                 edgecolor='none',
                                 facecolor=rgb))

            if idx < bed.block_count - 1:
                # plot small arrows using the character '<' or '>' over the back bone
                intron_length = bed.block_starts[idx + 1] - (bed.block_starts[idx] + bed.block_sizes[idx])
                marker = 5 if bed.strand == '+' else 4
                if intron_length > 3 * self.small_relative:
                    pos = np.arange(x1 + 1 * self.small_relative,
                                    x1 + intron_length + self.small_relative, int(2 * self.small_relative))
                    ax.plot(pos, np.zeros(len(pos)) + ypos + half_height, '.', marker=marker,
                            fillstyle='none', color='blue', markersize=3)

                elif intron_length > self.small_relative:
                    intron_center = x1 + int(intron_length) / 2
                    ax.plot([intron_center], [ypos + half_height], '.', marker=5,
                            fillstyle='none', color='blue', markersize=3)


class PlotArcs(TrackPlot):

    def __init__(self, *args, **kwarg):
        TrackPlot.__init__(self, *args, **kwarg)
        # the file format expected is similar to file format of links in
        # circos:
        # chr1 100 200 chr1 250 300 0.5
        # where the last value is a score.

        valid_intervals = 0
        interval_tree = {}
        line_number = 0
        with open(self.properties['file'], 'r') as file_h:
            for line in file_h.readlines():
                line_number += 1
                if line.startswith('browser') or line.startswith('track') or line.startswith('#') or self.__is_header(line):
                    continue
                try:
                    chrom1, start1, end1, chrom2, start2, end2, score, *other = line.strip().split('\t')
                except Exception as detail:
                    msg = 'File not valid. The format is chrom1 start1, end1, ' \
                          'chrom2, start2, end2, score\nError: {}\n in line\n {}'.format(detail, line)
                    raise IOError(msg)

                try:
                    start1 = int(start1)
                    end1 = int(end1)
                    start2 = int(start2)
                    end2 = int(end2)
                except ValueError as detail:
                    msg = "Error reading line: {}. One of the fields is not " \
                          "an integer.\nError message: {}".format(line_number, detail)
                    raise IOError(msg)

                assert start1 <= end1, "Error in line #{}, end1 larger than start1 in {}".format(line_number, line)
                assert start2 <= end2, "Error in line #{}, end2 larger than start2 in {}".format(line_number, line)
                try:
                    score = float(score)
                except ValueError as detail:
                    msg = "Error reading line: {}. The score is not valid {}. " \
                          "\nError message: {}".format(line_number, score, detail)
                    raise IOError(msg)

                if chrom1 != chrom2:
                    log.warning("Only links in same chromosome are used. Skipping line\n{}\n".format(line))
                    continue

                if chrom1 not in interval_tree:
                    interval_tree[chrom1] = IntervalTree()

                if start2 < start1:
                    start1, start2 = start2, start1
                    end1, end2 = end2, end1

                # each interval spans from the smallest start to the largest end
                interval_tree[chrom1].add(Interval(start1, end2, score))
                valid_intervals += 1

        if valid_intervals == 0:
            log.warning("No valid intervals were found in file {}".format(self.properties['file']))

        file_h.close()
        self.interval_tree = interval_tree

        if 'color' not in self.properties:
            self.properties['color'] = 'blue'

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 0.8

    def plot(self, ax, label_ax, chrom_region, region_start, region_end):
        """
        Makes and arc connecting two points on a linear scale representing
        interactions between Hi-C bins.
        :param ax: matplotlib axis
        :param label_ax: matplotlib axis for labels
        """
        from matplotlib.patches import Arc
        height = 1
        max_diameter = 0
        count = 0
        if chrom_region not in list(self.interval_tree):
            chrom_region = change_chrom_names(chrom_region)
        arcs_in_region = sorted(self.interval_tree[chrom_region][region_start:region_end])

        for idx, interval in enumerate(arcs_in_region):
            # skip arcs whose start and end are outside the plotted region
            if interval.begin < region_start and interval.end > region_end:
                continue

            if 'line_width' in self.properties:
                line_width = float(self.properties['line_width'])
            else:
                line_width = 0.5 * np.sqrt(interval.data)

            diameter = (interval.end - interval.begin)
            center = (interval.begin + interval.end) / 2
            if diameter > max_diameter:
                max_diameter = diameter
            count += 1
            ax.plot([center], [diameter])
            ax.add_patch(Arc((center, 0), diameter,
                             height*2, 0, 0, 180, color=self.properties['color'], lw=line_width))

        # increase max_diameter slightly to avoid cropping of the arcs.
        #max_diameter += max_diameter * 0.05
        height += height * 0.05
        log.debug("{} were arcs plotted".format(count))
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(height, 0.001)
        else:
            ax.set_ylim(-0.001, height)

        ax.set_xlim(region_start, region_end)
        log.debug('title is {}'.format(self.properties['title']))
        label_ax.text(0.15, 0.5, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='center')

    def __is_header(self, line):
        fields = line.split()
        for idx in [1, 2, 4, 5]:
            if not fields[idx].isdigit():
                return True
        return False


class PlotTADs(PlotBed):

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        """
        Plots the boundaries as triangles in the given ax.
        """
        from matplotlib.patches import Polygon
        ymax = 0.001
        valid_regions = 0
        if chrom_region not in self.interval_tree:
            orig = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            log.info('Chromosome name: {} does not exists. Changing name to {}'.format(orig, chrom_region))

        for region in sorted(self.interval_tree[chrom_region][start_region:end_region]):
            """
                  /\
                 /  \
                /    \
            _____________________
               x1 x2 x3
            """
            x1 = region.begin
            x2 = x1 + float(region.end - region.begin) / 2
            x3 = region.end
            y1 = 0
            y2 = (region.end - region.begin)

            rgb, edgecolor = self.__get_rgb_and_edge_color(region.data)

            triangle = Polygon(np.array([[x1, y1], [x2, y2], [x3, y1]]), closed=True,
                               facecolor=rgb, edgecolor=edgecolor)
            ax.add_artist(triangle)
            valid_regions += 1

            if y2 > ymax:
                ymax = y2

        if valid_regions == 0:
            log.warning("No regions found for Track {}.".format(self.properties['name']))

        ax.set_xlim(start_region, end_region)
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(ymax, 0)
        else:
            ax.set_ylim(0, ymax)

        label_ax.text(0.15, 0.5, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='center')


class PlotCool(TrackPlot):

    DEFAULT_COLOR = 'YlOrRd'

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)

        import cooler
        self.cool = cooler.Cooler(self.properties['file'])

        self.__set_default_properties()

        self.small_value = 1e-12
        self.ax = None
        self.label_ax = None
        self.matrix = None

    def __set_default_properties(self):
        self.properties['height'] = 'cool_auto'

        if 'color' not in self.properties:
            self.properties['color'] = PlotCool.DEFAULT_COLOR
        if 'triangular' not in self.properties:
            self.properties['triangular'] = 'yes'
        if 'balance' not in self.properties:
            self.properties['balance'] = 'no'
        if 'color_bar' not in self.properties:
            self.properties['color_bar'] = 'yes'
        if 'transform' not in self.properties:
            self.properties['transform'] = 'no'
        if 'title' not in self.properties:
            self.properties['title'] = ''
        if 'depth_ratio' not in self.properties:
            self.properties['depth_ratio'] = 'full'
        if 'norm' not in self.properties:
            self.properties['norm'] = 'log'

    @property
    def is_inverted(self):
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            return True
        else:
            # default: not inverted
            return False

    @property
    def is_triangular(self):
        if 'triangular' in self.properties and self.properties['triangular'] == 'no':
            return False
        else:
            # default: triangular form
            return True

    @property
    def is_balance(self):
        if 'balance' in self.properties and self.properties['balance'] == 'no':
            return False
        else:
            # default: balance
            return True

    def __transform_matrix(self, arr):
        if self.properties['transform'] == 'log10':
            arr = np.log10(arr)
        elif self.properties['transform'] == 'log2':
            arr = np.log2(arr)
        elif self.properties['transform'] == 'log':
            arr = np.log(arr)
        return arr

    @property
    def matrix_val_range(self):
        small = 1e-4
        arr = self.matrix
        arr_no_nan = arr[np.logical_not(np.isnan(arr))]

        if self.properties['min_value'] == 'auto':
            min_ = small
        else:
            min_ = self.properties['min_value']

        if self.properties['max_value'] == 'auto':
            max_ = arr_no_nan.max()
        else:
            max_ = self.properties['max_value']

        return min_, max_

    def __fetch_matrix(self, genome_range):
        chrom, start, end = genome_range.chrom, genome_range.start, genome_range.end
        if chrom not in self.cool.chromnames:
            chrom = change_chrom_names(chrom)

        range_str = str(GenomeRange(chrom, start, end))
        arr = self.cool.matrix(balance=self.is_balance).fetch(range_str)

        # fill zero and nan with small value
        small = self.small_value
        arr[arr == 0] = small
        arr[np.isnan(arr)] = small
        return arr

    def __get_triangular_matrix(self, arr):
        small = self.small_value
        tri_matrix = ndimage.rotate(arr, 45, prefilter=False, cval=small)

        rows = tri_matrix.shape[0]

        tri_matrix = tri_matrix[0:(rows//2 + 1), :]

        # cut depth
        if self.properties['depth_ratio'] != 'auto' and self.properties['depth_ratio'] != 'full':
            depth_ratio = float(self.properties['depth_ratio'])
            depth = int(tri_matrix.shape[0] * depth_ratio)
            tri_matrix = tri_matrix[-depth:, :]

        return tri_matrix

    def __plot_matrix(self, genome_range):
        start, end = genome_range.start, genome_range.end
        ax = self.ax
        arr = self.matrix
        cmap = plt.get_cmap(self.properties['color'])
        cmap.set_bad("white")
        cmap.set_under("white")
        c_min, c_max = self.matrix_val_range
        if self.is_triangular:
            tri_matrix = self.__get_triangular_matrix(arr)
            img = ax.matshow(tri_matrix, cmap=cmap,
                             extent=(start, end, 0, (end - start)/2),
                             aspect='auto')
        else:
            img = ax.matshow(arr, cmap=cmap,
                             extent=(start, end, end, start),
                             aspect='auto')

        if self.properties['norm'] == 'log':
            img.set_norm(colors.LogNorm(vmin=c_min, vmax=c_max))
        else:
            img.set_norm(colors.Normalize(vmin=c_min, vmax=c_max))

        return img

    def __adjust_figure(self, genome_range):
        ax = self.ax
        start, end = genome_range.start, genome_range.end
        if self.is_triangular:
            if self.is_inverted:
                ax.set_ylim(genome_range.length / 2, 0)
            else:
                ax.set_ylim(0, genome_range.length / 2)
        else:
            ax.set_ylim(end, start)
        ax.set_xlim(start, end)

    def __plot_colorbar(self, img):
        ax_divider = make_axes_locatable(self.ax)
        if self.is_inverted:
            cax = ax_divider.append_axes("top", size=0.09, pad=0.2)
        else:
            cax = ax_divider.append_axes("bottom", size=0.09, pad=0.2)
        colorbar(img, cax=cax, orientation='horizontal')

    def __plot_label(self):
        self.label_ax.text(0.15, 0.5, self.properties['title'],
                           horizontalalignment='left', size='large',
                           verticalalignment='center', transform=self.label_ax.transAxes)

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):

        log.debug("plotting {}".format(self.properties['file']))

        genome_range = GenomeRange(chrom_region, start_region, end_region)

        self.ax = ax
        self.label_ax = label_ax

        # fetch matrix and perform transform process
        arr = self.__fetch_matrix(genome_range)
        if 'transform' in self.properties and self.properties['transform'] != 'no':
            arr = self.__transform_matrix(arr)

        self.matrix = arr

        # plot matrix
        img = self.__plot_matrix(genome_range)
        self.__adjust_figure(genome_range)

        # plot colorbar
        if self.properties['color_bar'] == 'yes':
            self.__plot_colorbar(img)
        else:
            pass

        # plot label
        self.__plot_label()

    def get_tracks_height(self, frame_width):
        """
        calculate track height dynamically.
        """
        if self.properties['triangular'] == 'yes':
            cool_height = frame_width * 0.5
        else:
            cool_height = frame_width

        if 'depth_ratio' in self.properties and self.properties['depth_ratio'] != 'full':
            cool_height = cool_height * self.properties['depth_ratio']

        if 'color_bar' in self.properties and self.properties['color_bar'] != 'no':
            cool_height += 1.5

        return cool_height


class CoveragePlot(object):

    """
    Coverage plot holder, similar to TrackPlot.
    """

    def __init__(self, *args, **kwargs):
        if not hasattr(self, 'properties'):
            self.properties = args[0]
        super().__init__()


class PlotVlines(CoveragePlot):

    DEFAULT_LINE_WIDTH = 0.5
    DEFAULT_LINE_STYLE = 'dashed'
    DEFAULT_COLOR = '#1e1e1e'
    DEFAULT_ALPHA = 0.8

    def __init__(self, *args, **kwargs):

        CoveragePlot.__init__(self, *args, **kwargs)
        if 'line_width' not in self.properties:
            self.properties['line_width'] = PlotVlines.DEFAULT_LINE_WIDTH
        if 'line_style' not in self.properties:
            self.properties['line_style'] = PlotVlines.DEFAULT_LINE_STYLE
        if 'color' not in self.properties:
            self.properties['color'] = PlotVlines.DEFAULT_COLOR
        if 'alpha' not in self.properties:
            self.properties['alpha'] = PlotVlines.DEFAULT_ALPHA
        if 'chr' not in self.properties:
            self.properties['chr'] = None

        if 'file' in self.properties:
            # plot vlines from file
            self.vlines_intval_tree, _, _ = file_to_intervaltree(self.properties['file'])

    def __extract_vlines_from_file(self, chrom, start, end):
        vlines_list = []

        if chrom not in list(self.vlines_intval_tree):
            chrom = change_chrom_names(chrom)

        for region in sorted(self.vlines_intval_tree[chrom][start-1:end+1]):
            vlines_list.append(region.begin)

        return vlines_list

    def __get_vlines_from_properties(self, chrom):
        if self.properties['chr'] is not None:
            chr = self.properties['chr']
            if chr == chrom:
                vlines_list = self.properties['vlines_list']
            else:
                vlines_list = []
        else:
            vlines_list = self.properties['vlines_list']

        return vlines_list

    def plot(self, ax, chrom_region, start_region, end_region):
        if 'file' in self.properties:
            vlines_list = self.__extract_vlines_from_file(chrom_region, start_region, end_region)
        else:
            vlines_list = self.__get_vlines_from_properties(chrom_region)

        ymin, ymax = ax.get_ylim()        

        ax.vlines(vlines_list, ymin, ymax,
                  linestyle=self.properties['line_style'],
                  linewidth=self.properties['line_width'],
                  color=self.properties['color'],
                  alpha=self.properties['alpha'])


class PlotHighLightRegions(CoveragePlot):

    DEFAULT_COLOR = '#ff5d0f'
    DEFAULT_ALPHA = 0.6
    DEFAULT_BORDER_LINE_COLOR = '#000000'
    DEFAULT_BORDER_LINE_ALPHA = 0.8
    DEFAULT_BORDER_LINE_WIDTH = 0.5
    DEFAULT_BORDER_LINE_STYLE = 'dashed'

    def __init__(self, *args, **kwargs):

        CoveragePlot.__init__(self, *args, **kwargs)

        if 'color' not in self.properties:
            self.properties['color'] = PlotHighLightRegions.DEFAULT_COLOR
        if 'alpha' not in self.properties:
            self.properties['alpha'] = PlotHighLightRegions.DEFAULT_BORDER_LINE_ALPHA
        if 'border_line' not in self.properties:
            self.properties['border_line'] = 'yes'
        if 'border_line_color' not in self.properties:
            self.properties['border_line_color'] = PlotHighLightRegions.DEFAULT_BORDER_LINE_COLOR
        if 'border_line_alpha' not in self.properties:
            self.properties['border_line_alpha'] = PlotHighLightRegions.DEFAULT_BORDER_LINE_ALPHA
        if 'border_line_width' not in self.properties:
            self.properties['border_line_width'] = PlotHighLightRegions.DEFAULT_BORDER_LINE_WIDTH
        if 'border_line_style' not in self.properties:
            self.properties['border_line_style'] = PlotHighLightRegions.DEFAULT_BORDER_LINE_STYLE

        if 'file' in self.properties:
            # from bed file
            self.interval_tree = self.__process_bed()
        else:
            # from self.properties['regions']
            self.regions = self.properties['highlight_regions']

    def __extract_regions_from_file(self, chrom, start, end):
        regions = []

        if chrom not in list(self.interval_tree):
            chrom = change_chrom_names(chrom)

        for region in sorted(self.interval_tree[chrom][start-10000:end+10000]):
            regions.append((region.begin, region.end, region.data))

        return regions

    def __get_regions_from_properties(self, chrom):

        def regions_with_color(regions):
            color = self.properties['color']
            return [(start, end, color) for (start, end) in self.regions]

        if self.properties['chr'] is not None:
            if chrom == self.properties['chr']:
                regions = regions_with_color(self.regions)
            else:
                regions = []
        else:
            regions = regions_with_color(self.regions)

        return regions

    def plot(self, ax, chrom_region, start_region, end_region):
        if 'file' in self.properties:
            regions = self.__extract_regions_from_file(chrom_region, start_region, end_region)
        else:
            regions = self.__get_regions_from_properties(chrom_region)

        for (start, end, color) in regions:
            if self.properties['color'] != 'bed_rgb':
                color = self.properties['color']
            if type(color) is not str:
                color = rgb2hex(*color)

            ax.axvspan(start, end, color=color, alpha=self.properties['alpha'])

            if self.properties['border_line'] == 'yes':
                # plot border line
                ymin, ymax = ax.get_ylim()
                ax.vlines([start, end], ymin, ymax,
                          linestyle=self.properties['border_line_style'],
                          linewidth=self.properties['border_line_width'],
                          color=self.properties['border_line_color'],
                          alpha=self.properties['border_line_alpha'])

    def __process_bed(self):
        bed = ReadBed(opener(self.properties['file']))

        if self.properties['color'] == 'bed_rgb' and bed.file_type not in ['bed12', 'bed9']:
            log.warning("*WARNING* Color set to 'bed_rgb', but bed file does not have the rgb field. The color has "
                        "been set to {}".format(PlotBed.DEFAULT_COLOR))
            self.properties['color'] = PlotBed.DEFAULT_COLOR

        interval_tree = {}

        for intval in bed:

            if intval.chromosome not in interval_tree:
                interval_tree[intval.chromosome] = IntervalTree()

            if self.properties['color'] == 'bed_rgb':
                color = intval.rgb
            else:
                color = self.properties['color']
            interval_tree[intval.chromosome].add(Interval(intval.start, intval.end, color))

        return interval_tree


class PlotHiCPeaks(CoveragePlot):

    FILL = "no"
    DEFAULT_FILL_COLOR = "#ff0f39"
    DEFAULT_COLOR = "#ff0f39"
    DEFAULT_ALPHA = 0.8
    DEFAULT_LINEWIDTH = 1.5
    DEFAULT_LINESTYLE = "solid"

    def __init__(self, *args, **kwargs):

        CoveragePlot.__init__(self, *args, **kwargs)

        if 'color' not in self.properties:
            self.properties['color'] = PlotHiCPeaks.DEFAULT_COLOR
        if 'alpha' not in self.properties:
            self.properties['alpha'] = PlotHiCPeaks.DEFAULT_ALPHA
        if 'line_width' not in self.properties:
            self.properties['line_width'] = PlotHiCPeaks.DEFAULT_LINEWIDTH
        if 'line_style' not in self.properties:
            self.properties['line_style'] = PlotHiCPeaks.DEFAULT_LINESTYLE
        if 'fill' not in self.properties:
            self.properties['fill'] = PlotHiCPeaks.FILL
        if 'fill_color' not in self.properties:
            self.properties['fill_color'] = PlotHiCPeaks.DEFAULT_FILL_COLOR

        self.track = None

        self.LoopInverval = collections.namedtuple("LoopInterval",
            ("chr1", "x1", "x2", "chr2", "y1", "y2", "color"))
        self.interval_tree = self.__process_loop_file()

    def plot(self, ax, chrom_region, start_region, end_region):

        if chrom_region not in self.interval_tree:
            chrom_region = change_chrom_names(chrom_region)

        for intval in sorted(self.interval_tree[chrom_region][start_region:end_region]):
            loop = intval.data

            if self.properties['color'] == 'rgb' or 'bed_rgb':
                color = loop.color
            else:
                color = self.properties['color']

            if self.properties['fill_color'] == 'rgb' or 'bed_rgb':
                fill_color = loop.color
            else:
                fill_color = self.properties['fill_color']

            fill = True if self.properties['fill'] == 'yes' else False

            try:
                self.properties['triangular'] = self.track.properties['triangular']
            except KeyError:
                self.properties['triangular'] = 'yes'
                log.warning("*WARNING* 'self.track' attribute not set, "
                            "use default setting(self.properties['triangular'] = 'yes')")

            if self.properties['triangular'] != 'no':

                x, y, (w, h) = self.__get_position_and_size(loop.x1, loop.x2, loop.y1, loop.y2)
                rec = Rectangle((x, y), w, h, angle=45,
                                fill=fill,
                                alpha=self.properties['alpha'],
                                facecolor=fill_color,
                                edgecolor=color,
                                linewidth=self.properties['line_width'],
                                linestyle=self.properties['line_style'])
                ax.add_patch(rec)
            
            else:

                # plot upper rectangle
                x, y, (w, h) = self.__get_position_and_size(loop.x1, loop.x2, loop.y1, loop.y2,
                                                          triangular=False, pos="upper")
                rec = Rectangle((x, y), w, h,
                                fill=fill,
                                alpha=self.properties['alpha'],
                                facecolor=fill_color,
                                edgecolor=color,
                                linewidth=self.properties['line_width'],
                                linestyle=self.properties['line_style'])
                ax.add_patch(rec)

                # plot lower rectangle
                x, y, (w, h) = self.__get_position_and_size(loop.x1, loop.x2, loop.y1, loop.y2,
                                                          triangular=False, pos="lower")
                rec = Rectangle((x, y), w, h,
                                fill=fill,
                                alpha=self.properties['alpha'],
                                facecolor=fill_color,
                                edgecolor=color,
                                linewidth=self.properties['line_width'],
                                linestyle=self.properties['line_style'])
                ax.add_patch(rec)

    def __get_position_and_size(self, start1, end1, start2, end2,
                                triangular=True, pos="upper"):
        """
        Calculate the position and size of the box from loop's start and end postion.
        """
        if triangular:
            m1 = (start1 + end1) / 2
            m2 = (start2 + end2) / 2
            x = (m1 + m2) / 2
            y = x - m1
            w = ( (end2 - start2) / 2 ) / math.cos(math.pi/4)
            h = ( (end1 - start1) / 2 ) / math.cos(math.pi/4)
        else:
            if pos == "upper":
                x = (start2 + end2) / 2
                y = (start1 + end1) / 2
                w = end2 - start2
                h = end1 - start1
            else:
                x = (start1 + end1) / 2
                y = (start2 + end2) / 2
                w = end1 - start1
                h = end2 - start2
        return x, y, (w, h)

    def __process_loop_file(self):
        interval_tree = {}

        with opener(self.properties['file']) as f:
            for idx, line in enumerate(f):
                line = to_string(line)
                # skip header line
                if idx == 0 and self.__is_header(line):
                    continue

                fields = line.split()
                chr1, x1, x2, chr2, y1, y2, *other = fields
                x1, x2, y1, y2 = list(map(int, [x1, x2, y1, y2]))

                # skip inter-chromosome interaction
                if chr1 != chr2:
                    continue
                chromosome = chr1

                if not chromosome.startswith("chr"):
                    chromosome = change_chrom_names(chromosome)
                if chromosome not in interval_tree:
                    interval_tree[chromosome] = IntervalTree()

                if len(other) == 0:
                    color = PlotHiCPeaks.DEFAULT_COLOR
                else:
                    rgb = other[0].split(",")
                    rgb = list(map(int, rgb))
                    color = rgb2hex(*rgb)

                loop = self.LoopInverval(chr1, x1, x2, chr2, y1, y2, color)
                interval_tree[chromosome].add(Interval(x1, y2, loop))

        return interval_tree

    def __is_header(self, line):
        fields = line.split()
        for idx in [1, 2, 4, 5]:
            if not fields[idx].isdigit():
                return True
        return False
