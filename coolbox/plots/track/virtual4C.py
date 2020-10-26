import numpy as np

from coolbox.plots.track.base import TrackPlot
from coolbox.utilities import (
    get_logger, GenomeRange
)


log = get_logger(__name__)


class PlotVirtual4C(TrackPlot):

    DEFAULT_COLOR = '#2855d8'

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)
        self.hic = self.properties['hic']
        self.position = GenomeRange(self.properties['genome_position'])
        self.bin_width = self.properties['bin_width']

        if 'color' not in self.properties:
            self.properties['color'] = PlotVirtual4C.DEFAULT_COLOR
        if 'data_range_style' not in self.properties:
            # choices: {'text', 'y-axis'}
            # default 'text'
            # 'y-axis' style need set the .y_ax attribute
            self.properties['data_range_style'] = 'text'
        self.properties['type'] = self.properties['style']

    def plot(self, ax, chrom_region, start_region, end_region):
        self.ax = ax
        genome_range = GenomeRange(chrom_region, start_region, end_region)
        scores_per_bin = self.fetch_mean_arr(genome_range)
        num_bins = scores_per_bin.shape[0]

        x_values = np.linspace(genome_range.start, genome_range.end, num_bins)

        if 'type' in self.properties and self.properties['type'] != 'fill':
            self.__plot_line_or_points(scores_per_bin, x_values)
        else:
            self.__plot_fill(scores_per_bin, x_values)

        ymin, ymax = self.__adjust_plot(genome_range)

        if "show_data_range" in self.properties and self.properties["show_data_range"] == 'no':
            pass
        else:
            self.genome_range = genome_range
            self.plot_data_range(ymin, ymax, self.properties['data_range_style'])

        self.plot_label()

    def fetch_mean_arr(self, genome_range):
        from copy import copy
        bin_width = self.bin_width
        position = self.position
        binsize = self.hic.fetched_binsize
        assert binsize is not None, "Related HiC track should plot firstly"
        window_range = copy(position)
        offset_ = (bin_width - 1) // 2
        assert offset_ >= 0, "bin width must >= 1"
        window_range.start = window_range.start - offset_ * binsize
        window_range.end = window_range.end + offset_ * binsize
        arr = self.hic.fetch_array(window_range, genome_range)
        #mean_arr = arr.mean(axis=0)
        mean_arr = np.nanmean(arr, axis=0)
        return mean_arr

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

    def plot_data_range(self, ymin, ymax, data_range_style):

        if data_range_style == 'text':
            assert hasattr(self, 'genome_range'), \
                "If use text style data range must, must set the .genome_range attribute"
            self.__plot_range_text(ymin, ymax)

        else:  # 'y-axis' style
            try:
                y_ax = self.y_ax
                self.plot_y_axis(y_ax)
            except AttributeError as e:
                log.exception(e)
                msg = "If use y-axis style data range must, must set the .y_ax attribute, switch to text style."
                log.warn(msg)
                self.plot_data_range(ymin, ymax, data_range_style='text')

    def __plot_range_text(self, ymin, ymax):
        genome_range = self.genome_range
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
