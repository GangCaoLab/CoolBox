import numpy as np

from coolbox.plots.track.base import TrackPlot
from coolbox.utilities import get_logger, GenomeRange

log = get_logger(__name__)


class CoveragePlot(object):
    """Mixin for plot Coverage plot(BigWig, BedGraph, BAM(coverage))."""
    def plot_coverage(self, ax, genome_range, scores_per_bin, x_values=None):

        num_bins = scores_per_bin.shape[0]
        if x_values is None:
            x_values = np.linspace(genome_range.start, genome_range.end, num_bins)

        if 'style' in self.properties and self.properties['style'] != 'fill':
            self.__plot_line_or_points(ax, scores_per_bin, x_values)
        else:
            self.__plot_fill(ax, scores_per_bin, x_values)

        ymin, ymax = self.__adjust_plot(ax, genome_range)

        if "show_data_range" in self.properties and self.properties["show_data_range"] == 'no':
            pass
        else:
            self.genome_range = genome_range
            self.plot_data_range(ax, ymin, ymax, self.properties['data_range_style'])

    def plot_data_range(self, ax, ymin, ymax, data_range_style):
        if data_range_style == 'text':
            assert hasattr(self, 'genome_range'), \
                "If use text style data range must, must set the .genome_range attribute"
            self.__plot_range_text(ax, ymin, ymax)
        else:  # 'y-axis' style
            try:
                y_ax = self.y_ax
                self.plot_y_axis(y_ax)
            except AttributeError as e:
                log.exception(e)
                msg = "If use y-axis style data range must, must set the .y_ax attribute, switch to text style."
                log.warn(msg)
                self.plot_data_range(ax, ymin, ymax, data_range_style='text')

    def __plot_range_text(self, ax, ymin, ymax):
        genome_range = self.genome_range
        ydelta = ymax - ymin

        # set min max
        format_lim = lambda lim: int(lim) if float(lim) %1 == 0 else "{:.2f}".format(lim)
        ymax_print = format_lim(ymax)
        ymin_print = format_lim(ymin)
        small_x = 0.01 * genome_range.length
        # by default show the data range
        ax.text(genome_range.start - small_x, ymax - ydelta * 0.2,
                "[ {} ~ {} ]".format(ymin_print, ymax_print),
                horizontalalignment='left',
                verticalalignment='top')

    def __plot_line_or_points(self, ax, scores_per_bin, x_values):
        if self.properties['style'].find(":") > 0:
            plot_type, size = self.properties['style'].split(":")
            try:
                size = float(size)
            except ValueError:
                raise ValueError("Invalid value: 'type = {}' in Track: {}\n"
                                 "A number was expected and found '{}'".format(
                    self.properties['style'],
                    self.properties['name'],
                    size))
        else:
            plot_type = self.properties['style']
            size = None

        alpha = self.properties.get("alpha", 1.0)
        if plot_type == 'line':
            ax.plot(x_values, scores_per_bin, '-',
                    linewidth=size, color=self.properties['color'],
                    alpha=alpha)

        elif plot_type == 'points':
            ax.plot(x_values, scores_per_bin,
                    '.', markersize=size,
                    color=self.properties['color'],
                    alpha=alpha)

        else:
            raise ValueError("Invalid: 'type = {}' in Track: {}\n".format(
                self.properties['style'], self.properties['name']))

    def __plot_fill(self, ax, scores_per_bin, x_values):
        alpha = self.properties.get('alpha', 1.0)
        if 'positive_color' not in self.properties or 'negative_color' not in self.properties:
            ax.fill_between(x_values, scores_per_bin, linewidth=0.1,
                            color=self.properties['color'],
                            facecolor=self.properties['color'],
                            alpha=alpha)
        else:
            ax.fill_between(x_values, 0, scores_per_bin, where=(scores_per_bin > 0),
                            linewidth=0.1, color=self.properties['positive_color'],
                            facecolor=self.properties['positive_color'],
                            alpha=alpha)
            ax.fill_between(x_values, 0, scores_per_bin, where=(scores_per_bin < 0),
                            linewidth=0.1, color=self.properties['negative_color'],
                            facecolor=self.properties['negative_color'],
                            alpha=alpha)

    def __adjust_plot(self, ax, genome_range):
        ax.set_xlim(genome_range.start, genome_range.end)
        ymin, ymax = ax.get_ylim()
        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
            ymax = self.properties['max_value']
        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            ymin = self.properties['min_value']

        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(ymax, ymin)
        else:
            ax.set_ylim(ymin, ymax)
        return ymin, ymax


class PlotBigWig(TrackPlot, CoveragePlot):

    DEFAULT_COLOR = "#33a02c"

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)
        import pyBigWig
        self.bw = pyBigWig.open(self.properties['file'])
        if 'color' not in self.properties:
            self.properties['color'] = PlotBigWig.DEFAULT_COLOR
        if 'data_range_style' not in self.properties:
            # choices: {'text', 'y-axis'}
            # default 'text'
            # 'y-axis' style need set the .y_ax attribute
            self.properties['data_range_style'] = 'text'

    def plot(self, ax, chrom_region, start_region, end_region):
        self.ax = ax

        genome_range = GenomeRange(chrom_region, start_region, end_region)
        log.debug("plotting {}".format(self.properties['file']))

        num_bins = self.__get_bins_num()
        self.__check_chrom_name(genome_range)
        scores_per_bin = self.__get_scores_per_bin(genome_range, num_bins)

        self.plot_coverage(ax, genome_range, scores_per_bin)
        self.plot_label()

        return ax

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
        scores_per_bin = np.zeros(num_bins)
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


class PlotABCompartment(PlotBigWig):

    DEFAULT_POSITIVE_COLOR = "#ff9c9c"
    DEFAULT_NEGATIVE_COLOR = "#66ccff"

    def __init__(self, *args, **kwargs):
        PlotBigWig.__init__(self, *args, **kwargs)
        if 'positive_color' not in self.properties:
            self.properties['positive_color'] = PlotABCompartment.DEFAULT_POSITIVE_COLOR
        if 'negative_color' not in self.properties:
            self.properties['negative_color'] = PlotABCompartment.DEFAULT_NEGATIVE_COLOR



