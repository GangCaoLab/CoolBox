import numpy as np

from coolbox.plots.track.base import TrackPlot
from coolbox.utilities import get_logger, GenomeRange

log = get_logger(__name__)


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



