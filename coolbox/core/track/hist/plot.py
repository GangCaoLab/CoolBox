import numpy as np

from coolbox.utilities import GenomeRange, get_logger

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
                self.plot_y_axis(ax, y_ax)
            except AttributeError as e:
                log.exception(e)
                msg = "If use y-axis style data range must, must set the .y_ax attribute, switch to text style."
                log.warning(msg)
                self.plot_data_range(ax, ymin, ymax, data_range_style='text')

    def __plot_range_text(self, ax, ymin, ymax):
        genome_range = self.genome_range
        ydelta = ymax - ymin

        # set min max
        format_lim = lambda lim: int(lim) if float(lim) % 1 == 0 else "{:.2f}".format(lim)
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
