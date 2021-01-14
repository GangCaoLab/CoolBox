import numpy as np

from coolbox.utilities import get_logger
from coolbox.utilities.genome import GenomeRange

log = get_logger(__name__)


class PlotHist(object):
    """Mixin for plot Coverage plot(BigWig, BedGraph, BAM(coverage))."""

    def plot_hist(self, ax, gr: GenomeRange, indexes: np.ndarray, values: np.ndarray):
        assert len(indexes.shape) == 1 and 1 <= len(values.shape) <= 2
        style = self.properties.get('style', 'line')
        if 'heatmap' in style:
            self.plot_heatmap(ax, gr, values)
        elif 'fill' in style:
            self.plot_fill(ax, indexes, values)
        elif 'scatter' in style:
            self.plot_scatter(ax, indexes, values)
        elif 'line' in style:
            self.plot_line(ax, indexes, values)

        ymin, ymax = self.adjust_plot(ax, gr)
        # disable plot data range in coverage mode
        if self.properties['data_range_style'] != 'no' and 'track' not in self.__dict__:
            self.plot_data_range(ax, ymin, ymax, self.properties['data_range_style'], gr)

    def plot_heatmap(self, ax, gr: GenomeRange, mat: np.ndarray):
        # no threshold supported
        cmap = self.properties.get('cmap', None)
        if len(mat.shape) == 1:
            height = self.properties.get('height', 3)
            mat = np.array([mat] * height)
        ax.matshow(mat, cmap=cmap, aspect="auto", extent=(gr.start, gr.end, 0, mat.shape[0]))

    def plot_fill(self, ax, indexes, values):
        alpha = self.properties.get('alpha', 1.0)
        threshold = self.properties.get("threshold", 0)
        ax.fill_between(indexes, 0, values, where=(values > threshold),
                        linewidth=0.1, color=self.properties['threshold_color'],
                        facecolor=self.properties['threshold_color'],
                        alpha=alpha)
        ax.fill_between(indexes, 0, values, where=(values < threshold),
                        linewidth=0.1, color=self.properties['color'],
                        facecolor=self.properties['color'],
                        alpha=alpha)

    def plot_scatter(self, ax, indexes: np.ndarray, values: np.ndarray):
        threshold = self.properties.get("threshold", float("inf"))
        color = self.properties['color']
        size = self.properties['size']
        alpha = self.properties['alpha']
        mask = values > threshold
        ax.scatter(indexes[mask], values[mask], s=size, alpha=alpha, c=self.properties.get('threshold_color', color))
        mask = ~mask
        ax.scatter(indexes[mask], values[mask], s=size, alpha=alpha, c=color)

    def plot_line(self, ax, indexes, values):
        # reference for plotting with threshold: https://stackoverflow.com/a/30122991/10336496
        # if the data is 2d matrix, transpose to comply wtih matplotlib
        if len(values.shape) == 2:
            values = values.T
        fmt = self.properties.get("fmt")
        line_width = self.properties.get('line_width', 1)
        color = self.properties.get('color')
        alpha = self.properties.get("alpha", 1.0)
        threshold = float(self.properties.get("threshold"))
        threshold_color = self.properties.get("threshold_color")
        ax.plot(indexes, values, fmt, linewidth=line_width, color=color, alpha=alpha)
        if threshold and np.sum(values > threshold) > 0:
            masked_values = np.ma.masked_greater_equal(values, threshold)
            ax.plot(indexes, masked_values, fmt, linewidth=line_width, color=threshold_color, alpha=alpha)

    def plot_data_range(self, ax, ymin, ymax, data_range_style, gr: GenomeRange):
        if data_range_style == 'text':
            self.plot_text_range(ax, ymin, ymax, gr)
        else:  # 'y-axis' style
            try:
                y_ax = self.y_ax
                self.plot_yaxis_range(ax, y_ax)
            except AttributeError as e:
                log.exception(e)
                msg = "If use y-axis style data range must, must set the .y_ax attribute, switch to text style."
                log.warning(msg)
                self.plot_data_range(ax, ymin, ymax, 'text', gr)

    def plot_yaxis_range(self, plot_axis, y_ax):
        """
        Plot the scale of the y axis with respect to the plot_axis

        plot something that looks like this:
        ymax ┐
             │
             │
        ymin ┘

        Parameters
        ----------
        plot_axis : matplotlib.axes.Axes
            Main plot axis.

        y_ax : matplotlib.axes.Axes
            Axis to use to plot the scale
        """

        if 'show_data_range' in self.properties and self.properties['show_data_range'] == 'no':
            return

        def value_to_str(value):
            if value % 1 == 0:
                str_value = str(int(value))
            else:
                if value < 0.01:
                    str_value = "{:.4f}".format(value)
                else:
                    str_value = "{:.2f}".format(value)
            return str_value

        ymin, ymax = plot_axis.get_ylim()

        ymax_str = value_to_str(ymax)
        ymin_str = value_to_str(ymin)
        x_pos = [0, 0.5, 0.5, 0]
        y_pos = [0.01, 0.01, 0.99, 0.99]
        y_ax.plot(x_pos, y_pos, color='black', linewidth=1, transform=y_ax.transAxes)
        y_ax.text(-0.2, -0.01, ymin_str, verticalalignment='bottom', horizontalalignment='right',
                  transform=y_ax.transAxes)
        y_ax.text(-0.2, 1, ymax_str, verticalalignment='top', horizontalalignment='right', transform=y_ax.transAxes)
        y_ax.patch.set_visible(False)

    def plot_text_range(self, ax, ymin, ymax, gr: GenomeRange):
        ydelta = ymax - ymin

        # set min max
        format_lim = lambda lim: int(lim) if float(lim) % 1 == 0 else "{:.2f}".format(lim)
        ymax_print = format_lim(ymax)
        ymin_print = format_lim(ymin)
        small_x = 0.01 * gr.length
        # by default show the data range
        ax.text(gr.start - small_x, ymax - ydelta * 0.2,
                "[ {} ~ {} ]".format(ymin_print, ymax_print),
                horizontalalignment='left',
                verticalalignment='top')

    def adjust_plot(self, ax, gr: GenomeRange):
        ax.set_xlim(gr.start, gr.end)
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
