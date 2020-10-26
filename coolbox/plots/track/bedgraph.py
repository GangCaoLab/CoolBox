from coolbox.plots.track.base import TrackPlot
from coolbox.utilities import file_to_intervaltree, change_chrom_names, get_logger

log = get_logger(__name__)


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

    def plot(self, ax, chrom_region, start_region, end_region):
        self.ax = ax

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

        if 'show_data_range' in self.properties and self.properties['show_data_range'] == 'no':
            pass
        else:
            # by default show the data range
            self.plot_data_range(ymin, ymax, self.properties['data_range_style'])

        self.plot_label()

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
