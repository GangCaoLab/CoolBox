import numpy as np

from coolbox.utilities import (
    GenomeRange,
    to_gr,
    get_logger,
)

from coolbox.plots.track.base import TrackPlot
from matplotlib.patches import Arc


log = get_logger(__name__)


class PlotArcs(object):

    def plot_arcs(self, ax, genome_range, intervals):
        """
        Parameters
        ----------
        intervals : List[Tuple(int, int, float)]
            List of intervals (start, end, score).
        """
        gr = to_gr(genome_range)
        max_height = self.properties.get('height', 1.0)
        alpha = self.properties.get('alpha', 1.0)
        self.__adjust_yaxis(ax, max_height)
        ax.set_xlim(gr.start, gr.end)
        color = self.properties['color']

        if len(intervals) == 0:
            return

        max_itv = max(intervals, key=lambda t: t[1]-t[0])
        max_diameter = max_itv[1] - max_itv[0]

        for itv in intervals:
            start, end, score = itv
            line_width = self.__get_linewidth(score)
            diameter = (end - start)
            height = 2 * self.__get_height(max_height, max_diameter, diameter)
            center = (start + end) / 2
            ax.plot([center], [diameter])
            arc = Arc(
                (center, 0), diameter,
                height, 0, 0, 180,
                color=color,
                alpha=alpha,
                lw=line_width,
            )
            ax.add_patch(arc)

    def __get_height(self, max_height, max_diameter, diameter):
        max_height = 0.97 * max_height
        key_ = 'diameter_to_height'
        if key_ in self.properties:
            try:
                return eval(self.properties[key_])
            except Exception:
                pass
        return max_height * diameter / max_diameter

    def __get_linewidth(self, score):
        if 'line_width' in self.properties:
            return float(self.properties['line_width'])
        elif 'score_to_width' in self.properties:
            try:
                return eval(self.properties['score_to_width'])
            except Exception:
                pass
        line_width = 0.5 * np.sqrt(score)
        return line_width

    def __adjust_yaxis(self, ax, height):
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(height, 0.001)
        else:
            ax.set_ylim(-0.001, height)


class PlotBEDPE(TrackPlot, PlotArcs):
    def __init__(self, *args, **kwarg):
        TrackPlot.__init__(self, *args, **kwarg)

        if 'color' not in self.properties:
            self.properties['color'] = 'blue'

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 0.8

    def plot(self, ax, region_chrom, region_start, region_end):
        """
        Makes and arc connecting two points on a linear scale representing
        interactions between bins.
        """
        self.ax = ax
        gr = GenomeRange(region_chrom, region_start, region_end)
        itv_df = self.fetch_intervals(gr)
        point_at = self.properties.get('point_at', 'mid')
        intervals = []
        for _, row in itv_df.iterrows():
            if point_at == 'start':
                s, e = row['start1'], row['start2']
            elif point_at == 'end':
                s, e = row['end1'], row['end2']
            else:
                s = (row['start1'] + row['end1']) // 2
                e = (row['start2'] + row['end2']) // 2
            s, e = sorted([s, e])
            intervals.append((s, e, row['score']))
        self.plot_arcs(ax, gr, intervals)
        self.plot_label()


class PlotPairs(TrackPlot, PlotArcs):

    def __init__(self, *args, **kwarg):
        TrackPlot.__init__(self, *args, **kwarg)

        if 'color' not in self.properties:
            self.properties['color'] = 'cyan'

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 0.7

    def plot(self, ax, region_chrom, region_start, region_end):
        """
        Makes and arc connecting two points on a linear scale representing
        interactions between bins.
        """
        self.ax = ax
        gr = GenomeRange(region_chrom, region_start, region_end)
        itv_df = self.fetch_intervals(gr)
        intervals = []
        for _, row in itv_df.iterrows():
            s, e = row['pos1'], row['pos2']
            s, e = sorted([s, e])
            intervals.append((s, e, 1))
        self.plot_arcs(ax, gr, intervals)
        self.plot_label()
