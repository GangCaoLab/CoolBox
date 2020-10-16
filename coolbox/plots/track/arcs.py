import numpy as np

from coolbox.utilities import (
    change_chrom_names,
    Interval, IntervalTree,
    get_logger
)

from coolbox.plots.track.base import TrackPlot


log = get_logger(__name__)


class PlotArcs(TrackPlot):

    def __init__(self, *args, **kwarg):
        TrackPlot.__init__(self, *args, **kwarg)
        # the file format expected is similar to file format of links in
        # circos:
        # chr1 100 200 chr1 250 300 0.5
        # where the last value is a score.

        if 'color' not in self.properties:
            self.properties['color'] = 'blue'

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 0.8

    def plot(self, ax, chrom_region, region_start, region_end):
        """
        Makes and arc connecting two points on a linear scale representing
        interactions between Hi-C bins.
        """
        self.ax = ax

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
#       max_diameter += max_diameter * 0.05
        height += height * 0.05
        log.debug("{} were arcs plotted".format(count))
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(height, 0.001)
        else:
            ax.set_ylim(-0.001, height)

        ax.set_xlim(region_start, region_end)
        log.debug('title is {}'.format(self.properties['title']))

        self.plot_label()


