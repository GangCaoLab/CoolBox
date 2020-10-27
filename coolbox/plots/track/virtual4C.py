import numpy as np

from coolbox.plots.track.base import TrackPlot
from coolbox.plots.track.bigwig import CoveragePlot
from coolbox.utilities import (
    get_logger, GenomeRange
)


log = get_logger(__name__)


class PlotVirtual4C(TrackPlot, CoveragePlot):

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
        self.plot_coverage(ax, genome_range, scores_per_bin)
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

