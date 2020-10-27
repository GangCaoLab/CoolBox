import numpy as np

from coolbox.plots.track.base import TrackPlot
from coolbox.plots.track.bigwig import CoveragePlot
from coolbox.utilities import (
    get_logger, GenomeRange
)

log = get_logger(__name__)


class PlotBedGraph(TrackPlot, CoveragePlot):

    DEFAULT_COLOR = '#a6cee3'

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)
        if 'color' not in self.properties:
            self.properties['color'] = PlotBedGraph.DEFAULT_COLOR

    def plot(self, ax, chrom_region, start_region, end_region):
        self.ax = ax
        genome_range = GenomeRange(chrom_region, start_region, end_region)
        self.genome_range = genome_range

        itv_df = self.fetch_intervals(genome_range)

        itv_df['pos'] = itv_df['start'] + (itv_df['end'] - itv_df['start'])/2

        self.plot_coverage(ax, genome_range,
                           np.array(itv_df['score']),
                           np.array(itv_df['pos']))
        self.plot_label()
        return ax

