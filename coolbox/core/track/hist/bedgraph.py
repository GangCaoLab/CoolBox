import numpy as np
import pandas as pd

from coolbox.utilities import (
    split_genome_range, change_chrom_names,
    GenomeRange, get_logger,
)
from coolbox.utilities.bed import tabix_query, build_bedgraph_bgz


from ..base import Track
from .plot import CoveragePlot


log = get_logger(__name__)


class BedGraph(Track, CoveragePlot):
    """
    BedGraph track.

    Parameters
    ----------
    file_ : str
        Path to bedgraph file.

    height : float, optional
        Height of track, default BigWig.DEFAULT_HEIGHT

    color : str, optional
        Track color, default BigWig.DEFAULT_COLOR

    style : str, optional
        Track graph type, format {'fill', 'line:`size`', 'points:`size`'},
        example: 'line:2', 'points:0.5'. default: 'fill'

    extra : optional

    show_data_range : bool, optional
        Show_data_range or not, default True.

    data_range_style : {'text', 'y-axis'}
        The style of the data range. default: 'y-axis'

    title : str, optional
        Label text, default ''.

    max_value : {float, 'auto'}, optional
        Max value of track, use 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of track, use 'auto' for specify min value automatically, default 'auto'.

    name : str, optional
        Track's name.
    """

    DEFAULT_COLOR = '#a6cee3'

    def __init__(self, file_, **kwargs):

        properties_dict = {
            'file': file_,
            'height': BedGraph.DEFAULT_HEIGHT,
            'color': BedGraph.DEFAULT_COLOR,
            'style': 'fill',
            'show_data_range': True,
            'data_range_style': 'y-axis',
            'title': '',
            'max_value': 'auto',
            'min_value': 'auto',
        }
        properties_dict.update(kwargs)
        self.bgz_file = build_bedgraph_bgz(file_)
        properties_dict['type'] = properties_dict['style']  # change key word

        super().__init__(properties_dict)

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

    def fetch_data(self, genome_range):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        Return
        ------
        intervals : pandas.core.frame.DataFrame
            BED graph interval table.
        """
        return self.fetch_intervals(genome_range)

    def fetch_intervals(self, genome_range):
        """
        Fetch intervals within input chromosome range.
        """
        chrom, start, end = split_genome_range(genome_range)
        gr = GenomeRange(chrom, start, end)

        rows = self.__load(gr)
        if len(rows) == 0:
            chrom = change_chrom_names(chrom)
            rows = self.__load(GenomeRange(chrom, start, end))

        intval_table = pd.DataFrame(
            rows,
            columns=['chromsome', 'start', 'end', 'score']
        )

        return intval_table

    def __load(self, genome_range):
        rows = []
        gr = genome_range
        for it in tabix_query(self.bgz_file, gr.chrom, gr.start, gr.end, split=True):
            rows.append([
                it[0], int(it[1]), int(it[2]), float(it[3])
            ])
        return rows
