from typing import Union, Tuple

import numpy as np
import pandas as pd

from coolbox.utilities import (
    split_genome_range, change_chrom_names,
    GenomeRange, get_logger,
)
from coolbox.utilities.bed import tabix_query, build_bedgraph_bgz
from .base import HistBase

log = get_logger(__name__)


class BedGraph(HistBase):
    """
    BedGraph track.

    Parameters
    ----------
    file_ : str
        Path to bedgraph file.

    height : float, optional
        Height of track, default BedGraph.DEFAULT_HEIGHT

    color : str, optional
        Track color, default BedGraph.DEFAULT_COLOR

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
            'color': BedGraph.DEFAULT_COLOR,
            'style': 'fill',
        }
        self.bgz_file = build_bedgraph_bgz(file_)
        self.genome_range = None
        properties_dict.update(kwargs)
        super().__init__(**properties_dict)

    def fetch_data(self, genome_range: Union[str, GenomeRange]) -> pd.DataFrame:
        chrom, start, end = split_genome_range(genome_range)
        gr = GenomeRange(chrom, start, end)
        rows = self.__load(gr)
        if len(rows) == 0:
            chrom = change_chrom_names(chrom)
            rows = self.__load(GenomeRange(chrom, start, end))

        intval_table = pd.DataFrame(rows, columns=['chromsome', 'start', 'end', 'score'])

        return intval_table

    def fetch_plot_data(self, genome_range: GenomeRange) -> Tuple[np.ndarray, np.ndarray]:
        itv_df = self.fetch_data(genome_range)
        itv_df['pos'] = itv_df['start'] + (itv_df['end'] - itv_df['start']) / 2
        return np.array(itv_df['score']), np.array(itv_df['pos'])

    def __load(self, genome_range):
        rows = []
        gr = genome_range
        for it in tabix_query(self.bgz_file, gr.chrom, gr.start, gr.end, split=True):
            rows.append([
                it[0], int(it[1]), int(it[2]), float(it[3])
            ])
        return rows
