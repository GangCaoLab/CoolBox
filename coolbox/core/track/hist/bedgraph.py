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
    -----------
    file : str
        File path of bedgraph file.

    """

    DEFAULT_PROPERTIES = {
        "style": HistBase.STYLE_FILL,
    }

    def __init__(self, file, **kwargs):
        properties = BedGraph.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(**properties)
        self.bgz_file = build_bedgraph_bgz(file)

    def fetch_plot_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        itv_df = self.fetch_data(gr, **kwargs)
        index_array = itv_df['start'] + (itv_df['end'] - itv_df['start']) / 2
        itv_df['pos'] = index_array
        return itv_df

    def fetch_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        rows = self.load(gr)
        if len(rows) == 0:
            gr.chrom = change_chrom_names(gr.chrom)
            rows = self.load(gr)

        intval_table = pd.DataFrame(rows, columns=['chromsome', 'start', 'end', 'score'])

        return intval_table

    def load(self, genome_range):
        rows = []
        gr = genome_range
        for it in tabix_query(self.bgz_file, gr.chrom, gr.start, gr.end, split=True):
            rows.append([
                it[0], int(it[1]), int(it[2]), float(it[3])
            ])
        return rows
