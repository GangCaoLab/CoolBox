import pandas as pd

from .base import ArcsBase
from .fetch import FetchParix
from coolbox.utilities.bed import process_pairs
from coolbox.utilities.genome import GenomeRange


class Pairs(ArcsBase, FetchParix):
    """
    Arcs track from .pairs file.

    Parameters
    ----------
    file: str
        Path of .pairs file

    """
    DEFAULT_PROPERTIES = {
        'color': '#dc9732'
    }
    FIELDS = ["name", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]

    def __init__(self, file, **kwargs):
        properties = Pairs.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(**properties)
        self.bgz_file = process_pairs(file)

    def fetch_data(self, gr: GenomeRange, **kwargs):
        # filter peaks manually in peaks style
        df = self.fetch_intervals(self.bgz_file, gr, kwargs.get('gr2'))
        # TODO the returned df has no named columns, may cause error
        if len(df) == 0:
            return df

        columns = list(df.columns)
        for i, col in enumerate(self.FIELDS):
            if i >= len(columns):
                break
            columns[i] = col
        df.columns = columns
        for col in ['pos1', 'pos2']:
            df[col] = df[col].astype(int)
        return df
