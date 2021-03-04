import pandas as pd

from .base import ArcsBase
from .fetch import FetchParix
from coolbox.utilities.bed import process_bedpe
from coolbox.utilities.genome import GenomeRange


class BEDPE(ArcsBase, FetchParix):
    """
    Arcs track from .bedpe file.

    Parameters
    ----------
    file: str
        Path of .bedpe file

    pos : str, optional
        Method for choosing arch anchor for bedpe data: 'start', 'end', 'mid', default 'mid'

    """
    DEFAULT_PROPERTIES = {
        'pos': 'mid'
    }
    FIELDS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2",
              "name", "score", "strand1", "strand2"]

    def __init__(self, file, **kwargs):
        properties = BEDPE.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(**properties)
        self.bgz_file = process_bedpe(file)

    def fetch_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        # filter peaks manually for hicpeaks style in fetch_plot_data
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
        for col in ['start1', 'end1', 'start2', 'end2']:
            df[col] = df[col].astype(int)
        if 'score' in df.columns:
            df['score'] = df['score'].astype(float)
        return df

    def fetch_plot_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        df = self.fetch_data(gr, **kwargs)
        if len(df) == 0:
            return df

        style = self.properties['style']
        if style == self.STYLE_ARCS:
            pos_at = self.properties['pos']
            if pos_at == 'start':
                pos1, pos2 = df['start1'], df['start2']
            elif pos_at == 'end':
                pos1, pos2 = df['end1'], df['end2']
            else:
                pos1, pos2 = (df['end1'] + df['start1']) / 2, (df['end2'] + df['start2']) / 2
            return pd.DataFrame({'pos1': pos1, 'pos2': pos2, 'score': df['score']})
        elif style == self.STYLE_HICPEAKS:
            gr2 = kwargs.get('gr2')
            if gr2 and gr2 != gr:
                mask = (df['start2'] >= gr2.start) & (df['end2'] <= gr2.end)
                df = df[mask]
            return df
        else:
            raise ValueError("The supported style for bedpe data are ['arcs', 'hicpeaks']")
