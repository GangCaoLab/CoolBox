import pandas as pd

from coolbox.utilities import GenomeRange, get_logger
from coolbox.utilities.reader.tab import get_indexed_tab_reader
from .base import HistBase

log = get_logger(__name__)


class BedGraph(HistBase):
    """
    BedGraph track.

    Parameters
    ----------
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
        self.reader = get_indexed_tab_reader(file)

    def fetch_plot_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        itv_df = self.fetch_data(gr, **kwargs)
        index_array = itv_df['start'] + (itv_df['end'] - itv_df['start']) / 2
        itv_df['pos'] = index_array
        itv_df['score'] = itv_df.pop('value')
        return itv_df

    def fetch_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        return self.reader.query_var_chr(gr)
