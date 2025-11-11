from coolbox.utilities import (
    GenomeRange, get_logger
)
import oxbow as ox
from .base import HistBase

log = get_logger(__name__)


class BigWig(HistBase):
    """
    BigWig track

    Parameters
    ----------
    file : str
        File path of bigwig file.

    num_bins : int, optional
        Number of bins to plot the hist in current range, default 700.


    """

    DEFAULT_PROPERTIES = {
        "color": "#dfccde",
        "style": HistBase.STYLE_FILL,
        "num_bins": 700,
        "threshold": "inf"
    }

    def __init__(self, file, **kwargs):
        properties = BigWig.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(**properties)
        self.ds = ox.from_bigwig(self.properties['file'])

    def fetch_plot_data(self, gr: GenomeRange, **kwargs):
        intervals = self.fetch_data(gr, **kwargs)
        values = intervals['value'].values
        return values

    def fetch_data(self, gr: GenomeRange, **kwargs):
        """
        Parameters
        ----------
        gr : GenomeRange

        Returns
        -------
        intervals : pandas.core.frame.DataFrame
            BigWig interval table.
        """
        gr = self.check_chrom_name(gr, self.ds.chrom_names)

        intervals = self.ds.regions(str(gr)).pd()
        return intervals
