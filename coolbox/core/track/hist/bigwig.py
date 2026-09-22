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
        intervals[0] = intervals["end"]
        intervals[1] = intervals["end"]
        intervals[2] = intervals["start"]
        intervals[3] = intervals["start"]
        intervals = intervals.melt(
            id_vars="value",
            value_vars=[0, 1, 2, 3],
            var_name="corner",
            value_name="position"
        ).sort_values(by=["position", "corner"], ignore_index=True).assign(
            value = lambda df: df["value"].where(df["corner"].isin([0, 3]), 0.0)
        )
        positions = intervals["position"].values
        values = intervals['value'].values
        return positions, values

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
