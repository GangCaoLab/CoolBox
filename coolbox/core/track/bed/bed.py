from coolbox.core.track.bed.fetch import FetchBed
from coolbox.utilities import (
    get_logger
)
from coolbox.utilities.genome import GenomeRange
from .base import BedBase
from .plot import PlotGenes

log = get_logger(__name__)


class BED(BedBase, PlotGenes, FetchBed):
    """
    Bed Track for plotting 1d intervals data from .bed file.
    The input bed file can be bed3/bed6/bed9/bed12

    Parameters
    ----------
    gene_style: {'flybase', 'normal'}

    display : {'stacked', 'interlaced', 'collapsed'}, optional
        Display mode. (Default: 'stacked')

    fontsize : int, optional
        Font size. (Default: BED.DEFAULT_FONTSIZE)

    labels : {True, False, 'auto'}, optional
        Draw bed name or not. 'auto' for automate decision according to density.
        (Default: 'auto')

    interval_height : int, optional
        The height of the interval. (Default: 100)

    num_rows : int, optional
        Set the max interval rows. (Default: unlimited interval rows)

    row_height : float
        Height of a row. default 0.5
    """

    DEFAULT_PROPERTIES = {
        'labels': 'auto',
        'height': 'auto',
        'gene_style': 'flybase',
        'display': 'stacked',
        'fontsize': 12,
        'interval_height': 100,
        'num_rows': None,
        'row_height': 0.5,
    }

    def __init__(self, file, **kwargs):
        properties = BED.DEFAULT_PROPERTIES.copy()
        properties.update(kwargs)
        super().__init__(file, **properties)
        PlotGenes.__init__(self)

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax
        ov_intervals: pd.DataFrame = self.fetch_plot_data(gr, **kwargs)
        self.plot_genes(ax, gr, ov_intervals)
        self.plot_label()
