import pandas as pd

from coolbox.utilities import (
    get_logger
)
from coolbox.utilities.genome import GenomeRange

from coolbox.core.track.base import Track
from .plot import PlotBed

log = get_logger(__name__)


class BedBase(Track, PlotBed):
    """
    BED Base track.

    Parameters
    ----------
    style : {'gene', 'tad'}

    gene_style: {'flybase', 'normal'}

    display : {'stacked', 'interlaced', 'collapsed'}, optional
        Display mode. (Default: 'stacked')

    color : str, optional
        Track color, 'bed_rgb' for auto specify color according to bed record.
        (Default: 'bed_rgb')

    border_color : str, optional
        Border_color of gene. (Default: 'black')

    fontsize : int, optional
        Font size. (Default: BED.DEFAULT_FONTSIZE)

    labels : {True, False, 'auto'}, optional
        Draw bed name or not. 'auto' for automate decision according to density.
        (Default: 'auto')

    interval_height : int, optional
        The height of the interval. (Default: 100)

    num_rows : int, optional
        Set the max interval rows. (Default: unlimited interval rows)

    max_value : float, optional
        Max score. (Default: inf)

    min_value : float, optional
        Min score. (Default: -inf)

    border_style: str, optional
        Border style of tad. (Default: 'solid')

    border_width: int, optional
        Border width of tad. (Default: '2.0')

    show_score : bool
        Show bed score or not.
        default False.

    score_font_size : {'auto', int}
        Score text font size.
        default 'auto'

    score_font_color : str
        Score text color.
        default '#000000'

    score_height_ratio : float
        (text tag height) / (TAD height). used for adjust the position of Score text.
        default 0.5

    border_only : bool
        Only show border, default False

    """

    STYLE_GENE = "gene"
    STYLE_TAD = "tad"

    COLOR = "#1f78b4"

    DEFAULT_PROPERTIES = {
        'style': STYLE_GENE,
        # gene
        'gene_style': 'flybase',
        'display': 'stacked',
        'color': "bed_rgb",
        'border_color': "#1f78b4",
        'fontsize': 12,
        'interval_height': 100,
        'num_rows': None,
        'labels': 'off',
        'min_score': '-inf',
        'max_score': 'inf',
        'bed_type': None,
        # tad
        'border_style': "--",
        'border_width': 2.0,
        "show_score": False,
        "score_font_size": 'auto',
        "score_font_color": "#000000",
        "score_height_ratio": 0.4,
        "border_only": False,
    }

    def __init__(self, **kwargs):
        properties = BedBase.DEFAULT_PROPERTIES.copy()
        properties.update(kwargs)
        super().__init__(properties)
        self.init_for_plot()

    def fetch_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        """

        Returns
        -------
        intervals : pandas.core.frame.DataFrame
            BED interval table. The table should be in format like:

            bed_fields = ['chromosome', 'start', 'end',
                          'name', 'score', 'strand',
                          'thick_start', 'thick_end',
                          'rgb', 'block_count',
                          'block_sizes', 'block_starts']
            The table can be in bed6/bed9/bed12 format and the trailing columns can be omited.

        """
        raise NotImplementedError

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax
        ov_intervals: pd.DataFrame = self.fetch_plot_data(gr, **kwargs)

        style = self.properties['style']
        if style == self.STYLE_TAD:
            self.plot_tads(ax, gr, ov_intervals)
        elif style == self.STYLE_GENE:
            self.plot_genes(ax, gr, ov_intervals)
        else:
            raise ValueError("style not supportted, should be one of 'gene' 'tad' ")
        self.plot_label()
