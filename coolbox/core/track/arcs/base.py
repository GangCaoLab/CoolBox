import pandas as pd

from coolbox.utilities import GenomeRange
from coolbox.core.track.base import Track
from .plot import PlotContacts


class ArcsBase(Track, PlotContacts):
    """
    Arcs(link) track.

    Parameters
    ----------
    style : str, optional
        Style of arcs plot: 'arcs', 'hicpeaks', default 'arcs'

    score_to_width : str, optional
        Mapping function of score to width, default: '0.5 + math.sqrt(score)'

    line_width : float, optional
        Width of arc line.

    open_region : bool, optional
        If specified to True, will fetch the contacts on side in the region,
        default True

    diameter_to_height : str, optional
        Mapping function of arc diameter(interval end - start) to height.
        You can specify to 'max_diameter' let all arcs has same height.
        default 'max_height * diameter / max_diameter'.

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot. default None

    line_style : str, optional
        Border line style, default 'solid'

    fill : bool, optional
        Fill center or not, default False.

    fill_color : str, optional
        Fill color, use 'bed_rgb' for specify color in file,
        default 'bed_rgb'.

    fill_alpha : float, optional
        Alpha value of fill region.
        default 0.2

    side : {'upper', 'lower', 'both'}
        Plot peak in which side of the matrix.
        NOTE: This parameters is useful only if the Cool track in matrix format.

    alpha : float, optional
        Alpha value of track, default 0.8.
    """

    STYLE_ARCS = 'arcs'
    STYLE_HICPEAKS = "hicpeaks"

    DEFAULT_PROPERTIES = {
        'style': STYLE_ARCS,
        # arcs
        'line_width': None,
        'score_to_width': '0.5 + math.sqrt(score)',
        'diameter_to_height': '(max_height - 0.5) * diameter / max_diameter',
        'orientation': None,
        'cmap': None,
        'vmin': None,
        'vmax': None,
        # hicpeaks
        "side": "both",
        # general
        "line_style": "solid",
        'open_region': True,
        'color': '#3297dc',
        'alpha': 0.8,
        "fill": False,
        "fill_color": "#3297dc",
        "fill_alpha": 0.2,
    }

    def __init__(self, **kwargs):
        properties = ArcsBase.DEFAULT_PROPERTIES.copy()
        properties.update(kwargs)
        super().__init__(properties)

    def fetch_plot_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        """

        Returns
        -------
        intervals : pandas.core.frame.DataFrame
            Can be two types:
            1: with columns: ['pos1', 'pos2', 'score'] 'score' is optional
            2: with columns: ['start1', 'end1', 'start2', 'end2', 'score'] 'score' is optional
        """
        return self.fetch_data(gr, **kwargs)

    def plot(self, ax, gr: GenomeRange, **kwargs):
        """
        Plot arc connecting two positions on a linear scale representing interactions between bins.
        """
        self.ax = ax
        df = self.fetch_plot_data(gr, **kwargs)
        self.plot_contacts(ax, gr, kwargs.get("gr2"), df)
        self.plot_label()
