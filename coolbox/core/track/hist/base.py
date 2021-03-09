from typing import Union, Tuple

import numpy as np
import pandas as pd

from coolbox.utilities import GenomeRange, get_logger
from ..base import Track
from .plot import PlotHist

log = get_logger(__name__)

HistData = Union[Tuple[np.ndarray, np.ndarray], pd.DataFrame, np.ndarray]


class HistBase(Track, PlotHist):
    """
    HistBase track class

    Parameters
    ----------
    style : str, optional
        Track graph type, format {'line', 'fill', 'heatmap', 'scatter'}

    fmt : str, optional
        Line styles used when the plot style is `line`, will be passed as `fmt` parameter in `matplotlib.pyplot.plot`

    line_width : int, optional
        Value used when the plot style is `line`, will be passed as `linewidth` parameter in `matplotlib.pyplot.plot`

    size : int, optional
        Value used when the plot style is `scatter`, will be passed as `s` parameter in `matplotlib.pyplot.scatter`

    color : str, optional
        Main color

    threshold_color : str, optional
        Track's color for values greater than specified threshold.

    threshold : float, optional
        Threshold used when the plot style is `line` or `scatter`, values greater than this thresh will be
        plotted with `color = threshold_color`

    cmap: str, optional
        Cmap used when the plot type is heatmap, will be passed as `cmap` paramerter in `matplotlib.pyplot.matshow`

    alpha : float, optional
        Alpha value of plot, default 1.0

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    data_range_style : {bool, 'text', 'y-axis'}, optional
        Show_data_range or not, default True.

    max_value : {float, 'auto'}, optional
        Max value of track. 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of track. 'auto' for specify max value automatically, default 'auto'.
    """

    STYLE_LINE = "line"
    STYLE_FILL = "fill"
    STYLE_HEATMAP = "heatmap"
    STYLE_SCATTER = "scatter"

    STYLES = (STYLE_LINE, STYLE_FILL, STYLE_HEATMAP, STYLE_SCATTER)

    DEFAULT_PROPERTIES = {
        "style": STYLE_LINE,
        "fmt": "-",
        "line_width": 2.0,
        "size": 10,
        "color": "#a6cee3",
        "threshold_color": "#ff9c9c",
        "threshold": "inf",
        "cmap": "bwr",
        "alpha": 1.0,
        "orientation": None,
        "data_range_style": "y-axis",
        "min_value": "auto",
        "max_value": "auto"
    }

    def __init__(self, **kwargs):
        properties = HistBase.DEFAULT_PROPERTIES.copy()
        properties.update({
            'type': properties['style'],
            **kwargs,
        })
        super().__init__(properties)
        self.genome_range = [None, None]

    def plot(self, ax, gr: GenomeRange, **kwargs):
        log.debug("plotting {}".format(self.properties.get('file', None)))

        self.ax = ax
        # if genome range would change?
        self.genome_range = gr

        data = self.fetch_plot_data(gr, **kwargs)
        if isinstance(data, pd.DataFrame):
            # ['pos', 'score']
            if 'pos' in data:
                data = data['pos'].to_numpy(), data['score'].to_numpy()
            # ['score']
            else:
                data = np.linspace(gr.start, gr.end, len(data)), data['score'].to_numpy()
        elif isinstance(data, np.ndarray):
            # 1d or 2d ndarray
            if len(data.shape) == 1:
                data = np.linspace(gr.start, gr.end, len(data)), data
            elif len(data.shape) == 2:
                data = np.linspace(gr.start, gr.end, data.shape[1]), data
            else:
                raise ValueError("The ndarray must be in 1d or 2d format")
        elif not isinstance(data, tuple):
            # not (indexes, values)
            raise ValueError(f"Data format not supported. should be one of {HistData}")

        indexes, values = data
        self.plot_hist(ax, gr, indexes, values)
        self.plot_label()

    def fetch_plot_data(self, gr: GenomeRange, **kwargs) -> HistData:
        """

        Returns
        -------
        hist_data: Union[Tuple[np.ndarray, np.ndarray], pd.DataFrame, np.ndarray]
            data used for plotting.
        """
        return self.fetch_data(gr, **kwargs)

    def fetch_data(self, gr: GenomeRange, **kwargs) -> HistData:
        raise NotImplementedError
