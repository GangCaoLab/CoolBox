from typing import Tuple

import numpy as np
import pandas as pd

from .process import ProcessHicMat
from .plot import PlotHiCMat
from ..base import Track
from coolbox.utilities.genome import GenomeRange


class HicMatBase(Track, PlotHiCMat, ProcessHicMat):
    """
    HicMatBase class for all track plotting matrix-like data.

    Parameters
    ----------
    style : {'triangular', 'window', 'matrix'}, optional
        Matrix style, default 'window'.

    cmap : str, optional
        Color map of hic matrix, default "JuiceBoxLike".

    color_bar : {'vertical', 'horizontal', 'no'}, optional
        Color bar style. default 'vertical'.

    depth_ratio : float, optional
        Depth ratio of triangular matrix, use 'full' for full depth. default 'full'.

    max_value : {float, 'auto'}, optional
        Max value of hic matrix, use 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of hic matrix, use 'auto' for specify min value automatically, default 'auto'.

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    resolution : {int, 'auto'}, optional
        Matrix resolution, default 'auto'.

    transform : {str, bool}, optional
        Transform for matrix, like 'log2', 'log10', default False.

    normalize : {'zscore', 'expect', 'total', False}
        Normalization method, default False.

    norm: {'log', 'no'}
        Method used when normalizing. default 'no'

    gaussian_sigma : {float, False}, optional
        Do gaussian filter(with sigma, for example 1.0) on matrix if specified. default False.

    process_func : {callable, str, False}, optional
        Process matrix with a user-defined function(receive a matrix, return a processed matrix). default False.


    """

    STYLE_TRIANGULAR = 'triangular'
    STYLE_MATRIX = 'matrix'
    STYLE_WINDOW = 'window'

    SMALL_VALUE = 1e-12
    DEPTH_FULL = 'full'

    DEFAULT_PROPERTIES = {
        "style": STYLE_WINDOW,
        # processing
        "resolution": "auto",
        "transform": 'log',
        "normalize": False,
        "norm": False,
        "gaussian_sigma": False,
        "process_func": False,
        # plotting
        'height': 'hic_auto',
        'cmap': "JuiceBoxLike",
        "color_bar": "vertical",
        "max_value": "auto",
        "min_value": "auto",
        "depth_ratio": DEPTH_FULL,
        "orientation": None,
    }

    def __init__(self, **kwargs):
        properties = HicMatBase.DEFAULT_PROPERTIES.copy()
        properties.update(kwargs)
        if 'cmap' in properties:
            properties['color'] = properties['cmap']
        super().__init__(properties)
        self.fetched_binsize = None
        self.fetched_gr = None
        self.fetched_gr2 = None

    def fetch_data(self, gr: GenomeRange, **kwargs) -> np.ndarray:
        """
        Fetch the raw matrix should be plotted. Normally it's a matrix with raw contacts

        Parameters
        ----------
        gr2 : GenomeRange, optional, keyword argument

        Returns
        -------
        matrix : np.array
            Hi-C raw contact matrix.

        """
        raise NotImplementedError

    def fetch_plot_data(self, gr: GenomeRange, **kwargs) -> np.ndarray:
        """
        Fetch 2d contact matrix, the matrix may be processed in case
        'transform', 'normalize', 'gaussian_sigma', 'process_func' exits in properties.

        Parameters
        ----------
        gr2 : GenomeRange, optional

        gr_updated: bool, optional
            If the input GenomeRange has been updated. default False
            Default False means that the input gr will be expanded in window mode

        Returns
        -------
        matrix : np.array
            Processed hic matrix used for plotting.
        """
        gr2 = kwargs.get('gr2')
        if self.properties['style'] == self.STYLE_WINDOW and not kwargs.get("gr_updated", False):
            gr, gr2 = self.fetch_window_genome_range(gr, gr2)
            kwargs.update({'gr2': gr2})
        arr = self.fetch_data(gr, **kwargs)
        # store fetched_gr
        self.fetched_gr = gr
        self.fetched_gr2 = gr2
        return self.process_matrix(arr)

    def plot(self, ax, gr: GenomeRange, **kwargs):
        """
        Plot matrix

        Parameters
        ----------
        ax
        gr2 : GenomeRange, optional
        """
        self.ax = ax
        # fetch processed plot_data
        self.matrix = self.fetch_plot_data(gr, **kwargs)
        # plot matrix
        img = self.plot_matrix(gr, kwargs.get('gr2'))
        self.adjust_figure(gr, kwargs.get('gr2'))
        self.draw_colorbar(img)
        self.plot_label()

    def fetch_pixels(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        raise NotImplementedError

    def infer_binsize(self, gr: GenomeRange, **kwargs) -> int:
        raise NotImplementedError

    def fetch_window_genome_range(self, gr: GenomeRange, gr2: GenomeRange = None) -> Tuple[GenomeRange, GenomeRange]:
        from copy import copy
        fetch_gr = copy(gr)
        dr = self.properties['depth_ratio']
        dr = 1.0 if dr == "full" else dr
        dr = min(1.0, dr + 0.05)
        x = int(gr.length * dr // 2)
        fetch_gr.start = gr.start - x
        fetch_gr.end = gr.end + x

        if fetch_gr.start < 0:
            fetch_gr.start = 0

        return fetch_gr, gr2
