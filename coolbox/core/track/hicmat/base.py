from abc import ABC

from .fetch import FetchHiC
from .plot import PlotHiCMatrix, STYLE_WINDOW
from ..base import Track

common_doc1 = """    Parameters
    ----------
    file_ : str
        Path to bed file.

    cmap : str, optional
        Color map of hic matrix, default "JuiceBoxLike".

    style : {'triangular', 'window', 'matrix'}, optional
        Matrix style, default 'window'."""

common_doc2 = """    resolution : {int, 'auto'}, optional
        Matrix resolution, default 'auto'.

    normalize : {'zscore', 'expect', 'total', False}
        Normalization method, default False.
        
    gaussian_sigma : {float, False}, optional
        Do gaussian filter(with sigma, for example 1.0) on matrix if specified. default False.
        
    process_func : {callable, str, False}, optional
        Process matrix with a user-defined function(receive a matrix, return a processed matrix). default False.

    depth_ratio : float, optional
        Depth ratio of triangular matrix, use 'full' for full depth. default 'full'.

    color_bar : {'vertical', 'horizontal', 'no'}, optional
        Color bar style. default 'vertical'.

    transform : {str, bool}, optional
        Transform for matrix, like 'log2', 'log10', default False.

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    title : str, optional
        Label text, default ''.

    max_value : {float, 'auto'}, optional
        Max value of hic matrix, use 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of hic matrix, use 'auto' for specify min value automatically, default 'auto'.

    name : str, optional
        Track's name."""

fetch_matrix = """Fetch the original matrix(not processed).

    Parameters
    ----------
    genome_range : coolbox.utilities.GenomeRange
    The genome range to fetch.

    genome_range2 : coolbox.utilities.GenomeRange, optional
    Second genome range to fetch.

    resolution : {'auto', int}
    The matrix resolution, for multi-resolution(.hic or multi-cool) file.
    Use 'auto' to infer the resolution automatically.
    default 'auto'

    Return
    ------
    matrix : np.array
    Hi-C contact matrix.
"""

hic_doc = {
    "doc1": common_doc1,
    "doc2": common_doc2,
    "fetch_matrix": fetch_matrix,
}


class HicMatBase(Track, PlotHiCMatrix, FetchHiC, ABC):
    def __init__(self, file_, **kwargs):
        properties_dict = {
            "file": file_,
            "style": STYLE_WINDOW,
            "balance": True,
            "resolution": "auto",
            "normalize": False,
            "gaussian_sigma": False,
            "process_func": False,
            "depth_ratio": "full",
            "color_bar": "vertical",
            "transform": False,
            "norm": 'log',
            "max_value": "auto",
            "min_value": "auto",
            "title": '',
        }
        properties_dict.update(kwargs)
        if 'cmap' in properties_dict:
            properties_dict['color'] = properties_dict['cmap']
        super().__init__(properties_dict)
        self.fetched_binsize = None
