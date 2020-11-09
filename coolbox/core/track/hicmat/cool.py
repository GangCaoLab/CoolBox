from ..base import Track
from .plot import PlotHiCMatrix
from .fetch import FetchHiC


class Cool(Track, PlotHiCMatrix, FetchHiC):
    """
    Cool Hi-C matrix (or triangular matrix) track.

    Parameters
    ----------
    file_ : str
        Path to bed file.

    cmap : str, optional
        Color map of hic matrix, default Cool.DEFAULT_COLOR.

    style : {'triangular', 'window', 'matrix'}, optional
        Matrix style, default 'triangular'.

    balance : bool, optional
        Show balanced matrix or not, default True

    depth_ratio : float, optional
        Depth ratio of triangular matrix, use 'full' for full depth. default 'full'.

    color_bar : bool, optional
        Show color_bar or not, default True.

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
        Track's name.
    """

    DEFAULT_COLOR = 'YlOrRd'

    def __init__(self, file_, **kwargs):

        properties_dict = {
            "file": file_,
            "cmap": Cool.DEFAULT_COLOR,
            "style": 'triangular',
            "balance": True,
            "depth_ratio": "full",
            "color_bar": True,
            "transform": False,
            "norm": 'log',
            "max_value": "auto",
            "min_value": "auto",
            "title": '',
        }
        properties_dict.update(kwargs)
        properties_dict['color'] = properties_dict['cmap']

        super().__init__(properties_dict)
        self.fetched_binsize = None

    def fetch_array(self, genome_range, genome_range2=None, balance=None, resolution='auto'):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}
            Intervals within input chromosome range.

        genome_range2 : {str, GenomeRange}, optional.

        balance : bool, optional
            balance matrix or not,
            default `self.is_balance`.

        resolution : {'auto', int}
            resolution of the data. for example 5000.
            'auto' for calculate resolution automatically.
            default 'auto'

        Return
        ------
        arr : numpy.ndarray
        """
        from coolbox.utilities.hic.wrap import CoolerWrap

        path = self.properties['file']
        if balance is None:
            balance = self.is_balance
        wrap = CoolerWrap(path, balance=balance, binsize=resolution)

        arr = wrap.fetch(genome_range, genome_range2)
        return arr

    def fetch_pixels(self, genome_range, genome_range2=None, balance=None, resolution='auto', join=True):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}
            Intervals within input chromosome range.

        genome_range2 : {str, GenomeRange}, optional.

        balance : bool, optional
            balance matrix or not,
            default `self.is_balance`.

        resolution : {'auto', int}
            resolution of the data. for example 5000.
            'auto' for calculate resolution automatically.
            default 'auto'

        join : bool
            whether to expand the bin ID columns
            into (chrom, start, end).
            default True

        Return
        ------
        pixels : pandas.core.frame.DataFrame
            Hi-C pixels table.
            The pixel table contains the non-zero upper triangle entries of the contact map.
        """
        from coolbox.utilities.hic.wrap import CoolerWrap

        path = self.properties['file']
        if balance is None:
            balance = self.is_balance
        wrap = CoolerWrap(path, balance=balance, binsize=resolution)

        pixels = wrap.fetch_pixels(genome_range, genome_range2, join=join)
        return pixels
