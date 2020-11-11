from ..base import Track
from .plot import PlotHiCMatrix
from .fetch import FetchHiC


class DotHiC(Track, PlotHiCMatrix, FetchHiC):

    """
    .hic Hi-C matrix (or triangular matrix) track.

    Parameters
    ----------
    file_ : str
        Path to bed file.

    cmap : str, optional
        Color map of hic matrix, default "JuiceBoxLike".

    style : {'triangular', 'window', 'matrix'}, optional
        Matrix style,
        default 'window'.

    balance : {bool, 'KR', 'VC', 'VC_SQRT'}, optional
        Matrix balance method,
        default True('KR' balance)

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
        Track's name.

    """
    def __init__(self, file_, **kwargs):

        properties_dict = {
            "file": file_,
            "cmap": "JuiceBoxLike",
            "style": 'window',
            "balance": True,
            "depth_ratio": "full",
            "color_bar": "vertical",
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

        balance : {bool, 'KR', 'VC', 'VC_SQRT'}, optional
            matrix balance method,
            default `self.balance`.

        resolution : {'auto', int}
            resolution of the data. for example 5000.
            'auto' for calculate resolution automatically.
            default 'auto'

        Return
        ------
        arr : numpy.ndarray
        """
        from coolbox.utilities.hic.wrap import StrawWrap

        path = self.properties['file']
        if balance is None:
            balance = self.balance
        wrap = StrawWrap(path, normalization=balance, binsize=resolution)

        arr = wrap.fetch(genome_range, genome_range2)
        return arr

    def fetch_pixels(self, genome_range, genome_range2=None, balance=None, resolution='auto'):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}
            Intervals within input chromosome range.

        balance : {bool, 'KR', 'VC', 'VC_SQRT'}, optional
            matrix balance method,
            default `self.balance`.

        resolution : {'auto', int}
            resolution of the data. for example 5000.
            'auto' for calculate resolution automatically.
            default 'auto'

        Return
        ------
        pixels : pandas.core.frame.DataFrame
        """
        from coolbox.utilities.hic.wrap import StrawWrap

        path = self.properties['file']
        if balance is None:
            balance = self.balance
        wrap = StrawWrap(path, normalization=balance, binsize=resolution)

        pixels = wrap.fetch_pixels(genome_range, genome_range2)
        return pixels


