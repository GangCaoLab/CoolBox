from ..base import Track
from .plot import PlotHiCMatrix
from .fetch import FetchHiC

from coolbox.utilities.doctool import paste_doc
from coolbox.utilities import to_gr
from .base import hic_doc


@paste_doc(hic_doc)
class DotHiC(Track, PlotHiCMatrix, FetchHiC):
    """
    .hic Hi-C matrix (or triangular matrix) track.

${doc1}

    balance : {bool, 'KR', 'VC', 'VC_SQRT'}, optional
        Matrix balance method,
        default True('KR' balance)

${doc2}
    """
    def __init__(self, file_, **kwargs):

        properties_dict = {
            "file": file_,
            "cmap": "JuiceBoxLike2",
            "style": 'window',
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
        properties_dict['color'] = properties_dict['cmap']

        super().__init__(properties_dict)
        self.fetched_binsize = None

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
            Hi-C pixels table.
            The pixel table contains the non-zero upper triangle entries of the contact map.
        """
        from coolbox.utilities.hic.wrap import StrawWrap

        genome_range = to_gr(genome_range)
        if genome_range2 is not None:
            genome_range2 = to_gr(genome_range2)

        path = self.properties['file']
        if balance is None:
            balance = self.balance
        wrap = StrawWrap(path, normalization=balance, binsize=resolution)

        pixels = wrap.fetch_pixels(genome_range, genome_range2)
        return pixels


