from ..base import Track
from .plot import PlotHiCMatrix
from .fetch import FetchHiC

from coolbox.utilities.doctool import paste_doc
from coolbox.utilities import to_gr
from .base import hic_doc


@paste_doc(hic_doc)
class Cool(Track, PlotHiCMatrix, FetchHiC):
    """
    Cool Hi-C matrix (or triangular matrix) track.

${doc1}

    balance : bool, optional
        Show balanced matrix or not, default True

${doc2}
    """

    DEFAULT_COLOR = 'JuiceBoxLike'

    def __init__(self, file_, **kwargs):

        properties_dict = {
            "file": file_,
            "cmap": Cool.DEFAULT_COLOR,
            "style": 'window',
            "balance": True,
            "resolution": "auto",
            "normalize": False,
            "gaussian_sigma": False,
            "process_func": False,
            "depth_ratio": "full",
            "color_bar": 'vertical',
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

        genome_range = to_gr(genome_range)
        if genome_range2 is not None:
            genome_range2 = to_gr(genome_range2)

        path = self.properties['file']
        if balance is None:
            balance = self.is_balance
        wrap = CoolerWrap(path, balance=balance, binsize=resolution)

        pixels = wrap.fetch_pixels(genome_range, genome_range2, join=join)
        return pixels
