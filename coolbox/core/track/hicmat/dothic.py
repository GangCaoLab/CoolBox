import numpy as np

from coolbox.utilities import to_gr
from coolbox.utilities.doctool import paste_doc
from .base import hic_doc, HicMatBase


@paste_doc(hic_doc)
class DotHiC(HicMatBase):
    """
    .hic Hi-C matrix (or triangular matrix) track.

${doc1}

    balance : {bool, 'KR', 'VC', 'VC_SQRT'}, optional
        Matrix balance method,
        default True('KR' balance)

${doc2}
    """
    DEFAULT_COLOR = "JuiceBoxLike2"

    def __init__(self, file_, **kwargs):
        properties_dict = {
            "cmap": self.DEFAULT_COLOR,
        }
        properties_dict.update(kwargs)

        super().__init__(file_, **properties_dict)

    def fetch_pixels(self, genome_range, genome_range2=None, balance=None, resolution='auto'):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}
            Intervals within input chromosome range.

        genome_range2 : {str, GenomeRange}
            Intervals within input chromsome range2.

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

    @paste_doc(hic_doc)
    def fetch_matrix(self, genome_range, genome_range2=None, resolution='auto') -> np.ndarray:
        """
        ${fetch_matrix}
        """
        from coolbox.utilities.hic.wrap import StrawWrap

        path = self.properties['file']
        wrap = StrawWrap(path, normalization=self.balance, binsize=resolution)

        arr = wrap.fetch(genome_range, genome_range2)

        self.fetched_binsize = wrap.fetched_binsize  # expose fetched binsize

        return self.fill_zero_nan(arr)

    def _infer_binsize(self, genome_range1, genome_range2=None, resolution=None) -> int:
        from coolbox.utilities.hic.wrap import StrawWrap

        path = self.properties['file']
        wrap = StrawWrap(path, normalization=self.balance, binsize=resolution)
        return wrap.infer_binsize(genome_range1)
