import numpy as np

from coolbox.utilities import to_gr
from coolbox.utilities.doctool import paste_doc
from .base import HicMatBase


class DotHiC(HicMatBase):
    """HicMat track from .hic file.

    Parameters
    ----------
    file: str
        The file path of .hic file.

    balance : {bool, 'KR', 'VC', 'VC_SQRT'}, optional
        Matrix balance method,
        default True('KR' balance)


    """
    DEFAULT_PROPERTIES = {
        'cmap': "JuiceBoxLike2",
        'balance': True,
    }

    def __init__(self, file, **kwargs):
        properties = DotHiC.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(**properties)

    def fetch_data(self, genome_range, genome_range2=None, **kwargs) -> np.ndarray:
        from coolbox.utilities.hic.wrap import StrawWrap

        path = self.properties['file']
        wrap = StrawWrap(path, normalization=self.balance, binsize=kwargs.get('resolution', 'auto'))

        arr = wrap.fetch(genome_range, genome_range2)

        self.fetched_binsize = wrap.fetched_binsize  # expose fetched binsize

        return self.fill_zero_nan(arr)

    def fetch_pixels(self, genome_range, genome_range2=None, balance=None, **kwargs):
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
        balance = kwargs.get('balance', self.is_balance)
        wrap = StrawWrap(path, normalization=balance, binsize=kwargs.get('resolution', 'auto'))

        pixels = wrap.fetch_pixels(genome_range, genome_range2)
        return pixels

    def infer_binsize(self, genome_range1, genome_range2=None, **kwargs) -> int:
        from coolbox.utilities.hic.wrap import StrawWrap

        path = self.properties['file']
        wrap = StrawWrap(path, normalization=self.balance, binsize=kwargs.get('resolution', 'auto'))
        gr1 = to_gr(genome_range1)
        return wrap.infer_binsize(gr1)
