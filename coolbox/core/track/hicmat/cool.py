import numpy as np

from coolbox.utilities import GenomeRange
from .base import HicMatBase


class Cool(HicMatBase):
    """
    Cool track from .cool, .mcool file.

    Parameters
    ----------
    file: str
        The file path of .cool, .mcool file

    balance: bool
        If use the balanced contact matrix.


    """

    DEFAULT_PROPERTIES = {
        'cmap': "JuiceBoxLike",
        'balance': True,
        "norm": "no",
        "transform": "log"
    }

    def __init__(self, file, **kwargs):
        properties = Cool.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(**properties)

    def fetch_data(self, gr: GenomeRange, **kwargs) -> np.ndarray:
        from coolbox.utilities.hic.wrap import CoolerWrap

        path = self.properties['file']
        wrap = CoolerWrap(path, balance=self.balance, binsize=kwargs.get('resolution', 'auto'))
        arr = wrap.fetch(gr, kwargs.get('gr2'))

        self.fetched_binsize = wrap.fetched_binsize  # expose fetched binsize

        return self.fill_zero_nan(arr)

    def fetch_pixels(self, gr: GenomeRange, **kwargs):
        """
        Fetch the pixels table of upper triangle of the original contact matrix(not processed).

        Parameters
        ----------
        gr2 : GenomeRange, optional.

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

        Returns
        -------
        pixels : pandas.core.frame.DataFrame
            Hi-C pixels table.
            The pixel table contains the non-zero upper triangle entries of the contact map.
        """
        from coolbox.utilities.hic.wrap import CoolerWrap

        path = self.properties['file']
        balance = kwargs.get('balance', self.is_balance)
        wrap = CoolerWrap(path, balance=balance, binsize=kwargs.get('resolution', 'auto'))

        pixels = wrap.fetch_pixels(gr, kwargs.get('gr2'), join=kwargs.get('join', True))
        return pixels

    def infer_binsize(self, gr: GenomeRange, **kwargs) -> int:
        from coolbox.utilities.hic.wrap import CoolerWrap

        path = self.properties['file']
        wrap = CoolerWrap(path, balance=self.balance, binsize=kwargs.get('resolution', 'auto'))
        return wrap.infer_binsize(gr)
