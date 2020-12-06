import numpy as np

from coolbox.utilities import to_gr
from coolbox.utilities.doctool import paste_doc
from .base import hic_doc, HicMatBase


@paste_doc(hic_doc)
class Cool(HicMatBase):
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
            "cmap": Cool.DEFAULT_COLOR,
        }
        properties_dict.update(kwargs)

        super().__init__(file_, **properties_dict)

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

    @paste_doc(hic_doc)
    def fetch_matrix(self, genome_range, genome_range2=None, resolution='auto') -> np.ndarray:
        """
        ${fetch_matrix}
        """
        from coolbox.utilities.hic.wrap import CoolerWrap

        path = self.properties['file']
        wrap = CoolerWrap(path, balance=self.balance, binsize=resolution)

        arr = wrap.fetch(genome_range, genome_range2)

        self.fetched_binsize = wrap.fetched_binsize  # expose fetched binsize

        return self.fill_zero_nan(arr)
