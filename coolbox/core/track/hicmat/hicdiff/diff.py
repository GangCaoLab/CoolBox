from abc import ABC

from coolbox.utilities.genome import GenomeRange
import numpy as np

from coolbox.core.track.hicmat.base import HicMatBase
from coolbox.core.track.hicmat.hicmat import HiCMat


class HiCDiff(HicMatBase):
    """
    Track for express the comparison between two HiC Track.

    Parameters
    ----------
    hic1 : coolbox.api.track.Cool
        First HiC Track or hic file path(.cool, .mcool, .hic).

    hic2 : coolbox.api.track.Cool
        Second HiC Track or hic file path(.cool, .mcool, .hic).

    args_hic: dict
        Key word arguments send to create Cool/DoctHic instance if the input hic1/hic2 is file.

    diff_method : str
        Difference method ('diff', 'log2fc'), default 'diff'


    """

    DEFAULT_PROPERTIES = {
        'style': HicMatBase.STYLE_WINDOW,
        "cmap": "RdYlBu",
        "transform": False,
        "norm": "no",
        'diff_method': "diff",
    }

    def __init__(self, hic1, hic2, args_hic=None, **kwargs):
        # TODO how to set hic_args in cli mode ?
        if isinstance(hic1, str):
            hic1 = HiCMat(hic1, **(args_hic or {}))
        if isinstance(hic2, str):
            hic2 = HiCMat(hic2, **(args_hic or {}))

        properties = HiCDiff.DEFAULT_PROPERTIES.copy()
        properties.update({
            "hic1": hic1,
            "hic2": hic2,
            **kwargs
        })
        super().__init__(**properties)

        for hic in hic1, hic2:  # update related hic track
            hic.properties.update({
                "normalize": self.properties["normalize"],
                "resolution": self.properties['resolution'],
            })

    def fetch_data(self, gr: GenomeRange, **kwargs) -> np.ndarray:
        hic1, hic2 = self.properties['hic1'], self.properties['hic1']
        mat1 = hic1.fetch_data(gr, **kwargs)
        mat2 = hic2.fetch_data(gr, **kwargs)
        diff_mat = self.diff_matrix(mat1, mat2)
        try:
            self.small_value = diff_mat[diff_mat > 0].min()
        except Exception:
            self.small_value = 1e-12
        return diff_mat

    def diff_matrix(self, mat1, mat2):
        diff_mth = self.properties['diff_method']
        if diff_mth == 'log2fc':
            return np.log2((mat1 + 1) / (mat2 + 1))
        else:
            return mat1 - mat2
