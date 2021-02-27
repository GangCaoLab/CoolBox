import numpy as np

from coolbox.utilities.genome import GenomeRange
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

    DEFAULT_ARGS_HIC = {
        "transform": "no",
        "normalize": "expect",
        "gaussian_sigma": "no",
        "process_func": "no",
    }

    DEFAULT_PROPERTIES = {
        'style': HicMatBase.STYLE_WINDOW,
        "cmap": "bwr",
        "transform": False,
        'diff_method': "diff",
        'args_hic': DEFAULT_ARGS_HIC,
    }

    def __init__(self, hic1, hic2, **kwargs):
        properties = HiCDiff.DEFAULT_PROPERTIES.copy()
        properties.update(kwargs)
        super().__init__(**properties)

        # TODO how to set hic_args in cli mode ?
        if isinstance(hic1, str):
            hic1 = HiCMat(hic1)
        if isinstance(hic2, str):
            hic2 = HiCMat(hic2)
        self.properties.update({
            'hic1': hic1,
            'hic2': hic2,
        })

        self.mat1 = None
        self.mat2 = None

    def fetch_data(self, gr: GenomeRange, **kwargs) -> np.ndarray:
        # TODO find a better way
        # For cases when updating of Frame's depth_ratio leads to un equal matrix sizes for different tracks).
        args_hic = self.properties.copy()
        args_hic.update(self.properties.get('args_hic', {}))
        for hic in ('hic1', 'hic2'):
            self.properties[hic].properties.update(args_hic)

        hic1: HicMatBase = self.properties['hic1']
        hic2: HicMatBase = self.properties['hic2']
        # must set to avoid re-change GenomeRange in window style
        kwargs['gr_updated'] = True
        # use transformed gr for cases in 'window' style, other wise the gr would be transformed two times
        self.mat1 = mat1 = hic1.fetch_plot_data(gr, **kwargs)
        self.mat2 = mat2 = hic2.fetch_plot_data(gr, **kwargs)
        diff_mat = self.diff_matrix(mat1, mat2)
        try:
            self.SMALL_VALUE = diff_mat[diff_mat > 0].min()
        except Exception:
            pass
        return diff_mat

    def diff_matrix(self, mat1, mat2):
        diff_mth = self.properties['diff_method']
        if diff_mth == 'log2fc':
            return np.log2((mat1 + 1) / (mat2 + 1))
        else:
            return mat1 - mat2
