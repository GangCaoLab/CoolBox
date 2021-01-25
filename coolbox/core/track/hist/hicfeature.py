from typing import Union, Callable

from scipy import sparse
import numpy as np

from coolbox.utilities import GenomeRange, get_logger

from ..hicmat import HicMatBase, HiCMat
from .base import HistBase

log = get_logger(__name__)


class HicFeature(HistBase):
    """
    HicFeature base track

    Parameters
    -----------
    hicmat: {str, HicMatBase}
        The input hicmat file or HicMatBase object used to calculate di score.

    args_hic : dict, optional
        Key words arguments used to create Cool/DotHic tack when the input hicmat is a file path

    """
    IGNORE_DIAGS = 3

    def __init__(self, hicmat: Union[str, HicMatBase], args_hic=None, **kwargs):
        self.hicmat = HiCMat(hicmat, **(args_hic or {}))
        properties = {
            'file': self.hicmat.properties.get('file', ""),
        }
        properties.update(kwargs)
        super().__init__(**properties)


class DiScore(HicFeature):
    """
    Directionnality index track.

    Parameters
    ----------
    hicmat: {str, HicMatBase}
        The input hicmat file or HicMatBase object used to calculate di score.

    window_size: int, optional
        Width of the diamond region along the matrix used to calculate di. default: 40

    method: {'standard', 'adaptive'}, optional
        Method used for calculating di. default: 'adaptive'

    args_hic : dict, optional
        Key words arguments used to create Cool/DotHic tack when the input hicmat is a file path
    """

    DEFAULT_PROPERTIES = {
        'style': HistBase.STYLE_FILL
    }

    def __init__(self,
                 hicmat: Union[str, HicMatBase],
                 window_size: int = 40,
                 method="adaptive",
                 args_hic: dict = None,
                 **kwargs):
        properties = DiScore.DEFAULT_PROPERTIES.copy()
        properties.update(kwargs)
        super().__init__(hicmat, args_hic, **properties)
        self.window_size = window_size
        self.method = method

    def fetch_data(self, gr: GenomeRange, **kwargs) -> np.ndarray:
        ds_mat: np.ndarray = self.hicmat.fetch_data(gr)
        mlen = ds_mat.shape[0]
        if mlen < self.window_size * 2:
            return np.zeros(mlen)
        mat = sparse.csr_matrix(ds_mat)
        window_size = self.window_size
        ignore_diags = self.IGNORE_DIAGS

        x, y = mat.nonzero()
        dis = y - x
        # filter for non_nan pixels in valid region
        max_len = ignore_diags + window_size
        available = ((dis >= ignore_diags)
                     & (dis < max_len)
                     & (x >= max_len - 1)
                     & (y <= mlen - max_len))
        x, y = x[available], y[available]
        values = np.array(mat[x, y]).ravel()
        non_nan = ~np.isnan(values)
        x, y, values = x[non_nan], y[non_nan], values[non_nan]
        dis = y - x

        contacts_up = np.zeros((mlen, window_size), dtype=mat.dtype)
        contacts_down = np.zeros((mlen, window_size), dtype=mat.dtype)
        for shift in range(ignore_diags, max_len):
            window_pos = shift - ignore_diags
            mask = dis == shift
            _x, _values = x[mask], values[mask]
            contacts_up[_x + shift, window_pos] = _values
            contacts_down[_x, window_pos] = _values
        contacts_up[:max_len, :] = 0
        contacts_down[:max_len, :] = 0

        return self.di_methods(self.method)(contacts_up, contacts_down)

    @classmethod
    def di_methods(cls, method: str) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
        def standard(up, down):
            """
            Compute directionality index described in:\n
            Jesse R.Dixon 2012. Topological domains in mammalian genomes identified by analysis of chromatin interactions.
            """
            up = up.sum(axis=1)
            down = down.sum(axis=1)
            expected = (up + down) / 2.0
            di_array = (np.sign(down - up) *
                        ((up - expected) ** 2 + (down - expected) ** 2)
                        / expected)

            return di_array

        def adaptive(up, down):
            """
            Compute directionality index described in:\n
            Xiao-Tao Wang 2017.HiTAD: Detecting the structural and functional hierarchies of topologically associating
            domains from chromatin interactions.
            """
            window_size = up.shape[1]
            mean_up = up.mean(axis=1)
            mean_down = down.mean(axis=1)
            var_up = np.square(up - mean_up[:, None]).sum(axis=1)
            var_down = np.square(down - mean_down[:, None]).sum(axis=1)
            denom = np.sqrt((var_up + var_down) / (window_size * (window_size - 1)))
            denom[denom == 0] = 1
            di_array = (mean_down - mean_up) / denom

            return di_array

        if method not in locals():
            raise RuntimeError("Only support methods {}".format(list(locals().keys())))

        return locals()[method]


class InsuScore(HicFeature):
    """
    Insulation score track.

    Parameters
    ----------
    hicmat: {str, HicMatBase}
        The input hicmat file or HicMatBase object used to calculate di score.

    window_size: {int, str}, optional
        Width of the diamond region along the matrix used to calculate di. default: 20
        window_size can also be a str with format like "20-40" representing a range of window_sizes.

    normalize: bool, optional
        Weather to log-nomalize the insulation score array. default: true

    method: {'standard', 'adaptive'}, optional
        Method used for calculating di. default: 'adaptive'

    args_hic : dict, optional
        Key words arguments used to create Cool/DotHic tack when the input hicmat is a file path
    """

    DEFAULT_PROPERTIES = {
        "style": "line",
        "cmap": "coolwarm_r",
    }

    def __init__(self,
                 hicmat: Union[str, HicMatBase],
                 window_size: Union[int, str] = 20,
                 normalize: bool = True,
                 args_hic: dict = None,
                 **kwargs):
        properties = InsuScore.DEFAULT_PROPERTIES.copy()
        properties.update(kwargs)
        super().__init__(hicmat, args_hic, **properties)
        self.window_size = window_size
        self.normalize = normalize

    def fetch_data(self, gr: GenomeRange, **kwargs) -> np.ndarray:
        window_size = self.window_size
        try:
            if isinstance(window_size, int):
                window_size = [window_size]
            elif isinstance(window_size, str):
                st, ed = window_size.split("-")
                window_size = list(range(int(st), int(ed)))
            else:
                raise ValueError
        except Exception:
            raise ValueError("window_size should be like int: 12 or str: '10-50'.")

        insus = []
        # for multiple window sizes
        for ws in window_size:
            mat: np.ndarray = self.hicmat.fetch_data(gr)
            mlen = mat.shape[0]
            if mlen < ws * 2:
                insus.append(np.zeros(mlen))
            else:
                insu = np.full(mlen, np.nan, dtype=mat.dtype)
                for row in range(ws, mlen - ws):
                    sub_mat = mat[row - ws: row, row + 1: row + ws + 1]
                    insu[row] = np.nanmean(sub_mat)

                if self.normalize:
                    insu = np.log2(insu / np.nanmean(insu))
                    insu[~np.isfinite(insu)] = 0
                insus.append(insu)

        return insus[0] if len(insus) == 1 else np.array(insus)[::-1]


class Virtual4C(HicFeature):
    """
    Track for view virtual 4C related to a certain genome position,
    and a HiC Track (include `Cool` and `DotHiC`).

    Parameters
    ----------
    hicmat : {str, `Cool`, `DotHiC`}
        Related hic track or Hi-C file path.

    genome_position : str
        related genome position, like: 'chr1:2000000-2000000'

    bin_width : int, optional
        How many bin used for calculate the mean value.
        default 3

    args_hic : dict, optional
        Argument for create hic track, needed only if first argument is a path.
    """

    DEFAULT_PROPERTIES = {
        "style": HistBase.STYLE_LINE,
        "line_width": 1,
        "color": "#2855d8",
        "bin_width": 3,
    }

    def __init__(self,
                 hicmat: Union[str, HicMatBase],
                 genome_position: str,
                 args_hic: dict = None,
                 **kwargs):
        properties = Virtual4C.DEFAULT_PROPERTIES.copy()
        properties.update({
            "genome_position": genome_position,
            **kwargs,
        })
        super().__init__(hicmat, args_hic, **properties)
        self.position = GenomeRange(self.properties['genome_position'])
        self.bin_width = self.properties['bin_width']

    def fetch_data(self, gr: GenomeRange, **kwargs) -> np.ndarray:
        # fetch mean array
        from copy import copy
        bin_width = self.bin_width
        position = self.position
        binsize = self.hicmat.infer_binsize(gr)
        window_range = copy(position)
        offset_ = (bin_width - 1) // 2
        assert offset_ >= 0, "bin width must >= 1"
        window_range.start = window_range.start - offset_ * binsize
        window_range.end = window_range.end + offset_ * binsize
        arr = self.hicmat.fetch_data(window_range, gr2=gr)
        mean_arr = np.nanmean(arr, axis=0)
        return mean_arr
