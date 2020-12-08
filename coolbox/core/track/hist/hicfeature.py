from abc import ABC
from typing import Union, Callable

from scipy import sparse
import numpy as np

from coolbox.utilities import GenomeRange, get_logger

from ..hicmat import HicMatBase, HiCMat
from .base import HistBase

log = get_logger(__name__)


class HicFeature(HistBase, ABC):
    IGNORE_DIAGS = 3

    def __init__(self, hicmat: Union[str, HicMatBase], args_hic=None, **kwargs):
        self.hicmat = HiCMat(hicmat, **(args_hic or {}))
        properties_dict = {
            'file': self.hicmat.properties.get('file', ""),
            'height': self.DEFAULT_HEIGHT,
            'color': self.DEFAULT_COLOR,
        }
        properties_dict.update(kwargs)
        super().__init__(**properties_dict)


class DiScore(HicFeature):
    """
    Directionnality index track.

    Parameters
    ----------
    hicmat: {str, Cool, DotHic}
        The input hicmat file or HicMatBase object used to calculate di.

    window_size: int, optional
        Width of the diamond region along the matrix used to calculate di. default: 40

    method: {'standard', 'adaptive'}, optional
        Method used for calculating di. default: 'adaptive'

    args_hic: dict, optional
        The args passed to Cool or DotHic to create a HicMatBase object. default: None

    height : float, optional
        Height of track, default DiScore.DEFAULT_HEIGHT


    color : str, optional
        Track color, default DiScore.DEFAULT_COLOR

    style : str, optional
        Track graph type, format {'fill', 'line:`size`', 'points:`size`'},
        example: 'line:2', 'points:0.5'. default: 'fill'

    extra : optional

    show_data_range : bool, optional
        Show_data_range or not, default True.

    data_range_style : {'text', 'y-axis'}, optional
        The style of the data range. default: 'y-axis'

    title : str, optional
        Label text, default ''.

    max_value : {float, 'auto'}, optional
        Max value of track, use 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of track, use 'auto' for specify min value automatically, default 'auto'.

    name : str, optional
        Track's name.
    """

    def __init__(self,
                 hicmat: Union[str, HicMatBase],
                 window_size: int = 40,
                 method="adaptive",
                 args_hic: dict = None,
                 **kwargs):
        properties_dict = {
            "style": 'fill',
        }
        properties_dict.update(kwargs)
        super().__init__(hicmat, args_hic, **properties_dict)
        self.di_scorer = self.di_methods(method)
        self.window_size = window_size

    def fetch_plot_data(self, genome_range: Union[str, GenomeRange]) -> np.ndarray:
        return self.fetch_data(genome_range)

    def fetch_data(self, genome_range: Union[str, GenomeRange]) -> np.ndarray:
        ds_mat: np.ndarray = self.hicmat.fetch_data(genome_range)
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

        return self.di_scorer(contacts_up, contacts_down)

    @classmethod
    def di_methods(cls, method: str) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
        def standard(up, down):
            """Compute directionality index described in:\n
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
            """Compute directionality index described in:\n
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
    hicmat: {str, Cool, DotHic}
        The input hicmat file or HicMatBase object used to calculate insulation score.

    window_size: int, optional
        Width of the diamond region along the matrix used to calculate di. default: 20

    normalize: bool, optional
        Weather to log-nomalize the insulation score array. default: true

    args_hic: dict, optional
        The args passed to Cool or DotHic to create a HicMatBase object. default: None

    height : float, optional
        Height of track, default DiScore.DEFAULT_HEIGHT

    color : str, optional
        Track color, default DiScore.DEFAULT_COLOR

    style : str, optional
        Track graph type, format {'fill', 'line:`size`', 'points:`size`'},
        example: 'line:2', 'points:0.5'. default: 'line'

    extra : optional

    show_data_range : bool, optional
        Show_data_range or not, default True.

    data_range_style : {'text', 'y-axis'}, optional
        The style of the data range. default: 'y-axis'

    title : str, optional
        Label text, default ''.

    max_value : {float, 'auto'}, optional
        Max value of track, use 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of track, use 'auto' for specify min value automatically, default 'auto'.

    name : str, optional
        Track's name.
    """

    def __init__(self,
                 hicmat: Union[str, HicMatBase],
                 window_size: int = 20,
                 normalize: bool = True,
                 args_hic: dict = None,
                 **kwargs):
        properties_dict = {
            "style": "line"
        }
        properties_dict.update(kwargs)
        super().__init__(hicmat, args_hic, **properties_dict)
        self.window_size = window_size
        self.normalize = normalize

    def fetch_plot_data(self, genome_range: Union[str, GenomeRange]) -> np.ndarray:
        return self.fetch_data(genome_range)

    def fetch_data(self, genome_range: Union[str, GenomeRange]) -> np.ndarray:
        mat: np.ndarray = self.hicmat.fetch_data(genome_range)
        mlen = mat.shape[0]
        if mlen < self.window_size * 2:
            return np.zeros(mlen)
        window_size = self.window_size

        insu = np.full(mlen, np.nan, dtype=mat.dtype)
        for row in range(window_size, mlen - window_size):
            sub_mat = mat[row - window_size: row, row + 1: row + window_size + 1]
            insu[row] = np.nanmean(sub_mat)

        if self.normalize:
            insu = np.log2(insu / np.nanmean(insu))
            insu[~np.isfinite(insu)] = 0

        return insu


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

    args_hic : dict, optional
        Argument for create hic track, needed only if first argument is a path.

    bin_width : int, optional
        How many bin used for calculate the mean value.
        default 3

    color : str, optional
        Track color.

    height : int, optional
        Track height

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    max_value : {float, 'auto'}, optional
        Max value of track, use 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of track, use 'auto' for specify min value automatically, default 'auto'.

    show_data_range : bool, optional
        Show_data_range or not, default True.

    data_range_style : {'text', 'y-axis'}
        The style of the data range. default: 'y-axis'

    style : str, optional
        Track graph type, format {'fill', 'line:`size`', 'points:`size`'},
        example: 'line:2', 'points:0.5'. default: 'line:2'

    title : str, optional
        Label text, default ''.

    name : str, optional
        Track's name.

    """

    DEFAULT_COLOR = '#2855d8'

    def __init__(self,
                 hicmat: Union[str, HicMatBase],
                 genome_position: str,
                 args_hic: dict = None,
                 **kwargs):
        properties_dict = {
            'color': self.DEFAULT_COLOR,
            'genome_position': genome_position,
            'bin_width': 3,
            'style': 'line:1',
        }
        properties_dict.update(kwargs)
        super().__init__(hicmat, args_hic, **properties_dict)
        self.position = GenomeRange(self.properties['genome_position'])
        self.bin_width = self.properties['bin_width']

    def fetch_data(self, genome_range: Union[str, GenomeRange]) -> np.ndarray:
        return self.fetch_plot_data(genome_range)

    def fetch_plot_data(self, genome_range: Union[str, GenomeRange]) -> np.ndarray:
        # fetch mean array
        from copy import copy
        bin_width = self.bin_width
        position = self.position
        binsize = self.hicmat.fetched_binsize
        if binsize is None:
            self.hicmat.fetch_data(genome_range)
            binsize = self.hicmat.fetched_binsize
        window_range = copy(position)
        offset_ = (bin_width - 1) // 2
        assert offset_ >= 0, "bin width must >= 1"
        window_range.start = window_range.start - offset_ * binsize
        window_range.end = window_range.end + offset_ * binsize
        arr = self.hicmat.fetch_matrix(window_range, genome_range)
        mean_arr = np.nanmean(arr, axis=0)
        return mean_arr
