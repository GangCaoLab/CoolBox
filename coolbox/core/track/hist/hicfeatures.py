import abc
from abc import ABC
from typing import Union, Callable

from scipy import sparse
import numpy as np

from coolbox.utilities import GenomeRange, get_logger

from ..base import Track
from ..hicmat import HicMatBase
from .plot import CoveragePlot

log = get_logger(__name__)


class HicFeature(Track, CoveragePlot, ABC):
    IGNORE_DIAGS = 3

    def __init__(self, hicmat: HicMatBase, **kwargs):
        properties_dict = {
            'file': hicmat.properties.get('file', ""),
            'height': self.DEFAULT_HEIGHT,
            'color': self.DEFAULT_COLOR,
            'style': 'fill',
            'show_data_range': True,
            'data_range_style': 'y-axis',
            'title': '',
            'max_value': 'auto',
            'min_value': 'auto',
        }
        properties_dict.update(kwargs)
        properties_dict['type'] = properties_dict['style']  # change key word
        super().__init__(properties_dict)
        self.genome_range = None
        self.hicmat = hicmat

    @abc.abstractmethod
    def fetch_data(self, genome_range: Union[str, GenomeRange]):
        pass

    @abc.abstractmethod
    def fetch_plot_data(self, genome_range: Union[str, GenomeRange]):
        pass

    def plot(self, ax, chrom_region, start_region, end_region):
        self.ax = ax
        genome_range = GenomeRange(chrom_region, start_region, end_region)
        self.genome_range = genome_range
        plot_data = self.fetch_plot_data(genome_range)
        if plot_data is not None:
            if isinstance(plot_data, tuple):
                scores_per_bin, x_values = plot_data
            else:
                scores_per_bin, x_values = plot_data, None
            self.plot_coverage(ax, genome_range, scores_per_bin, x_values)
        self.plot_label()


class DiScore(HicFeature):

    def __init__(self, hicmat: HicMatBase, window_size: int = 20, method="adaptive", **kwargs):
        properties_dict = {
        }
        properties_dict.update(kwargs)
        super().__init__(hicmat, **properties_dict)
        self.di_scorer = self.di_methods(method)
        self.window_size = window_size

    def fetch_plot_data(self, genome_range: Union[str, GenomeRange]):
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
            di_array = ((mean_down - mean_up) /
                        np.sqrt((var_up + var_down) / (window_size * (window_size - 1))))

            return di_array

        if method not in locals():
            raise RuntimeError("Only support methods {}".format(list(locals().keys())))

        return locals()[method]


class InsuScore(HicFeature):
    def __init__(self, hicmat: HicMatBase, window_size: int = 20, normalize: bool = True, **kwargs):
        properties_dict = {
            "style": "line"
        }
        properties_dict.update(kwargs)
        super().__init__(hicmat, **properties_dict)
        self.window_size = window_size
        self.normalize = normalize

    def fetch_plot_data(self, genome_range: Union[str, GenomeRange]):
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
