from abc import ABC
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as smm
from scipy.ndimage import gaussian_filter
from scipy.stats import norm

from coolbox.utilities.genome import GenomeRange
from coolbox.core.track.hicmat.base import HicMatBase
from coolbox.core.track.hicmat.hicmat import HiCMat


class Selfish(HicMatBase):
    """
    Differential chromatin interaction.

    Parameters
    ----------
    hic1 : coolbox.api.track.Cool
        First HiC Track or hic file path(.cool, .mcool, .hic).

    hic2 : coolbox.api.track.Cool
        Second HiC Track or hic file path(.cool, .mcool, .hic).

    hic_args : dict, optional
        Argument to create Hi-C instance, only in use
        when first or second argument is a path.

    sigma0 : float, optional
        Initial sigma, parameter of SELFISH method. Default: 1.6

    s : int, optional
        Iteration count parameter of SELFISH method. Default: 10



    Reference
    ---------
    Abbas Roayaei Ardakany, Ferhat Ay, Stefano Lonardi,
    Selfish: discovery of differential chromatin interactions via a self-similarity measure,
    Bioinformatics, Volume 35, Issue 14, July 2019, Pages i145–i153,
    https://doi.org/10.1093/bioinformatics/btz362

    """

    DEFAULT_PROPERTIES = {
        "style": HicMatBase.STYLE_WINDOW,
        "cmap": "RdPu_r",
        "transform": False,
        "norm": "log",
        "s": 10,
        "sigma0": 1.6,
    }

    def __init__(self, hic1, hic2, args_hic=None, **kwargs):
        args_hic_ = {
            "normalize": "zscore",
        }
        if args_hic:
            args_hic_.update(args_hic)
        if isinstance(hic1, str):
            hic1 = HiCMat(hic1, **args_hic_)
        if isinstance(hic2, str):
            hic2 = HiCMat(hic2, **args_hic_)

        properties = Selfish.DEFAULT_PROPERTIES.copy()
        properties.update({
            'hic1': hic1,
            'hic2': hic2,
            **kwargs
        })
        super().__init__(**properties)
        for hic in hic1, hic2:  # update related hic track
            hic.properties.update({
                "resolution": self.properties["resolution"],
            })

        self.zero_indices = None
        self.mat1 = None
        self.mat2 = None

    def fetch_pixels(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        pvals = self.fetch_data(gr, **kwargs)
        mat1 = self.mat1
        mat2 = self.mat2
        binsize = self.properties['hic1'].fetched_binsize
        idx_y, idx_x = np.where(pvals <= kwargs.get('threshold', 1e-4))
        ix_up = idx_y <= idx_x
        idx_y, idx_x = idx_y[ix_up], idx_x[ix_up]
        start = (gr.start // binsize) * binsize
        start1 = start + idx_y * binsize
        end1 = start1 + binsize
        start2 = start + idx_x * binsize
        end2 = start2 + binsize
        df = pd.DataFrame({
            'chrom': gr.chrom,
            'start1': start1,
            'end1': end1,
            'start2': start2,
            'end2': end2,
            'value1': mat1[idx_y, idx_x],
            'value2': mat2[idx_y, idx_x],
            'qvalue': pvals[idx_y, idx_x],
        })
        return df

    def _get_scales(self):
        sigma0 = self.properties['sigma0']
        s = self.properties['s']
        return [(sigma0 * (2 ** (i / s))) for i in range(1, s + 3)]

    def diff_matrix(self, mat1: np.ndarray, mat2: np.ndarray) -> np.ndarray:
        diff = mat2 - mat1
        scales = self._get_scales()
        np.nan_to_num(diff, copy=False, posinf=0, neginf=0, nan=0)
        d_pre = gaussian_filter(diff, scales[0])
        final_p = np.ones(mat1.shape, dtype=mat1.dtype)
        for scale in scales[1:]:
            d_post = gaussian_filter(diff, scale)
            d_diff = d_post - d_pre
            params = norm.fit(d_diff)
            pvals = norm.cdf(d_diff, loc=params[0], scale=params[1])
            np.nan_to_num(pvals, copy=False, posinf=1, neginf=1, nan=1)
            pvals[pvals > 0.5] = 1 - pvals[pvals > 0.5]
            pvals *= 2
            lt_idx = pvals < final_p
            final_p[lt_idx] = pvals[lt_idx]
            d_pre = d_post.copy()
        _, out_p = smm.multipletests(final_p.ravel(), method='fdr_bh')[:2]
        out_p = out_p.reshape(*diff.shape)
        thr = 12
        out_p[out_p == 0] = 10 ** -1 * thr
        return out_p

    def fetch_data(self, gr: GenomeRange, **kwargs) -> np.ndarray:
        hic1, hic2 = self.properties['hic1'], self.properties['hic2']
        kwargs.update({
            'resolution': self.properties['resolution']
        })
        self.mat1 = mat1 = hic1.fetch_data(gr, **kwargs)
        self.mat2 = mat2 = hic2.fetch_data(gr, **kwargs)
        # what does this used for ?
        # self.zero_indices = (hic1.zero_indices | hic1.nan_indices) | (hic2.zero_indices | hic2.nan_indices)
        pvals = self.diff_matrix(mat1, mat2)

        return pvals
