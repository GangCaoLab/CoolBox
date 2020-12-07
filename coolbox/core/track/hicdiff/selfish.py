import numpy as np
import pandas as pd
import statsmodels.stats.multitest as smm
from scipy.ndimage import gaussian_filter
from scipy.stats import norm

from coolbox.utilities import to_gr
from ..base import Track
from ..hicmat import HiCMat
from ..hicmat.plot import PlotHiCMatrix


class Selfish(Track, PlotHiCMatrix):
    """
    Differential chromatin interaction.

    Parameters
    ----------
    hic1 : coolbox.api.track.Cool
        First HiC Track or hic file path(.cool, .mcool, .hic).

    hic2 : coolbox.api.track.Cool
        Second HiC Track or hic file path(.cool, .mcool, .hic).

    args_hic : dict, optional
        Argument to create Hi-C instance, only in use
        when first or second argument is a path.

    sigma0 : float, optional
        Initial sigma, parameter of SELFISH method. Default: 1.6

    s : int, optional
        Iteration count parameter of SELFISH method. Default: 10

    style : {'triangular', 'window', 'matrix'}, optional
        Matrix style, default 'triangular'.

    depth_ratio : float, optional
        Depth ratio of triangular matrix, use 'full' for full depth. default 'full'.

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    resolution : int, str
        Resolution of sub two sample. default 'auto'

    cmap : {str, matplotlib.colors.Colormap}, optional
        A diverging colormap, positive color represent the first HiC file,
        and negative represent the second HiC file.

    color_bar : bool, optional
        Show color bar or not.

    max_value : {float, 'auto'}, optional
        Max value of hic matrix, use 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of hic matrix, use 'auto' for specify min value automatically, default 'auto'.

    title : str, optional
        Label text, default ''.

    name : str, optional
        Track's name


    Reference
    ---------
    Abbas Roayaei Ardakany, Ferhat Ay, Stefano Lonardi,
    Selfish: discovery of differential chromatin interactions via a self-similarity measure,
    Bioinformatics, Volume 35, Issue 14, July 2019, Pages i145–i153,
    https://doi.org/10.1093/bioinformatics/btz362

    """

    DEFAULT_COLOR = "RdPu_r"

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

        properties_dict = {
            "hic1": hic1,
            "hic2": hic2,
            "resolution": "auto",
            "sigma0": 1.6,
            "s": 10,
            "style": "window",
            "depth_ratio": "full",
            "cmap": Selfish.DEFAULT_COLOR,
            "color_bar": "vertical",
            "max_value": "auto",
            "min_value": "auto",
            "title": '',
        }
        properties_dict.update(kwargs)
        for hic in hic1, hic2:  # update related hic track
            hic.properties.update({
                "resolution": properties_dict["resolution"],
            })
        properties_dict['color'] = properties_dict['cmap']  # change key word

        super().__init__(properties_dict)

        self.properties['transform'] = 'no'
        self.properties['norm'] = 'log'

        self.zero_indices = None
        self.mat1 = None
        self.mat2 = None

    def fetch_matrix(self, genome_range, resolution=None):
        pval = self.fetch_data(genome_range, resolution)
        return pval

    def fetch_related_tracks(self, genome_range, resolution=None):
        if resolution:
            reso = resolution
        else:
            reso = self.properties['resolution']
        hic1 = self.properties['hic1']
        hic2 = self.properties['hic2']
        mat1 = hic1.fetch_matrix(genome_range, resolution=reso)
        mat2 = hic2.fetch_matrix(genome_range, resolution=reso)
        self.mat1 = mat1
        self.mat2 = mat2
        zero_indices1 = hic1.zero_indices | hic1.nan_indices
        zero_indices2 = hic2.zero_indices | hic2.nan_indices
        self.zero_indices = zero_indices1 | zero_indices2
        return mat1, mat2

    def fetch_pixels(self, genome_range, threshold=1e-4, resolution=None):
        gr = to_gr(genome_range)
        pvals = self.fetch_matrix(genome_range, resolution)
        mat1 = self.mat1
        mat2 = self.mat2
        binsize = self.properties['hic1'].fetched_binsize
        idx_y, idx_x = np.where(pvals <= threshold)
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

    def __get_scales(self):
        sigma0 = self.properties['sigma0']
        s = self.properties['s']
        return [(sigma0 * (2 ** (i / s))) for i in range(1, s + 3)]

    def __selfish(self, a, b):
        diff = b - a
        scales = self.__get_scales()
        np.nan_to_num(diff, copy=False, posinf=0, neginf=0, nan=0)
        d_pre = gaussian_filter(diff, scales[0])
        final_p = np.ones(a.shape, dtype=a.dtype)
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

    def fetch_data(self, genome_range, resolution=None):
        mat1, mat2 = self.fetch_related_tracks(genome_range, resolution)
        p_val = self.__selfish(mat1, mat2)
        return p_val
