import numpy as np

from ..base import Track
from ..hicmat import HiCMat
from ..hicmat.plot import PlotHiCMatrix


class HiCDiff(Track, PlotHiCMatrix):
    """
    Track for express the comparison between two HiC Track.

    Parameters
    ----------
    hic1 : coolbox.api.track.Cool
        First HiC Track or hic file path(.cool, .mcool, .hic).

    hic2 : coolbox.api.track.Cool
        Second HiC Track or hic file path(.cool, .mcool, .hic).

    args_hic : dict, optional
        Argument to create Hi-C instance, only in use
        when first or second argument is a path.

    style : {'triangular', 'window', 'matrix'}, optional
        Matrix style, default 'triangular'.

    depth_ratio : float, optional
        Depth ratio of triangular matrix, use 'full' for full depth. default 'full'.

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    normalize : str
        Normalization method ('none', 'zscore', 'total', 'expect'), default 'expect'

    diff_method : str
        Difference method ('diff', 'log2fc'), default 'diff'

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

    """

    DEFAULT_COLOR = "RdYlBu"

    def __init__(self, hic1, hic2, args_hic=None, **kwargs):
        args_hic = args_hic or {}
        if isinstance(hic1, str):
            hic1 = HiCMat(hic1, **args_hic)
        if isinstance(hic2, str):
            hic2 = HiCMat(hic2, **args_hic)
        properties_dict = {
            "hic1": hic1,
            "hic2": hic2,
            "resolution": "auto",
            "normalize": "expect",
            "diff_method": "diff",
            "style": "triangular",
            "depth_ratio": "full",
            "cmap": HiCDiff.DEFAULT_COLOR,
            "color_bar": "vertical",
            "max_value": "auto",
            "min_value": "auto",
            "title": '',
        }
        properties_dict.update(kwargs)
        for hic in hic1, hic2:  # update related hic track
            hic.properties.update({
                "normalize":  properties_dict["normalize"],
                "resolution": properties_dict['resolution'],
            })
        properties_dict['color'] = properties_dict['cmap']  # change key word

        super().__init__(properties_dict)

        self.properties['transform'] = 'no'
        self.properties['norm'] = 'no'

    def fetch_matrix(self, genome_range, resolution=None):
        diff = self.fetch_data(genome_range, resolution)
        try:
            self.small_value = diff[diff > 0].min()
        except:
            self.small_value = 1e-12
        return diff

    def fetch_related_tracks(self, genome_range, resolution=None):
        if resolution:
            reso = resolution
        else:
            reso = self.properties['resolution']
        hic1 = self.properties['hic1']
        hic2 = self.properties['hic2']
        mat1 = hic1.fetch_matrix(genome_range, resolution=reso)
        mat2 = hic2.fetch_matrix(genome_range, resolution=reso)
        return mat1, mat2

    def __diff_data(self, mat1, mat2):
        diff_mth = self.properties['diff_method']
        if diff_mth == 'log2fc':
            return np.log2((mat1 + 1)/(mat2 + 1))
        else:
            return mat1 - mat2

    def fetch_data(self, genome_range, resolution=None):
        mat1, mat2 = self.fetch_related_tracks(genome_range, resolution)
        diff = self.__diff_data(mat1, mat2)
        return diff
