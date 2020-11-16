import numpy as np

from .hist.plot import CoveragePlot
from coolbox.utilities import (
    get_logger, GenomeRange
)
from .base import Track
from .hicmat import HiCMat


log = get_logger(__name__)


class Virtual4C(Track, CoveragePlot):
    """
    Track for view virtual 4C related to a certain genome position,
    and a HiC Track (include `Cool` and `DotHiC`).

    Parameters
    ----------
    hic_track_or_file : {`Cool`, `DotHiC`}
        related hic track or Hi-C file path.

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

    def __init__(self, hic_track_or_file, genome_position, args_hic=None, **kwargs):
        if isinstance(hic_track_or_file, str):
            args_hic = args_hic or {}
            hic_track = HiCMat(hic_track_or_file, **args_hic)
        else:
            hic_track = hic_track_or_file
        properties_dict = {
            'hic': hic_track,
            'color': Virtual4C.DEFAULT_COLOR,
            'height': Virtual4C.DEFAULT_HEIGHT,
            'genome_position': genome_position,
            'bin_width': 3,
            'max_value': 'auto',
            'min_value': 'auto',
            'show_data_range': True,
            'data_range_style': 'y-axis',
            'style': 'line:1',
            'title': '',
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)
        self.hic = self.properties['hic']
        self.position = GenomeRange(self.properties['genome_position'])
        self.bin_width = self.properties['bin_width']
        self.properties['type'] = self.properties['style']

    def fetch_data(self, genome_range):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        Return
        ------
        mean_arr : numpy.ndarray
        """
        return self.fetch_mean_arr(genome_range)

    def plot(self, ax, chrom_region, start_region, end_region):
        self.ax = ax
        genome_range = GenomeRange(chrom_region, start_region, end_region)
        scores_per_bin = self.fetch_mean_arr(genome_range)
        self.plot_coverage(ax, genome_range, scores_per_bin)
        self.plot_label()

    def fetch_mean_arr(self, genome_range):
        from copy import copy
        bin_width = self.bin_width
        position = self.position
        binsize = self.hic.fetched_binsize
        if binsize is None:
            self.hic.fetch_matrix(genome_range)
            binsize = self.hic.fetched_binsize
        window_range = copy(position)
        offset_ = (bin_width - 1) // 2
        assert offset_ >= 0, "bin width must >= 1"
        window_range.start = window_range.start - offset_ * binsize
        window_range.end = window_range.end + offset_ * binsize
        arr = self.hic.fetch_array(window_range, genome_range)
        mean_arr = np.nanmean(arr, axis=0)
        return mean_arr
