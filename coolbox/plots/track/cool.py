from collections import OrderedDict

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

import numpy as np
from scipy import ndimage

from coolbox.utilities import (
    GenomeRange,
    change_chrom_names,
    get_logger
)

from coolbox.plots.track.base import TrackPlot


log = get_logger(__name__)


STYLE_TRIANGULAR = 'triangular'
STYLE_MATRIX = 'matrix'
STYLE_WINDOW = 'window'

DEPTH_FULL = 'full'


def is_muliti_cool(cooler_file):
    """
    Judge a cooler is muliti-resolution cool or not.

    Parameters
    ----------
    cooler_file : str
        Path to cooler file.
    """
    import re
    if re.match(".*::/resolutions/[0-9]+$", cooler_file):
        return False

    import h5py
    h5_file = h5py.File(cooler_file)
    return 'resolutions' in h5_file


def get_cooler_resolutions(cooler_file):
    """
    Get the resolutions of a muliti-cooler file

    Parameters
    ----------
    cooler_file : str
        Path to cooler file.
    """
    import h5py
    h5_file = h5py.File(cooler_file)
    resolutions = list(h5_file['resolutions'])
    resolutions = [int(res) for res in resolutions]
    resolutions.sort()
    return resolutions


class PlotCool(TrackPlot):

    DEFAULT_COLOR = 'YlOrRd'

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)

        file_ = self.properties['file']
        import cooler
        self.is_multi = is_muliti_cool(file_)
        if self.is_multi:
            resolutions = get_cooler_resolutions(file_)

            def cooler_reso(resolution):
                return cooler.Cooler(file_ + "::/resolutions/{}".format(resolution))

            self.coolers = OrderedDict([(reso, cooler_reso(reso)) for reso in resolutions])
            self.cool = self.coolers[resolutions[0]] # set self.cool to lowest resolution
        else:
            self.cool = cooler.Cooler(file_)

        self.__set_default_properties()

        self.small_value = 1e-12
        self.ax = None
        self.label_ax = None
        self.matrix = None
        self._out_of_bound = False

    def __set_default_properties(self):
        self.properties['height'] = 'cool_auto'

        if 'color' not in self.properties:
            self.properties['color'] = PlotCool.DEFAULT_COLOR
        if 'style' not in self.properties:
            self.properties['style'] = STYLE_TRIANGULAR
        if 'balance' not in self.properties:
            self.properties['balance'] = 'no'
        if 'color_bar' not in self.properties:
            self.properties['color_bar'] = 'yes'
        if 'transform' not in self.properties:
            self.properties['transform'] = 'no'
        if 'title' not in self.properties:
            self.properties['title'] = ''
        if 'depth_ratio' not in self.properties:
            self.properties['depth_ratio'] = DEPTH_FULL
        if 'norm' not in self.properties:
            self.properties['norm'] = 'log'

    @property
    def is_inverted(self):
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            return True
        else:
            # default: not inverted
            return False

    @property
    def style(self):
        if 'style' in self.properties:
            return self.properties['style']
        else:
            # default triangular style
            return STYLE_TRIANGULAR

    @property
    def is_balance(self):
        if 'balance' in self.properties and self.properties['balance'] == 'no':
            return False
        else:
            # default: balance
            return True

    def __transform_matrix(self, arr):
        if self.properties['transform'] == 'log10':
            arr = np.log10(arr)
        elif self.properties['transform'] == 'log2':
            arr = np.log2(arr)
        elif self.properties['transform'] == 'log':
            arr = np.log(arr)
        return arr

    @property
    def matrix_val_range(self):
        small = 1e-4
        arr = self.matrix
        arr_no_nan = arr[np.logical_not(np.isnan(arr))]

        if self.properties['min_value'] == 'auto':
            min_ = small
        else:
            min_ = self.properties['min_value']

        if self.properties['max_value'] == 'auto':
            max_ = arr_no_nan.max()
        else:
            max_ = self.properties['max_value']

        return min_, max_

    def fetch_matrix(self, genome_range, resolution='auto'):
        """
        Fetch the matrix from the cooler file.

        Parameters
        ----------
        genome_range : coolbox.utilities.GenomeRange
            The genome range to fetch.

        resolution : {'auto', int}
            The matrix resolution, for multi-cooler file.
            Use 'auto' to infer the resolution automatically.
        """
        if self.is_multi:
            resolutions = [k for k in self.coolers.keys()]
            if resolution == 'auto':
                reso = self.__infer_resolution(genome_range, resolutions)
            else:
                assert resolution in resolutions, \
                    "Multi-Cooler file not contain the resolution {}.".format(resolution)
                reso = int(resolution)
            cool = self.coolers[reso]
        else:
            cool = self.cool

        chrom, start, end = genome_range.chrom, genome_range.start, genome_range.end
        if chrom not in cool.chromnames:
            chrom = change_chrom_names(chrom)

        range_str = str(GenomeRange(chrom, start, end))
        arr = cool.matrix(balance=self.is_balance).fetch(range_str)

        # fill zero and nan with small value
        small = self.small_value
        arr[arr == 0] = small
        arr[np.isnan(arr)] = small

        if 'transform' in self.properties and self.properties['transform'] != 'no':
            arr = self.__transform_matrix(arr)

        return arr

    def __infer_resolution(self, genome_range, resolutions):
        BIN_THRESH = 1000

        resolutions.sort()
        reso = resolutions[0]
        for r in resolutions:
            num_bins = genome_range.length // reso
            if num_bins >= BIN_THRESH:
                reso = r
            else:
                break

        return reso

    def __get_triangular_matrix(self, arr):
        small = self.small_value
        tri_matrix = ndimage.rotate(arr, 45, prefilter=False, cval=small)

        rows = tri_matrix.shape[0]

        tri_matrix = tri_matrix[0:(rows//2 + 1), :]

        # cut depth
        if self.properties['depth_ratio'] != 'auto' and self.properties['depth_ratio'] != DEPTH_FULL:
            depth_ratio = float(self.properties['depth_ratio'])
            depth = int(tri_matrix.shape[0] * depth_ratio)
            tri_matrix = tri_matrix[-depth:, :]

        return tri_matrix

    def __get_window_matrix(self, arr):
        small = self.small_value
        window_matrix = ndimage.rotate(arr, 45, prefilter=False, cval=small)
        rows, cols = window_matrix.shape
        if self._out_of_bound == 'left':
            # left side out of bound
            x = cols // 3
            window_matrix = window_matrix[(rows//6):((rows//2) + 1), :(2*x+1)]
        elif self._out_of_bound == 'right':
            # right side out of bound
            x = cols // 3
            window_matrix = window_matrix[(rows//6):((rows//2) + 1), :(2*x+1)]
        elif self._out_of_bound == 'both':
            # double side out of bound
            x = cols // 3
            window_matrix = window_matrix[(rows//6):((rows//2) + 1), :]
        else:
            # normal
            x = cols // 4
            window_matrix = window_matrix[(rows//4):(rows//2 + 1), x:(3*x + 1)]

        # cut depth
        if self.properties['depth_ratio'] != 'auto' and self.properties['depth_ratio'] != DEPTH_FULL:
            depth_ratio = float(self.properties['depth_ratio'])
            depth = int(window_matrix.shape[0] * depth_ratio)
            window_matrix = window_matrix[-depth:, :]

        return window_matrix

    def __plot_matrix(self, genome_range):
        start, end = genome_range.start, genome_range.end
        ax = self.ax
        arr = self.matrix
        cmap = plt.get_cmap(self.properties['color'])
        cmap.set_bad("white")
        cmap.set_under("white")
        c_min, c_max = self.matrix_val_range

        depth_ratio = 1.0 if self.properties['depth_ratio'] == DEPTH_FULL else self.properties['depth_ratio']

        if self.style == STYLE_TRIANGULAR:
            # triangular style
            tri_matrix = self.__get_triangular_matrix(arr)
            img = ax.matshow(tri_matrix, cmap=cmap,
                             extent=(start, end, 0, depth_ratio * (end - start)/2),
                             aspect='auto')
        elif self.style == STYLE_WINDOW:
            # window style
            window_matrix = self.__get_window_matrix(arr)
            img = ax.matshow(window_matrix, cmap=cmap,
                             extent=(start, end, 0, depth_ratio * (end - start)/2),
                             aspect='auto')
        else:
            # matrix style
            img = ax.matshow(arr, cmap=cmap,
                             extent=(start, end, end, start),
                             aspect='auto')

        if self.properties['norm'] == 'log':
            img.set_norm(colors.LogNorm(vmin=c_min, vmax=c_max))
        else:
            img.set_norm(colors.Normalize(vmin=c_min, vmax=c_max))

        return img

    def __adjust_figure(self, genome_range):
        ax = self.ax
        start, end = genome_range.start, genome_range.end
        if self.style == STYLE_TRIANGULAR or self.style == STYLE_WINDOW:

            if self.properties['depth_ratio'] == DEPTH_FULL:
                depth = genome_range.length / 2
            else:
                depth = (genome_range.length / 2) * self.properties['depth_ratio']

            if self.is_inverted:
                ax.set_ylim(depth, 0)
            else:
                ax.set_ylim(0, depth)
        else:
            ax.set_ylim(end, start)
        ax.set_xlim(start, end)

    def __plot_colorbar(self, img):
        ax_divider = make_axes_locatable(self.ax)
        if self.is_inverted:
            cax = ax_divider.append_axes("top", size=0.09, pad=0.2)
        else:
            cax = ax_divider.append_axes("bottom", size=0.09, pad=0.2)
        colorbar(img, cax=cax, orientation='horizontal')

    def __plot_label(self):
        self.label_ax.text(0.15, 0.5, self.properties['title'],
                           horizontalalignment='left', size='large',
                           verticalalignment='center', transform=self.label_ax.transAxes)

    def __fetch_window_matrix(self, genome_range):
        from copy import copy
        fetch_range = copy(genome_range)
        x = (genome_range.end - genome_range.start) // 2
        fetch_range.start = genome_range.start - x
        fetch_range.end = genome_range.end + x

        if fetch_range.start < 0:
            fetch_range.start = genome_range.start
            self._out_of_bound = 'left'

        try:
            arr = self.fetch_matrix(fetch_range)
        except ValueError as e:
            if self._out_of_bound == 'left':
                self._out_of_bound = 'both'
                arr = self.fetch_matrix(genome_range)
            else:
                self._out_of_bound = 'right'
                fetch_range.end = genome_range.end
                arr = self.fetch_matrix(fetch_range)
        return arr, fetch_range

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):

        self._out_of_bound = False

        log.debug("plotting {}".format(self.properties['file']))

        genome_range = GenomeRange(chrom_region, start_region, end_region)

        self.ax = ax
        self.label_ax = label_ax

        # fetch matrix and perform transform process
        if self.style == STYLE_WINDOW:
            arr, fetch_region = self.__fetch_window_matrix(genome_range)
            self.fetch_region = fetch_region
        else:
            arr = self.fetch_matrix(genome_range)

        self.matrix = arr

        # plot matrix
        img = self.__plot_matrix(genome_range)
        self.__adjust_figure(genome_range)

        # plot colorbar
        if self.properties['color_bar'] == 'yes':
            self.__plot_colorbar(img)
        else:
            pass

        # plot label
        self.__plot_label()

    def get_tracks_height(self, frame_width):
        """
        calculate track height dynamically.
        """
        if self.style == STYLE_TRIANGULAR:
            cool_height = frame_width * 0.5
        elif self.style == STYLE_WINDOW:
            if 'height' in self.properties and self.properties['height'] != 'cool_auto':
                cool_height = self.properties['height']
            else:
                cool_height = frame_width * 0.3
        else:
            cool_height = frame_width * 0.8

        if 'depth_ratio' in self.properties and self.properties['depth_ratio'] != DEPTH_FULL:
            if self.properties['style'] != STYLE_MATRIX:
                cool_height = cool_height * self.properties['depth_ratio']

        if 'color_bar' in self.properties and self.properties['color_bar'] != 'no':
            cool_height += 1.5

        return cool_height
