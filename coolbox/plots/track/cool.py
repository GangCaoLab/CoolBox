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


class PlotCool(TrackPlot):

    DEFAULT_COLOR = 'YlOrRd'

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)

        import cooler
        self.cool = cooler.Cooler(self.properties['file'])

        self.__set_default_properties()

        self.small_value = 1e-12
        self.ax = None
        self.label_ax = None
        self.matrix = None

    def __set_default_properties(self):
        self.properties['height'] = 'cool_auto'

        if 'color' not in self.properties:
            self.properties['color'] = PlotCool.DEFAULT_COLOR
        if 'triangular' not in self.properties:
            self.properties['triangular'] = 'yes'
        if 'balance' not in self.properties:
            self.properties['balance'] = 'no'
        if 'color_bar' not in self.properties:
            self.properties['color_bar'] = 'yes'
        if 'transform' not in self.properties:
            self.properties['transform'] = 'no'
        if 'title' not in self.properties:
            self.properties['title'] = ''
        if 'depth_ratio' not in self.properties:
            self.properties['depth_ratio'] = 'full'
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
    def is_triangular(self):
        if 'triangular' in self.properties and self.properties['triangular'] == 'no':
            return False
        else:
            # default: triangular form
            return True

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

    def fetch_matrix(self, genome_range):
        chrom, start, end = genome_range.chrom, genome_range.start, genome_range.end
        if chrom not in self.cool.chromnames:
            chrom = change_chrom_names(chrom)

        range_str = str(GenomeRange(chrom, start, end))
        arr = self.cool.matrix(balance=self.is_balance).fetch(range_str)

        # fill zero and nan with small value
        small = self.small_value
        arr[arr == 0] = small
        arr[np.isnan(arr)] = small

        if 'transform' in self.properties and self.properties['transform'] != 'no':
            arr = self.__transform_matrix(arr)

        return arr

    def __get_triangular_matrix(self, arr):
        small = self.small_value
        tri_matrix = ndimage.rotate(arr, 45, prefilter=False, cval=small)

        rows = tri_matrix.shape[0]

        tri_matrix = tri_matrix[0:(rows//2 + 1), :]

        # cut depth
        if self.properties['depth_ratio'] != 'auto' and self.properties['depth_ratio'] != 'full':
            depth_ratio = float(self.properties['depth_ratio'])
            depth = int(tri_matrix.shape[0] * depth_ratio)
            tri_matrix = tri_matrix[-depth:, :]

        return tri_matrix

    def __plot_matrix(self, genome_range):
        start, end = genome_range.start, genome_range.end
        ax = self.ax
        arr = self.matrix
        cmap = plt.get_cmap(self.properties['color'])
        cmap.set_bad("white")
        cmap.set_under("white")
        c_min, c_max = self.matrix_val_range
        if self.is_triangular:
            tri_matrix = self.__get_triangular_matrix(arr)
            img = ax.matshow(tri_matrix, cmap=cmap,
                             extent=(start, end, 0, (end - start)/2),
                             aspect='auto')
        else:
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
        if self.is_triangular:
            if self.is_inverted:
                ax.set_ylim(genome_range.length / 2, 0)
            else:
                ax.set_ylim(0, genome_range.length / 2)
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

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):

        log.debug("plotting {}".format(self.properties['file']))

        genome_range = GenomeRange(chrom_region, start_region, end_region)

        self.ax = ax
        self.label_ax = label_ax

        # fetch matrix and perform transform process
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
        if self.properties['triangular'] == 'yes':
            cool_height = frame_width * 0.5
        else:
            cool_height = frame_width * 0.8

        if 'depth_ratio' in self.properties and self.properties['depth_ratio'] != 'full':
            cool_height = cool_height * self.properties['depth_ratio']

        if 'color_bar' in self.properties and self.properties['color_bar'] != 'no':
            cool_height += 1.5

        return cool_height
