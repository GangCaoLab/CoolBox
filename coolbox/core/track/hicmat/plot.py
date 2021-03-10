import copy
import math

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

from coolbox.utilities import get_logger, GenomeRange

log = get_logger(__name__)

JuiceBoxLikeColor = LinearSegmentedColormap.from_list(
    'interaction', ['#FFFFFF', '#FFDFDF', '#FF7575', '#FF2626', '#F70000'])
JuiceBoxLikeColor.set_bad("white")
JuiceBoxLikeColor.set_under("white")
JuiceBoxLikeColor2 = LinearSegmentedColormap.from_list(
    'interaction', ['#FFFFFF', '#FFDFAF', '#FF7555', '#FF2600', '#F70000'])
JuiceBoxLikeColor2.set_bad("white")
JuiceBoxLikeColor2.set_under("white")

cmaps = {
    "JuiceBoxLike": JuiceBoxLikeColor,
    "JuiceBoxLike2": JuiceBoxLikeColor2,
}


class PlotHiCMat(object):
    JuiceBoxLikeColor = cmaps.get('JuiceBoxLike')
    JuiceBoxLikeColor2 = cmaps.get('JuiceBoxLike2')

    def __init__(self, *args, **kwargs):
        self.ax = None
        self.label_ax = None
        self.matrix = None

    def plot_matrix(self, gr: GenomeRange, gr2: GenomeRange = None):
        gr = GenomeRange(gr)

        def get_cmap(color: str):
            cm = color
            if isinstance(cm, str):
                if cm in cmaps:
                    cmap = cmaps[cm]
                else:
                    cmap = plt.get_cmap(color)
                    cmap = copy.copy(cmap)
                    lowest = cmap(0)
                    cmap.set_bad(lowest)
                    cmap.set_under(lowest)
            else:
                cmap = cm
            return cmap

        cmap = get_cmap(self.properties['color'])
        ax = self.ax
        arr = self.matrix
        c_min, c_max = self.matrix_val_range
        if gr2 is None and self.style == self.STYLE_TRIANGULAR:
            # triangular style
            scale_r = 1 / math.sqrt(2)
            r_len = gr.end - gr.start
            # Rotate image using Affine2D, reference:
            #     https://stackoverflow.com/a/50920567/8500469
            tr = transforms.Affine2D().translate(-gr.start, -gr.start) \
                .rotate_deg_around(0, 0, 45) \
                .scale(scale_r) \
                .translate(gr.start + r_len / 2, -r_len / 2)
            img = ax.matshow(arr, cmap=cmap,
                             transform=tr + ax.transData,
                             extent=(gr.start, gr.end, gr.start, gr.end),
                             aspect='auto')
        elif gr2 is None and self.style == self.STYLE_WINDOW:
            # window style
            # exist in HicMatBase
            fgr = self.fetched_gr
            scale_factor = fgr.length / gr.length
            scale_r = scale_factor / math.sqrt(2)
            length_dialog = gr.length * scale_factor
            delta_x = length_dialog * (gr.start - fgr.start) / fgr.length
            delta_x = length_dialog / 2 - delta_x
            tr = transforms.Affine2D().translate(-gr.start, -gr.start) \
                .rotate_deg_around(0, 0, 45) \
                .scale(scale_r) \
                .translate(gr.start + delta_x, -fgr.length / 2)
            img = ax.matshow(arr, cmap=cmap,
                             transform=tr + ax.transData,
                             extent=(gr.start, gr.end, gr.start, gr.end),
                             aspect='auto')
        else:
            if gr2 is None:
                gr2 = gr
            # matrix style
            img = ax.matshow(arr, cmap=cmap,
                             extent=(gr.start, gr.end, gr2.end, gr2.start),
                             aspect='auto')

        if self.norm == 'log':
            img.set_norm(colors.LogNorm(vmin=c_min, vmax=c_max))
        else:
            img.set_norm(colors.Normalize(vmin=c_min, vmax=c_max))

        return img

    def adjust_figure(self, gr: GenomeRange, gr2: GenomeRange = None):
        ax = self.ax
        if gr2 is None:
            if self.style == self.STYLE_TRIANGULAR or self.style == self.STYLE_WINDOW:

                if self.properties['depth_ratio'] == self.DEPTH_FULL:
                    depth = gr.length / 2
                else:
                    depth = (gr.length / 2) * self.properties['depth_ratio']

                if self.is_inverted:
                    ax.set_ylim(depth, 0)
                else:
                    ax.set_ylim(0, depth)
            else:
                ax.set_ylim(gr.end, gr.start)
            ax.set_xlim(gr.start, gr.end)
        else:
            ax.set_xlim(gr.start, gr.end)
            ax.set_ylim(gr2.end, gr2.start)

    def draw_colorbar(self, img):
        # plot colorbar
        if self.properties['color_bar'] == 'no':
            pass
        else:
            if hasattr(self, 'y_ax') and self.properties['color_bar'] == 'vertical':
                self.plot_colorbar(img, orientation='vertical')
            else:
                self.plot_colorbar(img, orientation='horizontal')

    def plot_colorbar(self, img, orientation='vertical'):
        if orientation == 'horizontal':
            ax_divider = make_axes_locatable(self.ax)
            if self.is_inverted:
                cax = ax_divider.append_axes("top", size=0.09, pad=0.2)
            else:
                cax = ax_divider.append_axes("bottom", size=0.09, pad=0.2)
            plt.colorbar(img, cax=cax, orientation='horizontal')
        else:  # vertical
            y_ax = self.y_ax

            if self.norm == 'log':
                from matplotlib.ticker import LogFormatter
                formatter = LogFormatter(10, labelOnlyBase=False)
                aa = np.array([1, 2, 5])
                c_min, c_max = self.matrix_val_range

                def abs_inc(num):
                    if num != 0:
                        sign = num / abs(num)
                        return int(sign * abs(num + 1))
                    else:
                        return 1

                lower_ = int(np.log10(c_min))
                upper_ = abs_inc(int(np.log10(c_max)))
                tick_values = np.concatenate([aa * 10 ** x for x in range(lower_, upper_)])

                c_bar = plt.colorbar(img, ax=y_ax, ticks=tick_values, format=formatter, fraction=0.98)
            else:
                c_bar = plt.colorbar(img, ax=y_ax, fraction=0.98)

            c_bar.solids.set_edgecolor("face")
            c_bar.ax.tick_params(labelsize='smaller')

            c_bar.ax.yaxis.set_ticks_position('left')

    def get_track_height(self, frame_width):
        """
        calculate track height dynamically.
        """
        if self.style == self.STYLE_TRIANGULAR:
            height = frame_width * 0.5
        elif self.style == self.STYLE_WINDOW:
            if 'height' in self.properties and self.properties['height'] != 'hic_auto':
                height = self.properties['height']
            else:
                height = frame_width * 0.5
        else:
            height = frame_width * 0.8

        if 'depth_ratio' in self.properties and self.properties['depth_ratio'] != self.DEPTH_FULL:
            if self.properties['style'] != self.STYLE_MATRIX:
                height = height * self.properties['depth_ratio']

        if 'color_bar' in self.properties and self.properties['color_bar'] != 'no':
            height += 1.5

        return height

    @property
    def is_inverted(self):
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            return True
        else:
            # default: not inverted
            return False

    @property
    def balance(self):
        if self.properties['balance'] == 'no':
            return False
        else:
            file = self.properties['file']
            from coolbox.utilities.hic.tools import hicmat_filetype
            if hicmat_filetype(file) == '.hic':
                if self.properties['balance'] == 'yes':
                    return 'KR'  # default use KR balance
                else:
                    return self.properties['balance']
            else:
                return True

    @property
    def is_balance(self):
        return bool(self.balance)

    @property
    def style(self):
        if 'style' in self.properties:
            return self.properties['style']
        else:
            # default triangular style
            return self.STYLE_TRIANGULAR

    @property
    def matrix_val_range(self):
        small = 1e-4
        arr = self.matrix
        arr_no_nan = arr[np.isfinite(arr)]
        min_, max_ = 1e-4, 1.0

        try:
            if self.properties['min_value'] == 'auto':
                # set minimal value for color bar
                min_ = arr[arr > arr_no_nan.min()].min()
            else:
                min_ = self.properties['min_value']

            if self.properties['max_value'] == 'auto':
                max_ = arr_no_nan.max()
            else:
                max_ = self.properties['max_value']
            if max_ <= min_:
                max_ = min_ + small
        except ValueError as e:
            log.warning(e)
            return min_, max_

        return min_, max_

    @property
    def resolution(self):
        return self.properties['resolution']

    @property
    def norm(self):
        return self.properties['norm']
