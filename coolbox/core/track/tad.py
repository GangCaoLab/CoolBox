import numpy as np
import pandas as pd

from coolbox.utilities import (
    get_logger
)
from coolbox.utilities.genome import GenomeRange
from .bed.base import BedBase
from .bed.fetch import FetchBed

log = get_logger(__name__)


class PlotTAD(object):
    def __init__(self):
        self.init_colormap()
        self.cache_gr = None
        self.cache_res = None

    def fetch_plot_data(self, gr: GenomeRange, **kwargs):
        if gr == self.cache_gr:
            return self.cache_res
        else:
            self.cache_gr = gr
            self.cache_res = self.fetch_data(gr, **kwargs)
            return self.cache_res

    def get_depth_ratio(self, gr=None):
        if 'depth_ratio' not in self.properties:
            return 1.0
        dr = self.properties['depth_ratio']
        if dr == 'full':
            return 1.0
        if dr == 'auto':
            assert gr is not None
            min_dr = 0.1
            tads = self.fetch_plot_data(gr)
            tads = tads[(tads['start'] >= gr.start) & (tads['end'] <= gr.end)]
            if tads.shape[0] > 0:
                dr = (tads['end'] - tads['start']).max() / (gr.end - gr.start)
                return max(dr, min_dr)
            else:
                return min_dr
        return dr

    def get_track_height(self, frame_width, gr):
        height = frame_width * 0.5
        height *= self.get_depth_ratio(gr)
        return height

    def plot_tads(self, ax, gr: GenomeRange, tads: pd.DataFrame):
        """
        Plots the boundaries as triangles in the given ax.
        """
        self.set_colormap(tads)
        depth = (gr.end - gr.start) / 2
        ymax = (gr.end - gr.start)
        if 'track' in self.__dict__:
            from coolbox.core.track.hicmat import HicMatBase
            assert isinstance(self.track, HicMatBase), f"The parent track should be instance of {HicMatBase}"

            hicmat_tri_style = (HicMatBase.STYLE_WINDOW, HicMatBase.STYLE_TRIANGULAR)
            hicmat_ma_style = (HicMatBase.STYLE_MATRIX,)

            hictrack = self.track
            hicmat_style = hictrack.properties['style']

            # TODO Should we add plotting in BigWig, BedGraph, ABCCompartment, Arcs support?(The original codes supports)
            for region in tads.itertuples():
                if hicmat_style in hicmat_tri_style:
                    self.plot_triangular(ax, gr, region, ymax, depth)
                elif hicmat_style in hicmat_ma_style:
                    self.plot_box(ax, gr, region)
                else:
                    raise ValueError(f"unsupported hicmat style {hicmat_style}")
        else:
            for region in tads.itertuples():
                self.plot_triangular(ax, gr, region, ymax, depth)
            dr = self.get_depth_ratio(gr)
            if self.properties['orientation'] == 'inverted':
                ax.set_ylim(depth * dr, 0)
            else:
                ax.set_ylim(0, depth * dr)
            ax.set_xlim(gr.start, gr.end)

        if len(tads) == 0:
            log.debug("No regions found for Coverage {}.".format(self.properties['name']))

    def plot_triangular(self, ax, gr, region, ymax, depth):
        """
              /\
             /  \
            /    \
        _____________________
           x1 x2 x3
        """

        from matplotlib.patches import Polygon
        x1 = region.start
        x2 = x1 + float(region.end - region.start) / 2
        x3 = region.end
        y1 = 0
        y2 = (region.end - region.start)

        y = (y2 / ymax) * depth

        rgb, edgecolor = self.get_rgb_and_edge_color(region)

        triangle = Polygon(np.array([[x1, y1], [x2, y], [x3, y1]]), closed=True,
                           facecolor=rgb, edgecolor=edgecolor,
                           alpha=self.properties['alpha'],
                           linestyle=self.properties['border_style'],
                           linewidth=self.properties['border_width'])
        ax.add_artist(triangle)
        self.plot_score(ax, gr, region, 'triangular', ymax, depth)

    def plot_box(self, ax, gr, region):
        from matplotlib.patches import Rectangle

        x1 = region.start
        x2 = region.end
        x = y = x1
        w = h = (x2 - x1)

        rgb, edgecolor = self.get_rgb_and_edge_color(region)

        fill = self.properties['border_only'] == 'no'

        rec = Rectangle((x, y), w, h,
                        fill=fill,
                        facecolor=rgb,
                        edgecolor=edgecolor,
                        alpha=self.properties['alpha'],
                        linestyle=self.properties['border_style'],
                        linewidth=self.properties['border_width'])
        ax.add_patch(rec)
        self.plot_score(ax, gr, region, 'box')

    def plot_score(self, ax, gr, region, style, ymax=None, depth=None):
        properties = self.properties

        if properties['show_score'] != 'yes':
            return
        bed = region
        score = bed.score
        if not isinstance(score, (float, int)):
            # score is not number not plot
            return
        region_length = region.end - region.start
        if region_length / gr.length < 0.05:
            # region too small not plot score
            return
        font_size = properties['score_font_size']
        if font_size == 'auto':
            # inference the font size
            from math import log2
            base_size = 18
            s_ = (region_length / gr.length) * 10
            s_ = int(log2(s_))
            font_size = base_size + s_
        ratio = properties['score_height_ratio']
        color = properties['score_font_color']
        if style == 'box':
            x1 = region.start
            x2 = region.end
            w = x2 - x1
            x = x2 - w * ratio
            y = x1 + w * ratio
        else:  # triangular
            x = region.begin + region_length * 0.4
            y = (region_length / ymax) * depth * ratio
        ax.text(x, y, "{0:.3f}".format(score), fontsize=font_size, color=color)


class TAD(BedBase, PlotTAD, FetchBed):
    """
    Tad tack from bed file

    Parameters
    ----------
    border_style: str, optional
        Border style of tad. (Default: 'solid')

    border_width: int, optional
        Border width of tad. (Default: '2.0')

    show_score : bool
        Show bed score or not.
        default False.

    score_font_size : {'auto', int}
        Score text font size.
        default 'auto'

    score_font_color : str
        Score text color.
        default '#000000'

    score_height_ratio : float
        (text tag height) / (TAD height). used for adjust the position of Score text.
        default 0.5

    border_only : bool
        Only show border, default False

    depth_ratio : {float, 'auto', 'full'}
        Depth ratio of triangular, use 'full' for full depth, use 'auto' for calculate depth by current genome_range. default 'auto'.

    orientation : {'normal', 'inverted'}
        Invert y-axis or not, default 'normal'

    """

    DEFAULT_PROPERTIES = {
        'alpha': 0.3,
        'border_style': "--",
        'border_width': 2.0,
        "show_score": False,
        "score_font_size": 'auto',
        "score_font_color": "#000000",
        "score_height_ratio": 0.4,
        "border_only": False,
        "depth_ratio": 'auto',
        "orientation": 'inverted',
    }

    def __init__(self, file, **kwargs):
        properties = TAD.DEFAULT_PROPERTIES.copy()
        properties.update(kwargs)
        super().__init__(file, **properties)
        PlotTAD.__init__(self)

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax
        ov_intervals: pd.DataFrame = self.fetch_plot_data(gr, **kwargs)
        self.plot_tads(ax, gr, ov_intervals)
        self.plot_label()
