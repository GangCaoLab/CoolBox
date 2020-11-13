import numpy as np

from coolbox.utilities import (change_chrom_names, get_logger,
                               GenomeRange)
from .base import Coverage
from coolbox.mixin.bed import FetchBed
from coolbox.utilities.bed import build_bed_index

log = get_logger(__name__)


class TADCoverage(Coverage, FetchBed):
    """
    TAD Coverage is used for plot TAD on upper layer of a track.

    Parameters
    ----------
    file_ : str
        Path to the loop file, loop file is a tab splited text file, fields:
        chr1, x1, x2, chr2, y1, y2, [color], ... (other optional fields)

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

    color : str, optional
        Peak color, use 'bed_rgb' for specify color in file,
        default 'bed_rgb'.

    alpha : float, optional
        Peak alpha value, default 0.3.

    border_color : str, optional
        Border line color, default '#000000'.

    border_width : float, optional
        Border line width, default 1.0.

    border_style : str, optional
        Border line style, default 'solid'.

    """
    DEFAULT_COLOR = "#1f78b4"

    def __init__(self, file_, **kwargs):
        properties_dict = {
            "file": file_,
            "show_score": False,
            "score_font_size": 'auto',
            "score_font_color": "#000000",
            "score_height_ratio": 0.4,
            "border_only": False,
            "color": "bed_rgb",
            "alpha": 0.3,
            "border_color": "#000000",
            "border_width": 2.0,
            "border_style": "--",
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)
        self.track = None
        self.bed_type = None  # once the bed file is read, this is bed3, bed6 or bed12
        self.interval_tree = {}  # interval tree of the bed regions
        self.bgz_file = build_bed_index(file_)

    def check_track_type(self):
        from coolbox.core.track import BigWig, BedGraph, Cool, DotHiC, Arcs
        Coverage.check_track_type(allow=[BigWig, BedGraph, Cool, DotHiC, Arcs])

    @property
    def track_type(self):
        """
        return a string for indicate self.track's type.
        """
        type_ = self.track.__class__.__name__
        if type_ == 'Cool' or type_ == 'DotHiC':
            type_ = "HiC" + ":" + self.track.properties['style']
        return type_

    def plot(self, ax, chrom_region, start_region, end_region):
        """
        Plots the boundaries as triangles in the given ax.
        """
        self.load_range(f"{chrom_region}:{start_region}-{end_region}")
        self.ax = ax
        genome_range = GenomeRange(chrom_region, start_region, end_region)
        self._genome_range = genome_range

        if chrom_region not in self.interval_tree:
            orig = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            log.debug('Chromosome name: {} does not exists. Changing name to {}'.format(orig, chrom_region))

        current_regions = sorted(self.interval_tree[chrom_region][start_region:end_region])
        ymax = max([region.end - region.begin for region in current_regions])
        for region in current_regions:
            if self.track_type.startswith('HiC'):

                if self.track_type == 'HiC:window' or self.track_type == 'HiC:triangular':
                    depth = (end_region - start_region) / 2
                    ymax = (end_region - start_region)
                    self.__plot_triangular(region, ymax, depth)
                else:
                    self.__plot_box(region)

            elif self.track_type == 'HicCompare':
                self.__plot_box(region)

            elif self.track_type in ['BigWig', 'BedGraph', 'ABCompartment', 'Arcs']:
                depth_neg, depth_pos = ax.get_ylim()
                if ('orientation' in self.track.properties) and (self.track.properties['orientation'] == 'inverted'):
                    depth = depth_neg
                else:
                    depth = depth_pos
                self.__plot_triangular(region, ymax, depth)

        if len(current_regions) == 0:
            log.warning("No regions found for Coverage {}.".format(self.properties['name']))

    def __plot_triangular(self, region, ymax, depth):
        """
              /\
             /  \
            /    \
        _____________________
           x1 x2 x3
        """
        ax = self.ax

        from matplotlib.patches import Polygon
        x1 = region.begin
        x2 = x1 + float(region.end - region.begin) / 2
        x3 = region.end
        y1 = 0
        y2 = (region.end - region.begin)

        y = (y2 / ymax) * depth

        rgb, edgecolor = self.__get_rgb_and_edge_color(region.data)

        triangle = Polygon(np.array([[x1, y1], [x2, y], [x3, y1]]), closed=True,
                           facecolor=rgb, edgecolor=edgecolor,
                           alpha=self.properties['alpha'],
                           linestyle=self.properties['border_style'],
                           linewidth=self.properties['border_width'])
        ax.add_artist(triangle)
        self.__plot_score(region, 'triangular', ymax, depth)

    def __plot_box(self, region):
        ax = self.ax
        from matplotlib.patches import Rectangle

        x1 = region.begin
        x2 = region.end
        x = y = x1
        w = h = (x2 - x1)

        rgb, edgecolor = self.__get_rgb_and_edge_color(region.data)

        fill = True if self.properties['border_only'] == 'no' else False

        rec = Rectangle((x, y), w, h,
                        fill=fill,
                        facecolor=rgb,
                        edgecolor=edgecolor,
                        alpha=self.properties['alpha'],
                        linestyle=self.properties['border_style'],
                        linewidth=self.properties['border_width'])
        ax.add_patch(rec)
        self.__plot_score(region, 'box')

    def __plot_score(self, region, style, ymax=None, depth=None):
        ax = self.ax
        genome_range = self._genome_range
        if self.properties['show_score'] != 'yes':
            return
        bed = region.data
        score = bed.score
        if not (isinstance(score, float) or isinstance(score, int)):
            # score is not number not plot
            return
        region_length = region.end - region.begin
        if region_length / genome_range.length < 0.05:
            # region too small not plot score
            return
        font_size = self.properties.get('score_font_size')
        if font_size == 'auto':
            # inference the font size
            from math import log2
            base_size = 18
            s_ = (region_length / genome_range.length) * 10
            s_ = int(log2(s_))
            font_size = base_size + s_
        ratio = self.properties['score_height_ratio']
        color = self.properties['score_font_color']
        if style == 'box':
            x1 = region.begin; x2 = region.end
            w = x2 - x1
            x = x2 - w * ratio
            y = x1 + w * ratio
        else:  # triangular
            x = region.begin + region_length * 0.4
            y = (region_length / ymax) * depth * ratio
        ax.text(x, y, "{0:.3f}".format(score), fontsize=font_size, color=color)

    def __get_rgb_and_edge_color(self, bed):
        rgb = self.properties['color']
        edgecolor = self.properties['border_color']

        if self.properties['border_only'] == "yes":
            rgb = 'none'
        elif self.properties['color'] == 'bed_rgb':
            # if rgb is set in the bed line, this overrides the previously
            # defined colormap
            if self.bed_type in ['bed9', 'bed12'] and len(bed.rgb) == 3:
                try:
                    rgb = [float(x) / 255 for x in bed.rgb]
                    if 'border_color' in self.properties:
                        edgecolor = self.properties['border_color']
                    else:
                        edgecolor = self.properties['color']
                except IndexError:
                    rgb = TADCoverage.DEFAULT_COLOR
            else:
                rgb = TADCoverage.DEFAULT_COLOR
        return rgb, edgecolor
