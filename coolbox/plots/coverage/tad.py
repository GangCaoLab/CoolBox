import numpy as np
from intervaltree import IntervalTree, Interval

from coolbox.utilities import (change_chrom_names, get_logger,
                               opener, ReadBed, GenomeRange)
from coolbox.plots.coverage.base import CoveragePlot


log = get_logger(__name__)


class PlotTADCoverage(CoveragePlot):

    DEFAULT_COLOR = "#1f78b4"
    DEFAULT_ALPHA = 0.4
    DEFAULT_BORDER_COLOR = "#000000"
    DEFAULT_BORDER_WIDTH = 1
    DEFAULT_BORDER_STYLE = "-"

    def __init__(self, *args, **kwargs):
        CoveragePlot.__init__(self, *args, **kwargs)
        self.track = None
        self.bed_type = None  # once the bed file is read, this is bed3, bed6 or bed12
        self.interval_tree = {}  # interval tree of the bed regions

        if 'color' not in self.properties:
            self.properties['color'] = 'bed_rgb'
        if 'show_score' not in self.properties:
            self.properties['show_score'] = 'yes'
        if 'score_font_size' not in self.properties:
            self.properties['score_font_size'] = 'auto'
        if 'score_height_ratio' not in self.properties:
            self.properties['score_height_ratio'] = 0.4
        if 'alpha' not in self.properties:
            self.properties['alpha'] = PlotTADCoverage.DEFAULT_ALPHA
        if 'border_color' not in self.properties:
            self.properties['border_color'] = PlotTADCoverage.DEFAULT_BORDER_COLOR
        if 'border_width' not in self.properties:
            self.properties['border_width'] = PlotTADCoverage.DEFAULT_BORDER_COLOR
        if 'border_style' not in self.properties:
            self.properties['border_style'] = PlotTADCoverage.DEFAULT_BORDER_STYLE
        if 'fill' not in self.properties:
            self.properties['fill'] = 'yes'

        self.interval_tree, min_score, max_score = self.__process_bed()

    def check_track_type(self):
        from coolbox.core.track import BigWig, BedGraph, Cool, HicCompare, Arcs
        CoveragePlot.check_track_type(allow=[BigWig, BedGraph, Cool, HicCompare, Arcs])

    @property
    def track_type(self):
        """
        return a string for indicate self.track's type.
        """
        type_ = self.track.__class__.__name__
        if type_ == 'Cool' or type_ == 'DotHiC':
            type_ = "HiC" + ":" + self.track.properties['style']
        return type_

    def __process_bed(self):

        bed_file_h = ReadBed(opener(self.properties['file']))
        self.bed_type = bed_file_h.file_type

        if 'color' in self.properties and self.properties['color'] == 'bed_rgb' and \
                self.bed_type not in ['bed12', 'bed9']:
            log.warning("*WARNING* Color set to 'bed_rgb', but bed file does not have the rgb field. The color has "
                        "been set to {}".format(PlotTADCoverage.DEFAULT_COLOR))
            self.properties['color'] = PlotTADCoverage.DEFAULT_COLOR

        valid_intervals = 0
        interval_tree = {}

        max_score = float('-inf')
        min_score = float('inf')
        for bed in bed_file_h:
            if bed.score < min_score:
                min_score = bed.score
            if bed.score > max_score:
                max_score = bed.score

            interval_tree.setdefault(bed.chromosome, IntervalTree())
            interval_tree[bed.chromosome].add(Interval(bed.start, bed.end, bed))
            valid_intervals += 1

        if valid_intervals == 0:
            log.warning("No valid intervals were found in file {}".format(self.properties['file_name']))

        return interval_tree, min_score, max_score

    def plot(self, ax, chrom_region, start_region, end_region):
        """
        Plots the boundaries as triangles in the given ax.
        """
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

        fill = True if self.properties['fill'] == 'yes' else False

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

        if self.properties['color'] == 'bed_rgb':
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
                    rgb = PlotTADCoverage.DEFAULT_COLOR
            else:
                rgb = PlotTADCoverage.DEFAULT_COLOR
        return rgb, edgecolor
