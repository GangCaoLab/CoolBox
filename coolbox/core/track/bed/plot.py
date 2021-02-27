from typing import Union
import matplotlib
import matplotlib.colors as colors
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from coolbox.utilities import (
    get_logger
)
from coolbox.utilities.genome import GenomeRange

log = get_logger(__name__)


class PlotBed(object):
    def init_for_plot(self):
        from matplotlib import font_manager

        properties = self.properties
        self.len_w = None  # this is the length of the letter 'w' given the font size
        self.counter = None
        self.small_relative = None
        self.is_draw_labels = properties['labels'] == 'on'
        self.fp = font_manager.FontProperties(size=properties['fontsize'])
        # set the distance between rows
        self.row_scale = properties['interval_height'] * 2.3
        # set colormap if possible, else None, wait for setting when plot
        self.colormap = None
        if not matplotlib.colors.is_color_like(self.properties['color']) and self.properties['color'] != 'bed_rgb':
            # check if the color is a valid colormap name
            if self.properties['color'] not in matplotlib.cm.datad:
                log.warning("*WARNING* color: '{}' for Track {} is not valid. Color has "
                            "been set to {}".format(self.properties['color'], self.properties['name'],
                                                    self.COLOR))
                self.properties['color'] = self.COLOR
            else:
                self.colormap = self.properties['color']

    def plot_genes(self, ax, gr: GenomeRange, ov_genes: pd.DataFrame):
        properties = self.properties
        # bed_type
        self.properties['bed_type'] = properties['bed_type'] or self.infer_bed_type(ov_genes)
        # as min_score and max_score change every plot, we compute them for every plot
        min_score, max_score = properties['min_score'], properties['max_score']
        has_score_col = properties['bed_type'] in ('bed6', 'bed9', 'bed12')
        if has_score_col and len(ov_genes):
            min_score = (min_score != 'inf') or ov_genes['score'].min()
            max_score = (max_score != '-inf') or ov_genes['score'].max()
        min_score, max_score = float(min_score), float(max_score)

        # set colormap
        if self.colormap is not None:
            norm = matplotlib.colors.Normalize(vmin=min_score, vmax=max_score)
            cmap = matplotlib.cm.get_cmap(properties['color'])
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        if properties['color'] == 'bed_rgb' and properties['bed_type'] not in ['bed12', 'bed9']:
            log.warning("*WARNING* Color set to 'bed_rgb', but bed file does not have the rgb field. The color has "
                        "been set to {}".format(self.COLOR))
            self.properties['color'] = self.COLOR
            self.colormap = None


        self.counter = 0
        self.small_relative = 0.004 * (gr.end - gr.start)
        self.get_length_w(ax.get_figure().get_figwidth(), gr.start, gr.end)
        # turn labels off when too many intervals are visible.
        if properties['labels'] == 'on' and len(ov_genes) > 60:
            self.is_draw_labels = False

        num_rows = properties['num_rows']
        max_num_row_local = 1
        max_ypos = 0
        # check for the number of other intervals that overlap
        #    with the given interval
        #            1         2
        #  012345678901234567890123456
        #  1=========       4=========
        #       2=========
        #         3============
        #
        # for 1 row_last_position = [9]
        # for 2 row_last_position = [9, 14]
        # for 3 row_last_position = [9, 14, 19]
        # for 4 row_last_position = [26, 14, 19]

        row_last_position = []  # each entry in this list contains the end position
        # of genomic interval. The list index is the row
        # in which the genomic interval was plotted.
        # Any new genomic interval that wants to be plotted,
        # knows the row to use by finding the list index that
        # is larger than its start

        # check for overlapping genes including
        # label size (if plotted)

        for bed in ov_genes.itertuples():
            """
            BED12 gene format with exon locations at the end
            chrX    20850   23076   CG17636-RA      0       -       20850   23017   0       3       946,765,64,     0,1031,2162,

            BED9
            bed with rgb at end
            chr2L   0       70000   ID_5    0.26864549832   .       0       70000   51,160,44

            BED6
            bed without rgb
            chr2L   0       70000   ID_5    0.26864549832   .

            BED3
            bed with only intervals
            chr2L  0        70000
            """
            self.counter += 1

            if self.is_draw_labels:
                num_name_characters = len(bed.name) + 2  # +2 to account for an space before and after the name
                bed_extended_end = int(bed.end + (num_name_characters * self.len_w))
            else:
                bed_extended_end = (bed.end + 2 * self.small_relative)

            # get smallest free row
            if len(row_last_position) == 0:
                free_row = 0
                row_last_position.append(bed_extended_end)
            else:
                # get list of rows that are less than bed.start, then take the min
                idx_list = [idx for idx, value in enumerate(row_last_position) if value < bed.start]
                if len(idx_list):
                    free_row = min(idx_list)
                    row_last_position[free_row] = bed_extended_end
                else:
                    free_row = len(row_last_position)
                    row_last_position.append(bed_extended_end)

            rgb, edgecolor = self.get_rgb_and_edge_color(bed)

            ypos = self.get_y_pos(free_row)

            # do not plot if the maximum interval rows to plot is reached
            if num_rows and free_row >= float(num_rows):
                continue

            if free_row > max_num_row_local:
                max_num_row_local = free_row
            if ypos > max_ypos:
                max_ypos = ypos

            if properties['bed_type'] == 'bed12':
                if properties['gene_style'] == 'flybase':
                    self.draw_gene_with_introns_flybase_style(ax, bed, ypos, rgb, edgecolor)
                else:
                    self.draw_gene_with_introns(ax, bed, ypos, rgb, edgecolor)
            else:
                self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)

            if self.is_draw_labels and bed.start > gr.start and bed.end < gr.end:
                ax.text(bed.end + self.small_relative,
                        ypos + (float(properties['interval_height']) / 2),
                        bed.name,
                        horizontalalignment='left',
                        verticalalignment='center',
                        fontproperties=self.fp)

        if self.counter == 0:
            log.warning(f"*Warning* No intervals were found for file {properties['file']} "
                        f"in Track \'{properties['name']}\' for the interval plotted ({gr}).\n")

        ymax = 0
        if num_rows:
            ymin = float(num_rows) * self.row_scale
        else:
            ymin = max_ypos + properties['interval_height']

        log.debug("ylim {},{}".format(ymin, ymax))
        # the axis is inverted (thus, ymax < ymin)
        ax.set_ylim(ymin, ymax)

        if properties['display'] == 'domain':
            ax.set_ylim(-5, 205)
        elif properties['display'] == 'collapsed':
            ax.set_ylim(-5, 105)

        ax.set_xlim(gr.start, gr.end)

    def draw_gene_with_introns(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws a gene like in flybase gbrowse.
        """
        from matplotlib.patches import Polygon
        properties = self.properties
        height = float(properties['interval_height'])

        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
            return
        half_height = height / 2
        quarter_height = height / 4
        three_quarter_height = quarter_height * 3

        # draw 'backbone', a line from the start until the end of the gene
        ax.plot([bed.start, bed.end], [ypos + half_height, ypos + half_height], 'black', linewidth=0.5, zorder=-1)

        for idx in range(0, bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if x1 < bed.thick_start or x0 > bed.thick_end:
                y0 = ypos + quarter_height
                y1 = ypos + three_quarter_height
            else:
                y0 = ypos
                y1 = ypos + height

            if x0 < bed.thick_start < x1:
                vertices = ([(x0, ypos + quarter_height), (x0, ypos + three_quarter_height),
                             (bed.thick_start, ypos + three_quarter_height),
                             (bed.thick_start, ypos + height),
                             (bed.thick_start, ypos + height),
                             (x1, ypos + height), (x1, ypos),
                             (bed.thick_start, ypos), (bed.thick_start, ypos + quarter_height)])

            elif x0 < bed.thick_end < x1:
                vertices = ([(x0, ypos),
                             (x0, ypos + height),
                             (bed.thick_end, ypos + height),
                             (bed.thick_end, ypos + three_quarter_height),
                             (x1, ypos + three_quarter_height),
                             (x1, ypos + quarter_height),
                             (bed.thick_end, ypos + quarter_height),
                             (bed.thick_end, ypos)])
            else:
                vertices = ([(x0, y0), (x0, y1), (x1, y1), (x1, y0)])

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 linewidth=0.1,
                                 edgecolor='none',
                                 facecolor=rgb))

            if idx < bed.block_count - 1:
                # plot small arrows using the character '<' or '>' over the back bone
                intron_length = bed.block_starts[idx + 1] - (bed.block_starts[idx] + bed.block_sizes[idx])
                marker = 5 if bed.strand == '+' else 4
                if intron_length > 3 * self.small_relative:
                    pos = np.arange(x1 + 1 * self.small_relative,
                                    x1 + intron_length + self.small_relative, int(2 * self.small_relative))
                    ax.plot(pos, np.zeros(len(pos)) + ypos + half_height, '.', marker=marker,
                            fillstyle='none', color='blue', markersize=3)

                elif intron_length > self.small_relative:
                    intron_center = x1 + int(intron_length) / 2
                    ax.plot([intron_center], [ypos + half_height], '.', marker=5,
                            fillstyle='none', color='blue', markersize=3)

    def draw_gene_with_introns_flybase_style(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws a gene using different styles
        """
        from matplotlib.patches import Polygon
        properties = self.properties
        height = float(properties['interval_height'])

        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
            return
        half_height = height / 2
        # draw 'backbone', a line from the start until the end of the gene
        ax.plot([bed.start, bed.end], [ypos + half_height, ypos + half_height], 'black', linewidth=0.5, zorder=-1)

        # get start, end of all the blocks
        positions = []
        for idx in range(0, bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if x0 < bed.thick_start < x1:
                positions.append((x0, bed.thick_start, 'UTR'))
                positions.append((bed.thick_start, x1, 'coding'))

            elif x0 < bed.thick_end < x1:
                positions.append((x0, bed.thick_end, 'coding'))
                positions.append((bed.thick_end, x1, 'UTR'))

            else:
                if x1 < bed.thick_start or x0 > bed.thick_end:
                    seg_type = 'UTR'
                else:
                    seg_type = 'coding'

                positions.append((x0, x1, seg_type))

        # plot all blocks as rectangles except the last if the strand is + or
        # the first is the strand is -, which are drawn as arrows.
        if bed.strand == '-':
            positions = positions[::-1]

        first_pos = positions.pop()
        if first_pos[2] == 'UTR':
            _rgb = 'grey'
        else:
            _rgb = rgb

        vertices = self.draw_arrow(ax, first_pos[0], first_pos[1], bed.strand, ypos)

        ax.add_patch(Polygon(vertices, closed=True, fill=True,
                             edgecolor=edgecolor,
                             facecolor=_rgb,
                             linewidth=0.5))

        for start_pos, end_pos, _type in positions:
            if _type == 'UTR':
                _rgb = 'grey'
            else:
                _rgb = rgb
            vertices = [(start_pos, ypos), (start_pos, ypos + height),
                        (end_pos, ypos + height), (end_pos, ypos)]

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=_rgb,
                                 linewidth=0.5))

    def draw_gene_simple(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws an interval with direction (if given)
        """
        from matplotlib.patches import Polygon
        properties = self.properties
        height = float(properties['interval_height'])

        if bed.strand not in ['+', '-']:
            ax.add_patch(Rectangle((bed.start, ypos), bed.end - bed.start, height,
                                   edgecolor=edgecolor, facecolor=rgb, linewidth=0.5))
        else:
            vertices = self.draw_arrow(ax, bed.start, bed.end, bed.strand, ypos)
            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=rgb,
                                 linewidth=0.5))

    def draw_arrow(self, ax, start, end, strand, ypos):
        properties = self.properties
        height = float(properties['interval_height'])

        half_height = height / 2
        if strand == '+':
            x0 = start
            x1 = end  # - self.small_relative
            y0 = ypos
            y1 = ypos + height

            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

            -----------------\
            ---------------- /

            """

            vertices = [(x0, y0), (x0, y1), (x1, y1), (x1 + self.small_relative, y0 + half_height), (x1, y0)]

        else:
            x0 = start  # + self.small_relative
            x1 = end
            y0 = ypos
            y1 = ypos + height

            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

            /-----------------
            \\-----------------
            """

            vertices = [(x0, y0), (x0 - self.small_relative, y0 + half_height), (x0, y1), (x1, y1), (x1, y0)]

        return vertices

    def get_length_w(self, fig_width, region_start, region_end):
        """
        to improve the visualization of the genes it is good to have an estimation of the label
        length. In the following code I try to get the length of a 'W' in base pairs.
        """
        if self.is_draw_labels:
            # from http://scipy-cookbook.readthedocs.org/items/Matplotlib_LaTeX_Examples.html
            inches_per_pt = 1.0 / 72.27
            font_in_inches = self.properties['fontsize'] * inches_per_pt
            region_len = region_end - region_start
            bp_per_inch = region_len / fig_width
            font_in_bp = font_in_inches * bp_per_inch
            self.len_w = font_in_bp
            log.debug("len of w set to: {} bp".format(self.len_w))
        else:
            self.len_w = 1

        return self.len_w

    def get_y_pos(self, free_row):
        """
        The y_pos is set such that regions to be plotted do not overlap (stacked). To override this
        the properties['collapsed'] needs to be set.

        The algorithm uses a interval tree (self.region_interval) to check the overlaps
        and a sort of coverage vector 'rows used' to identify the row in which to plot

        Return
        ------
        ypos : int
            y position.
        """

        # if the domain directive is given, ypos simply oscilates between 0 and 100
        if self.properties['display'] == 'interlaced':
            ypos = self.properties['interval_height'] if self.counter % 2 == 0 else 1
        elif self.properties['display'] == 'collapsed':
            ypos = 0
        else:
            ypos = free_row * self.row_scale
        return ypos

    def get_rgb_and_edge_color(self, bed):
        # TODO need simplification
        rgb = self.properties['color']
        edgecolor = self.properties['border_color']

        if self.colormap:
            # translate value field (in the example above is 0 or 0.2686...) into a color
            rgb = self.colormap.to_rgba(bed.score)

        # for tad coverage
        if self.properties['style'] == 'tad' and self.properties['border_only'] == 'yes':
            rgb = 'none'
        elif self.properties['color'] == 'bed_rgb':
            # if rgb is set in the bed line, this overrides the previously
            # defined colormap
            if self.properties['bed_type'] in ['bed9', 'bed12'] and len(bed.rgb) == 3:
                try:
                    rgb = [float(x) / 255 for x in bed.rgb]
                    if 'border_color' in self.properties:
                        edgecolor = self.properties['border_color']
                    else:
                        edgecolor = self.properties['color']
                except IndexError:
                    rgb = self.COLOR
            else:
                rgb = self.COLOR
        return rgb, edgecolor

    @staticmethod
    def infer_bed_type(df: pd.DataFrame) -> Union[str, None]:
        #  bed_type of dataframe are store in dataframe's __dict__ in FetchBed.fetch_intervals
        if 'bed_type' in df.__dict__:
            bed_type = df.bed_type
        else:
            bed_types = {
                12: 'bed12',
                9: 'bed9',
                6: 'bed6',
                3: 'bed3'
            }
            num_col = len(df.columns)
            bed_type = bed_types[num_col] if num_col in bed_types else 'bed3'
            if bed_type == 'bed3' and num_col < 3:
                raise ValueError(f"Invalid dataframe for bed3 with columns: {df.columns}")
        return bed_type

    def plot_tads(self, ax, gr: GenomeRange, tads: pd.DataFrame):
        """
        Plots the boundaries as triangles in the given ax.
        """
        from coolbox.core.track.hicmat import HicMatBase
        # coverage only
        assert 'track' in self.__dict__ and isinstance(self.track, HicMatBase), \
            f"The parent track should be instance of {HicMatBase}"

        hicmat_tri_style = (HicMatBase.STYLE_WINDOW, HicMatBase.STYLE_TRIANGULAR)
        hicmat_ma_style = (HicMatBase.STYLE_MATRIX,)

        hictrack = self.track
        hicmat_style = hictrack.properties['style']

        # TODO Should we add plotting in BigWig, BedGraph, ABCCompartment, Arcs support?(The original codes supports)
        for region in tads.itertuples():
            if hicmat_style in hicmat_tri_style:
                depth = (gr.end - gr.start) / 2
                ymax = (gr.end - gr.start)
                self.plot_triangular(ax, gr, region, ymax, depth)
            elif hicmat_style in hicmat_ma_style:
                self.plot_box(ax, gr, region)
            else:
                raise ValueError(f"unsupported hicmat style {hicmat_style}")

        if len(tads) == 0:
            log.warning("No regions found for Coverage {}.".format(self.properties['name']))

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

        fill = True if self.properties['border_only'] == 'no' else False

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
        if not (isinstance(score, float) or isinstance(score, int)):
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
