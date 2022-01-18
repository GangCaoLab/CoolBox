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


class PlotGenes(object):

    def __init__(self, *args, **kwargs):
        self.init_colormap()
        from matplotlib import font_manager
        properties = self.properties
        self.len_w = None  # this is the length of the letter 'w' given the font size
        self.counter = None
        self.small_relative = None
        self.is_draw_labels = properties['labels'] == 'on'
        self.fp = font_manager.FontProperties(size=properties['fontsize'])
        self.row_scale = properties['interval_height'] * 2.3
        self.cache_gr = None
        self.cache_res = None

    def fetch_plot_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        if gr == self.cache_gr:
            return self.cache_res
        else:
            self.cache_gr = gr
            self.cache_res = self.fetch_data(gr, **kwargs)
            return self.cache_res

    def get_track_height(self, frame_width, current_range):
        props = self.properties
        if (props.get('height', 'auto') == 'auto') and\
           ('row_height' in props) and\
           (props.get('display', 'stacked') == 'stacked'):
            ov_genes = self.fetch_plot_data(current_range)
            self.plot_genes(None, current_range, ov_genes, dry_run=True, fig_width=frame_width)
            return max(props['row_height'] * self.current_row_num, props['row_height'])
        else:
            try:
                height = float(self.properties['height'])
            except:
                height = 1.0
            return height

    def __set_plot_params(self, gr: GenomeRange, ov_genes: pd.DataFrame):
        properties = self.properties
        # bed_type
        self.properties['bed_type'] = properties['bed_type'] or self.infer_bed_type(ov_genes)
        self.set_colormap(ov_genes)
        # turn labels off when too many intervals are visible.
        if properties['labels'] == 'auto':
            if len(ov_genes) > 60:
                self.is_draw_labels = False
            else:
                self.is_draw_labels = True
        self.small_relative = 0.004 * (gr.end - gr.start)
        self.counter = 0

    def plot_genes(self, ax, gr: GenomeRange, ov_genes: pd.DataFrame, dry_run = False, fig_width = None):
        properties = self.properties
        self.__set_plot_params(gr, ov_genes)

        assert (not dry_run) or (fig_width is not None)
        if dry_run:
            self.__get_length_w(fig_width, gr.start, gr.end)
        else:
            self.__get_length_w(ax.get_figure().get_figwidth(), gr.start, gr.end)

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
            if not row_last_position:
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
            
            if not dry_run:
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
            log.debug(f"*Warning* No intervals were found for file {properties['file']} "
                        f"in Track \'{properties['name']}\' for the interval plotted ({gr}).\n")

        ymax = 0
        if num_rows:
            ymin = float(num_rows) * self.row_scale
            self.current_row_num = num_rows
        else:
            ymin = max_ypos + properties['interval_height']
            self.current_row_num = len(row_last_position)

        log.debug("ylim {},{}".format(ymin, ymax))
        # the axis is inverted (thus, ymax < ymin)
        if not dry_run:
            ax.set_ylim(ymin, ymax)

            if properties['display'] == 'collapsed':
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

        for idx in range(bed.block_count):
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
        for idx in range(bed.block_count):
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
        _rgb = 'grey' if first_pos[2] == 'UTR' else rgb
        vertices = self.draw_arrow(ax, first_pos[0], first_pos[1], bed.strand, ypos)

        ax.add_patch(Polygon(vertices, closed=True, fill=True,
                             edgecolor=edgecolor,
                             facecolor=_rgb,
                             linewidth=0.5))

        for start_pos, end_pos, _type in positions:
            _rgb = 'grey' if _type == 'UTR' else rgb
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

    def __get_length_w(self, fig_width, region_start, region_end):
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
            return self.properties['interval_height'] if self.counter % 2 == 0 else 1
        elif self.properties['display'] == 'collapsed':
            return 0
        else:
            return free_row * self.row_scale

