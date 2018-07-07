import matplotlib
import matplotlib.colors as colors
from matplotlib.patches import Rectangle

import numpy as np

from coolbox.utilities import (
    change_chrom_names,
    Interval, IntervalTree,
    opener, ReadBed,
    get_logger
)

from coolbox.plots.track.base import TrackPlot

log = get_logger(__name__)


class PlotBed(TrackPlot):

    DEFAULT_COLOR = "#1f78b4"

    def __init__(self, *args, **kwarg):
        TrackPlot.__init__(self, *args, **kwarg)
        self.bed_type = None  # once the bed file is read, this is bed3, bed6 or bed12
        self.len_w = None  # this is the length of the letter 'w' given the font size
        self.interval_tree = {}  # interval tree of the bed regions

        from matplotlib import font_manager
        if 'fontsize' not in self.properties:
            self.properties['fontsize'] = 12
        else:
            self.properties['fontsize'] = float(self.properties['fontsize'])

        self.fp = font_manager.FontProperties(size=self.properties['fontsize'])

        if 'color' not in self.properties:
            self.properties['color'] = 'bed_rgb'
        if 'border_color' not in self.properties:
            self.properties['border_color'] = 'black'
        if 'labels' not in self.properties:
            self.properties['labels'] = 'auto'
        if 'style' not in self.properties:
            self.properties['style'] = 'flybase'
        if 'display' not in self.properties:
            self.properties['display'] = 'stacked'
        if 'interval_height' not in self.properties:
            self.properties['interval_height'] = 100

        if self.properties['labels'] != 'on':
            self.is_draw_labels = False
        else:
            self.is_draw_labels = True

        self.colormap = None
        # check if the color given is a color map
        if not matplotlib.colors.is_color_like(self.properties['color']) and self.properties['color'] != 'bed_rgb':
            # check if the color is a valid colormap name
            if self.properties['color'] not in matplotlib.cm.datad:
                log.warning("*WARNING* color: '{}' for Track {} is not valid. Color has "
                            "been set to {}".format(self.properties['color'], self.properties['name'],
                                                    PlotBed.DEFAULT_COLOR))
                self.properties['color'] = PlotBed.DEFAULT_COLOR
            else:
                self.colormap = self.properties['color']

        # to set the distance between rows
        self.row_scale = self.properties['interval_height'] * 2.3

        self.interval_tree, min_score, max_score = self.__process_bed()
        if self.colormap is not None:
            if 'min_value' in self.properties:
                min_score = self.properties['min_value']
            if 'max_value' in self.properties:
                max_score = self.properties['max_value']

            norm = matplotlib.colors.Normalize(vmin=min_score,
                                               vmax=max_score)

            cmap = matplotlib.cm.get_cmap(self.properties['color'])
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    def __get_length_w(self, fig_width, region_start, region_end):
        '''
        to improve the visualization of the genes it is good to have an estimation of the label
        length. In the following code I try to get the length of a 'W' in base pairs.
        '''
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

    def __process_bed(self):

        bed_file_h = ReadBed(opener(self.properties['file']))
        self.bed_type = bed_file_h.file_type

        if 'color' in self.properties and self.properties['color'] == 'bed_rgb' and \
                self.bed_type not in ['bed12', 'bed9']:
            log.warning("*WARNING* Color set to 'bed_rgb', but bed file does not have the rgb field. The color has "
                        "been set to {}".format(PlotBed.DEFAULT_COLOR))
            self.properties['color'] = PlotBed.DEFAULT_COLOR

        valid_intervals = 0
        interval_tree = {}

        max_score = float('-inf')
        min_score = float('inf')
        for bed in bed_file_h:
            if bed.score < min_score:
                min_score = bed.score
            if bed.score > max_score:
                max_score = bed.score

            if bed.chromosome not in interval_tree:
                interval_tree[bed.chromosome] = IntervalTree()

            interval_tree[bed.chromosome].add(Interval(bed.start, bed.end, bed))
            valid_intervals += 1

        if valid_intervals == 0:
            log.warning("No valid intervals were found in file {}".format(self.properties['file_name']))

        return interval_tree, min_score, max_score

    def __get_max_num_row(self, len_w, small_relative):
        ''' Process the whole bed regions at the given figure length and font size to
        determine the maximum number of rows required.
        :return:
        '''

        self.max_num_row = {}
        for chrom in self.interval_tree:
            row_last_position = []  # each entry in this list contains the end position
            self.max_num_row[chrom] = 0
            for region in sorted(self.interval_tree[chrom][0:500000000]):
                bed = region.data
                if self.is_draw_labels:
                    bed_extended_end = int(bed.end + (len(bed.name) * len_w))
                else:
                    bed_extended_end = (bed.end + 2 * small_relative)

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

                if free_row > self.max_num_row[bed.chromosome]:
                    self.max_num_row[bed.chromosome] = free_row

        log.debug("max number of rows set to {}".format(self.max_num_row))
        return self.max_num_row

    def __get_y_pos(self, free_row):
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

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.counter = 0
        self.small_relative = 0.004 * (end_region - start_region)
        self.__get_length_w(ax.get_figure().get_figwidth(), start_region, end_region)
        if 'global_max_row' in self.properties and self.properties['global_max_row'] == 'yes':
            self.__get_max_num_row(self.len_w, self.small_relative)

        if chrom_region not in self.interval_tree.keys():
            chrom_region = change_chrom_names(chrom_region)

        genes_overlap = sorted(self.interval_tree[chrom_region][start_region:end_region])

        if self.properties['labels'] == 'auto':
            if len(genes_overlap) > 60:
                # turn labels off when too many intervals are visible.
                self.is_draw_labels = False
            else:
                self.is_draw_labels = True

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

        for region in genes_overlap:
            """
            BED12 gene format with exon locations at the end
            chrX    20850   23076   CG17636-RA      0       -       20850   23017   0       3       946,765,64,     0,1031,2162,

            BED9
            bed with rgb at end
            chr2L   0       70000   ID_5    0.26864549832   .       0       70000   51,160,44

            BED6
            bed without rgb
            chr2L   0       70000   ID_5    0.26864549832   .
            """
            self.counter += 1
            bed = region.data

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

            ypos = self.__get_y_pos(free_row)

            # do not plot if the maximum interval rows to plot is reached
            if 'gene_rows' in self.properties and free_row >= int(self.properties['gene_rows']):
                continue

            if free_row > max_num_row_local:
                max_num_row_local = free_row
            if ypos > max_ypos:
                max_ypos = ypos

            if self.bed_type == 'bed12':
                if self.properties['style'] == 'flybase':
                    self.__draw_gene_with_introns_flybase_style(ax, bed, ypos, rgb, edgecolor)
                else:
                    self.__draw_gene_with_introns(ax, bed, ypos, rgb, edgecolor)
            else:
                self.__draw_gene_simple(ax, bed, ypos, rgb, edgecolor)

            if not self.is_draw_labels:
                pass
            elif bed.start > start_region and bed.end < end_region:
                ax.text(bed.end + self.small_relative, ypos + (float(self.properties['interval_height']) / 2),
                        bed.name, horizontalalignment='left',
                        verticalalignment='center', fontproperties=self.fp)

        if self.counter == 0:
            log.warning("*Warning* No intervals were found for file {} "
                        "in Track '{}' for the interval plotted ({}:{}-{}).\n".
                        format(self.properties['file'], self.properties['name'], chrom_region, start_region, end_region))

        ymax = 0

        if 'global_max_row' in self.properties and self.properties['global_max_row'] == 'yes':
            ymin = self.max_num_row[chrom_region] * self.row_scale

        elif 'gene_rows' in self.properties:
            ymin = int(self.properties['gene_rows']) * self.row_scale
        else:
            ymin = max_ypos + self.properties['interval_height']

        log.debug("ylim {},{}".format(ymin, ymax))
        # the axis is inverted (thus, ymax < ymin)
        ax.set_ylim(ymin, ymax)

        if 'display' in self.properties:
            if self.properties['display'] == 'domain':
                ax.set_ylim(-5, 205)
            elif self.properties['display'] == 'collapsed':
                ax.set_ylim(-5, 105)

        ax.set_xlim(start_region, end_region)

        label_ax.text(0.15, 1, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='top', transform=label_ax.transAxes)

    def get_rgb_and_edge_color(self, bed):
        rgb = self.properties['color']
        edgecolor = self.properties['border_color']

        if self.colormap:
            # translate value field (in the example above is 0 or 0.2686...) into a color
            rgb = self.colormap.to_rgba(bed.score)

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
                    rgb = PlotBed.DEFAULT_COLOR
            else:
                rgb = PlotBed.DEFAULT_COLOR
        return rgb, edgecolor

    def __draw_gene_simple(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws an interval with direction (if given)
        """
        from matplotlib.patches import Polygon

        if bed.strand not in ['+', '-']:
            ax.add_patch(Rectangle((bed.start, ypos), bed.end - bed.start, self.properties['interval_height'],
                                   edgecolor=edgecolor, facecolor=rgb, linewidth=0.5))
        else:
            vertices = self.__draw_arrow(ax, bed.start, bed.end, bed.strand, ypos)
            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=rgb,
                                 linewidth=0.5))

    def __draw_gene_with_introns_flybase_style(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws a gene using different styles
        """
        from matplotlib.patches import Polygon
        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self.__draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
            return
        half_height = float(self.properties['interval_height']) / 2
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
                    type = 'UTR'
                else:
                    type = 'coding'

                positions.append((x0, x1, type))

        # plot all blocks as rectangles except the last if the strand is + or
        # the first is the strand is -, which are drawn as arrows.
        if bed.strand == '-':
            positions = positions[::-1]

        first_pos = positions.pop()
        if first_pos[2] == 'UTR':
            _rgb = 'grey'
        else:
            _rgb = rgb

        vertices = self.__draw_arrow(ax, first_pos[0], first_pos[1], bed.strand, ypos)

        ax.add_patch(Polygon(vertices, closed=True, fill=True,
                             edgecolor=edgecolor,
                             facecolor=_rgb,
                             linewidth=0.5))

        for start_pos, end_pos, _type in positions:
            if _type == 'UTR':
                _rgb = 'grey'
            else:
                _rgb = rgb
            vertices = [(start_pos, ypos), (start_pos, ypos + self.properties['interval_height']),
                        (end_pos, ypos + self.properties['interval_height']), (end_pos, ypos)]

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=_rgb,
                                 linewidth=0.5))

    def __draw_arrow(self, ax, start, end, strand, ypos):
        half_height = float(self.properties['interval_height']) / 2
        if strand == '+':
            x0 = start
            x1 = end  # - self.small_relative
            y0 = ypos
            y1 = ypos + self.properties['interval_height']

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
            y1 = ypos + self.properties['interval_height']

            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

            /-----------------
            \-----------------
            """

            vertices = [(x0, y0), (x0 - self.small_relative, y0 + half_height), (x0, y1), (x1, y1), (x1, y0)]

        return vertices

    def __draw_gene_with_introns(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws a gene like in flybase gbrowse.
        """
        from matplotlib.patches import Polygon

        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self.__draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
            return
        half_height = float(self.properties['interval_height']) / 2
        quarter_height = float(self.properties['interval_height']) / 4
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
                y1 = ypos + self.properties['interval_height']

            if x0 < bed.thick_start < x1:
                vertices = ([(x0, ypos + quarter_height), (x0, ypos + three_quarter_height),
                             (bed.thick_start, ypos + three_quarter_height),
                             (bed.thick_start, ypos + self.properties['interval_height']),
                             (bed.thick_start, ypos + self.properties['interval_height']),
                             (x1, ypos + self.properties['interval_height']), (x1, ypos),
                             (bed.thick_start, ypos), (bed.thick_start, ypos + quarter_height)])

            elif x0 < bed.thick_end < x1:
                vertices = ([(x0, ypos),
                             (x0, ypos + self.properties['interval_height']),
                             (bed.thick_end, ypos + self.properties['interval_height']),
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


class PlotTADs(PlotBed):

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        """
        Plots the boundaries as triangles in the given ax.
        """
        from matplotlib.patches import Polygon
        ymax = 0.001
        valid_regions = 0
        if chrom_region not in self.interval_tree:
            orig = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            log.info('Chromosome name: {} does not exists. Changing name to {}'.format(orig, chrom_region))

        for region in sorted(self.interval_tree[chrom_region][start_region:end_region]):
            """
                  /\
                 /  \
                /    \
            _____________________
               x1 x2 x3
            """
            x1 = region.begin
            x2 = x1 + float(region.end - region.begin) / 2
            x3 = region.end
            y1 = 0
            y2 = (region.end - region.begin)

            rgb, edgecolor = self.get_rgb_and_edge_color(region.data)

            triangle = Polygon(np.array([[x1, y1], [x2, y2], [x3, y1]]), closed=True,
                               facecolor=rgb, edgecolor=edgecolor)
            ax.add_artist(triangle)
            valid_regions += 1

            if y2 > ymax:
                ymax = y2

        if valid_regions == 0:
            log.warning("No regions found for Track {}.".format(self.properties['name']))

        ax.set_xlim(start_region, end_region)
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(ymax, 0)
        else:
            ax.set_ylim(0, ymax)

        label_ax.text(0.15, 0.5, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='center')

