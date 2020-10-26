from coolbox.utilities import (
    opener, ReadBed,
    Interval, IntervalTree,
    change_chrom_names, rgb2hex,
    get_logger, GenomeRange
)

from coolbox.plots.coverage.base import CoveragePlot
from coolbox.plots.track.bed import PlotBed


log = get_logger(__name__)


class PlotHighLightRegions(CoveragePlot):

    DEFAULT_COLOR = '#ff5d0f'
    DEFAULT_ALPHA = 0.6
    DEFAULT_BORDER_LINE_COLOR = '#000000'
    DEFAULT_BORDER_LINE_ALPHA = 0.8
    DEFAULT_BORDER_LINE_WIDTH = 0.5
    DEFAULT_BORDER_LINE_STYLE = 'dashed'

    def __init__(self, *args, **kwargs):

        CoveragePlot.__init__(self, *args, **kwargs)

        if 'color' not in self.properties:
            self.properties['color'] = PlotHighLightRegions.DEFAULT_COLOR
        if 'alpha' not in self.properties:
            self.properties['alpha'] = PlotHighLightRegions.DEFAULT_BORDER_LINE_ALPHA
        if 'border_line' not in self.properties:
            self.properties['border_line'] = 'yes'
        if 'border_line_color' not in self.properties:
            self.properties['border_line_color'] = PlotHighLightRegions.DEFAULT_BORDER_LINE_COLOR
        if 'border_line_alpha' not in self.properties:
            self.properties['border_line_alpha'] = PlotHighLightRegions.DEFAULT_BORDER_LINE_ALPHA
        if 'border_line_width' not in self.properties:
            self.properties['border_line_width'] = PlotHighLightRegions.DEFAULT_BORDER_LINE_WIDTH
        if 'border_line_style' not in self.properties:
            self.properties['border_line_style'] = PlotHighLightRegions.DEFAULT_BORDER_LINE_STYLE

        if 'file' in self.properties:
            # from bed file
            self.interval_tree = self.__process_bed()
        else:
            # from self.properties['regions']
            self.interval_tree = self.__intervaltree_from_list(self.properties['highlight_regions'])

    def __intervaltree_from_list(self, region_list):
        itree = {}
        for r in region_list:
            if isinstance(r, str):
                grange = GenomeRange(r)
            elif isinstance(r, tuple):
                grange = GenomeRange(r[0], r[1], r[2])
            elif isinstance(r, GenomeRange):
                grange = r
            else:
                raise ValueError("position must be a tuple or string.")
            chr_ = grange.chrom
            itree.setdefault(chr_, IntervalTree())
            itree[chr_][grange.start:grange.end+1] = grange
        return itree

    def __get_regions(self, chrom, start, end):
        regions = []

        if chrom not in list(self.interval_tree):
            chrom = change_chrom_names(chrom)

        for region in sorted(self.interval_tree[chrom][start-10000:end+10000]):
            regions.append((region.begin, region.end, region.data))

        return regions

    def plot(self, ax, chrom_region, start_region, end_region):

        regions = self.__get_regions(chrom_region, start_region, end_region)

        for (start, end, color) in regions:
            if self.properties['color'] != 'bed_rgb':
                color = self.properties['color']
            if type(color) is not str:
                color = rgb2hex(*color)

            ax.axvspan(start, end, color=color, alpha=self.properties['alpha'])

            if self.properties['border_line'] == 'yes':
                # plot border line
                ymin, ymax = ax.get_ylim()
                ax.vlines([start, end], ymin, ymax,
                          linestyle=self.properties['border_line_style'],
                          linewidth=self.properties['border_line_width'],
                          color=self.properties['border_line_color'],
                          alpha=self.properties['border_line_alpha'])

    def __process_bed(self):
        bed = ReadBed(opener(self.properties['file']))

        if self.properties['color'] == 'bed_rgb' and bed.file_type not in ['bed12', 'bed9']:
            log.warning("*WARNING* Color set to 'bed_rgb', but bed file does not have the rgb field. The color has "
                        "been set to {}".format(PlotBed.DEFAULT_COLOR))
            self.properties['color'] = PlotBed.DEFAULT_COLOR

        interval_tree = {}

        for intval in bed:

            if intval.chromosome not in interval_tree:
                interval_tree[intval.chromosome] = IntervalTree()

            if self.properties['color'] == 'bed_rgb':
                color = intval.rgb
            else:
                color = self.properties['color']
            interval_tree[intval.chromosome].add(Interval(intval.start, intval.end, color))

        return interval_tree
