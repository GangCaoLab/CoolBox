from coolbox.utilities import (
    opener, ReadBed,
    Interval, IntervalTree,
    change_chrom_names, rgb2hex,
    get_logger
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
            self.regions = self.properties['highlight_regions']

    def __extract_regions_from_file(self, chrom, start, end):
        regions = []

        if chrom not in list(self.interval_tree):
            chrom = change_chrom_names(chrom)

        for region in sorted(self.interval_tree[chrom][start-10000:end+10000]):
            regions.append((region.begin, region.end, region.data))

        return regions

    def __get_regions_from_properties(self, chrom):

        def regions_with_color(regions):
            color = self.properties['color']
            res = []
            for r in self.regions:
                assert len(r) >= 2
                start, end = r[0], r[1]
                if len(r) == 2:
                    res.append((start, end, color))
                else:
                    res.append(r)
            return res

        if self.properties['chr'] is not None:
            if chrom == self.properties['chr']:
                regions = regions_with_color(self.regions)
            else:
                regions = []
        else:
            regions = regions_with_color(self.regions)

        return regions

    def plot(self, ax, chrom_region, start_region, end_region):
        if 'file' in self.properties:
            regions = self.__extract_regions_from_file(chrom_region, start_region, end_region)
        else:
            regions = self.__get_regions_from_properties(chrom_region)

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
