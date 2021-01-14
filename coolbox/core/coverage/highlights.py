from coolbox.utilities import (
    opener, ReadBed,
    Interval, IntervalTree,
    rgb2hex,
    get_logger, GenomeRange,
    to_gr
)
from .base import Coverage

log = get_logger(__name__)


class _Highlights(object):
    DEFAULT_COLOR = '#ff5d0f'

    def fetch_data(self, gr: GenomeRange, **kwargs):
        gr = to_gr(gr)
        regions = []

        if gr.chrom not in list(self.interval_tree):
            gr.change_chrom_names()

        for region in sorted(self.interval_tree[gr.chrom][gr.start - 10000:gr.end + 10000]):
            regions.append((region.begin, region.end, region.data))

        return regions

    def plot(self, ax, gr: GenomeRange, **kwargs):

        regions = self.fetch_data(gr, **kwargs)

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


class HighLightsFromFile(Coverage, _Highlights):
    """
    High light regions coverage, read the regions from the file.

    Parameters
    ----------
    file_ : str
        Path to the file.

    color : str, optional
        High light region color,
        use 'bed_rgb' for specify color from the file, default 'bed_rgb'.

    alpha : float, optional
        High light region alpha value, default 0.1.

    border_line : bool, optional
        Plot border line or not, default True.

    border_line_style : str, optional
        Border line style, default 'dashed'.

    border_line_width : float, optional
        Border line width, default 1.0.

    border_line_color : str, optional
        Border line color, default '#000000'

    border_line_alpha : float, optional
        Border line alpha value, default 0.8

    name : str, optional
        The name of thr Coverage.
    """

    def __init__(self, file_, **kwargs):
        properties_dict = {
            "file": file_,
            "color": "bed_rgb",
            "alpha": 0.1,
            "border_line": True,
            "border_line_style": "dashed",
            "border_line_width": 0,
            "border_line_color": "#000000",
            "border_line_alpha": 0.8,
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)
        self.interval_tree = self.__process_bed()

    def __process_bed(self):
        bed = ReadBed(opener(self.properties['file']))

        if self.properties['color'] == 'bed_rgb' and bed.file_type not in ['bed12', 'bed9']:
            log.warning("*WARNING* Color set to 'bed_rgb', but bed file does not have the rgb field. The color has "
                        "been set to {}".format(HighLights.DEFAULT_COLOR))
            self.properties['color'] = HighLights.DEFAULT_COLOR

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


class HighLights(Coverage, _Highlights):
    """
    High light region.

    Parameters
    ----------
    highlight_regions : list of {str, tuple}
        A list of regions for highlights, region can be expressed as a tuple or string.
        region tuple like:
        [('chr1', 100000, 120000), ('chr2', 130000, 150000)]
        region string format: `chr:start-end` like:
        ['chr1:100000-120000', 'chr2:130000-150000'].

    color : str, optional
        High light region color, default HighLights.DEFAULT_COLOR.

    alpha : float, optional
        High light region alpha value, default 0.5

    border_line : bool, optional
        Plot border line or not, default True.

    border_line_style : str, optional
        Border line style, default 'dashed'

    border_line_width : float, optional
        Border line width, default 1.0

    border_line_color : str, optional
        Border line color, default '#000000'

    border_line_alpha : float, optional
        Border line alpha value, default 0.8

    name : str, optional
        The name of thr Coverage.
    """

    DEFAULT_COLOR = "#ff9c9c"

    def __init__(self, highlight_regions, **kwargs):
        if not isinstance(highlight_regions, list):
            highlight_regions = [highlight_regions]
        properties_dict = {
            "highlight_regions": highlight_regions,
            "color": HighLights.DEFAULT_COLOR,
            "alpha": 0.25,
            "border_line": True,
            "border_line_style": "dashed",
            "border_line_width": 0,
            "border_line_color": "#000000",
            "border_line_alpha": 0.8,
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)
        self.interval_tree = self.__intervaltree_from_list(self.properties['highlight_regions'])

    # TODO  may be duplicate of vlines's method
    def __intervaltree_from_list(self, region_list):
        from intervaltree import IntervalTree
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
            itree[chr_][grange.start:grange.end + 1] = grange
        return itree
