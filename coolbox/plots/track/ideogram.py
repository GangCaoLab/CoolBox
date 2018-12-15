from coolbox.plots.track.base import TrackPlot
from coolbox.utilities import (
    get_logger, GenomeRange, file_to_intervaltree, hex2rgb
)
from matplotlib.collections import BrokenBarHCollection


log = get_logger(__name__)


class PlotIdeogram(TrackPlot):
    DEFAULT_HEIGHT = 1.2
    DEFAULT_COLOR_SCHEME = {
        'gneg':    '#ffffff',
        'gpos25':  '#999999',
        'gpos50':  '#666666',
        'gpos75':  '#333333',
        'gpos100': '#000000',
        'acen':    '#cc6666',
        'gvar':    '#cccccc',
        'stalk':   '#e5e5e5',
    }
    DEFAULT_FONT_SIZE = 12
    DEFAULT_BORDER_COLOR = '#000000'
    DEFAULT_BORDER_WIDTH = 1.2

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)

        if 'height' not in self.properties:
            self.properties['height'] = PlotIdeogram.DEFAULT_HEIGHT
        if 'show_band_name' not in self.properties:
            self.properties['show_band_name'] = 'yes'
        if 'font_size' not in self.properties:
            self.properties['font_size'] = PlotIdeogram.DEFAULT_FONT_SIZE
        if 'color_scheme' not in self.properties:
            self.properties['color_scheme'] = PlotIdeogram.DEFAULT_COLOR_SCHEME
        if 'border_color' not in self.properties:
            self.properties['border_color'] = PlotIdeogram.DEFAULT_BORDER_COLOR
        if 'border_width' not in self.properties:
            self.properties['border_width'] = PlotIdeogram.DEFAULT_BORDER_WIDTH

        self.file = self.properties['file']
        self.interval_tree, _, _ = file_to_intervaltree(self.file)

    def lookup_band_color(self, band_type):
        color_scheme = self.properties['color_scheme']
        if band_type in color_scheme:
            return color_scheme[band_type]
        else:
            return color_scheme['gneg']

    def plot(self, ax, chrom_region, region_start, region_end):
        self.ax = ax
        grange = GenomeRange(chrom_region, region_start, region_end)
        if grange.chrom not in self.interval_tree:
            grange.change_chrom_names()
        bands_in_region = sorted(self.interval_tree[grange.chrom][grange.start:grange.end])
        band_height = self.properties['height']
        xranges, colors = [], []
        for itv in bands_in_region:
            start, end = itv.begin, itv.end
            band_name, band_type = itv.data[:2]
            band_color = self.lookup_band_color(band_type)
            xranges.append((start, end))
            colors.append(band_color)
            if self.properties['show_band_name'] != 'no':
                if grange.length < 80_000_000:
                    self.plot_text(band_name, start, end, band_color)
        coll = BrokenBarHCollection(xranges, (0, band_height), facecolors=colors,
                                    linewidths=self.properties['border_width'],
                                    edgecolors=self.properties['border_color'])
        ax.add_collection(coll)
        ax.set_ylim(-0.1, band_height+0.1)
        ax.set_xlim(region_start, region_end)
        self.plot_label()

    def plot_text(self, band_name, start, end, band_color):
        band_height = self.properties['height']
        x_pos = start + (end - start) * 0.15
        y_pos = band_height / 2
        if isinstance(band_color, str):
            rgb = hex2rgb(band_color)
        else:
            rgb = band_color
        if sum(rgb) < 100:
            color = '#e2e2e2'
        else:
            color = '#000000'
        self.ax.text(x_pos, y_pos, band_name, fontsize=self.properties['font_size'], color=color)
