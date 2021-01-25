import pandas as pd
from matplotlib.collections import BrokenBarHCollection

from coolbox.utilities import (
    get_logger, GenomeRange, file_to_intervaltree, hex2rgb,
)
from .base import Track

log = get_logger(__name__)


# TODO change properties to new mode
class Ideogram(Track):
    """
    The chromosome ideograme track.

    Parameters
    ----------
    file_ : str
        Path to chromosome ideogram txt file,
        ideogram file is download from the UCSC Table Browser CytoBandIdeo table (in "all table" group).
        see: http://genome.ucsc.edu/cgi-bin/hgTables?hgta_group=allTables&hgta_table=cytoBandIdeo

    color_scheme : dict, optional
        Color scheme of ideogram, default: Ideogram.DEFAULT_COLOR_SCHEME

    show_name : bool, optional
        Show band name or not. default True.

    font_size : int, optional
        Band name font size.

    border_color : str, optional
        Border color. default: '#000000'

    border_width : float, optional
        Border width. default: 1.2
    """

    DEFAULT_HEIGHT = 1.2
    DEFAULT_COLOR_SCHEME = {
        'gneg': '#ffffff',
        'gpos25': '#999999',
        'gpos50': '#666666',
        'gpos75': '#333333',
        'gpos100': '#000000',
        'acen': '#cc6666',
        'gvar': '#cccccc',
        'stalk': '#e5e5e5',
    }
    DEFAULT_FONT_SIZE = 12
    DEFAULT_BORDER_WIDTH = 1.2

    def __init__(self, file_, **kwargs):
        properties_dict = {
            'file': file_,
            'color_scheme': Ideogram.DEFAULT_COLOR_SCHEME,
            'show_name': True,
            'font_size': Ideogram.DEFAULT_FONT_SIZE,
            'border_color': '#000000',
            'border_width': Ideogram.DEFAULT_BORDER_WIDTH,
            'height': Ideogram.DEFAULT_HEIGHT,
            'title': '',
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)
        self.file = self.properties['file']
        self.interval_tree, _, _ = file_to_intervaltree(self.file)

    def lookup_band_color(self, band_type):
        color_scheme = self.properties['color_scheme']
        if band_type in color_scheme:
            return color_scheme[band_type]
        else:
            return color_scheme['gneg']

    def fetch_data(self, gr: GenomeRange, **kwargs):
        if gr.chrom not in self.interval_tree:
            gr.change_chrom_names()
        bands_in_region = sorted(self.interval_tree[gr.chrom][gr.start:gr.end])
        rows = []
        for itv in bands_in_region:
            start, end = itv.begin, itv.end
            band_name, band_type = itv.data[:2]
            rows.append([gr.chrom, start, end, band_name, band_type])
        fields = ['chrom', 'start', 'end', 'name', 'gieStain']
        df = pd.DataFrame(rows, columns=fields)
        return df

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax
        df = self.fetch_data(gr)
        xranges, colors = [], []
        band_height = self.properties['height']
        for _, row in df.iterrows():
            start, end = row['start'], row['end']
            band_name, band_type = row['band_name'], row['band_type']
            band_color = self.lookup_band_color(band_type)
            xranges.append((start, end))
            colors.append(band_color)
            if self.properties['show_band_name'] != 'no':
                if gr.length < 80_000_000:
                    self.plot_text(band_name, start, end, band_color)
        coll = BrokenBarHCollection(xranges, (0, band_height), facecolors=colors,
                                    linewidths=self.properties['border_width'],
                                    edgecolors=self.properties['border_color'])
        ax.add_collection(coll)
        ax.set_ylim(-0.1, band_height + 0.1)
        ax.set_xlim(gr.start, gr.end)
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
