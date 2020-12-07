import numpy as np

from coolbox.utilities import (
    change_chrom_names, get_logger,
)
from .bed import BED

log = get_logger(__name__)


class TADs(BED):
    """
    TADs track.

    Parameters
    ----------
    file_ : str
        Path to bed file.

    height : float, optional
        Height of track, default TADs.DEFAULT_HEIGHT

    fontsize : int, optional
        Font size, default TADs.DEFAULT_FONTSIZE

    color : str, optional
        Track color, use 'bed_rgb' for specify color according to file, default 'bed_rgb'.

    border_color : str, optional
        Border_color of gene, default 'black'.

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    title : str, optional
        Label text, default ''.

    name : str, optional
        Track's name
    """

    def __init__(self, file_, **kwargs):
        super().__init__(file_)
        properties_dict = {
            "file": file_,
            "height": TADs.DEFAULT_HEIGHT,
            "color": TADs.DEFAULT_COLOR,
            "border_color": 'black',
            "title": '',
        }
        properties_dict.update(kwargs)
        self.properties.update(properties_dict)

    def plot(self, ax, chrom_region, start_region, end_region):
        """
        Plots the boundaries as triangles in the given ax.
        """
        self.load_range(f"{chrom_region}:{start_region}-{end_region}")
        self.ax = ax

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

        self.plot_label()
