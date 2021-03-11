from copy import copy
from collections import OrderedDict

import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist

from coolbox.utilities import (
    cm2inch,
    GenomeRange,
    get_logger,
)
from .base import FrameBase

log = get_logger(__name__)

plt.rcParams['svg.fonttype'] = 'none'


class Frame(FrameBase):
    """
    Frame for arrange and group plots.

    Parameters
    ----------
    genome_range : {GenomeRange, str}, optional
        The init genome range.

    width : {int, float}, optional
        Width of frame, the unit is in 'cm', 1cm = 2.54inch,
        default Frame.DEFAULT_WIDTH.

    width_ratios : tuple, optional
        Width ratios of track and track title,
        like `(0.9, 0.1)`. default Frame.DEFAULT_WIDTH_RATIOS.

    margins : dict, optional
        Margins of frame, default Frame.DEFAULT_MARGINS.

    title : str, optional
        The title of this frame, default ''.

    Attributes
    ----------
    properties : dict
        Frame properties dict.

    tracks : OrderedDict
        Container of all tracks.

    current_range : GenomeRange, optional
        Current frame range.

    Examples
    --------
    >>> frame_1 = Frame()
    >>> frame_2 = Frame(gr="chr1:1000-2000")
    >>> str(frame_2.current_range)
    'chr1:1000-2000'
    >>> frame_3 = Frame(gr=GenomeRange("chr1", 1000, 2000))
    >>> str(frame_3.current_range)
    'chr1:1000-2000'
    """

    DEFAULT_WIDTH = 40
    DEFAULT_WIDTH_RATIOS = (0.01, 0.93, 0.06)
    DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0, 'top': 1}

    def __init__(self, *args, **kwargs):

        tracks = OrderedDict()
        self._tracks = tracks
        properties_dict = {
            "width": Frame.DEFAULT_WIDTH,
            "width_ratios": Frame.DEFAULT_WIDTH_RATIOS,
            "margins": Frame.DEFAULT_MARGINS,
            "title": ""
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict, *args, **kwargs)

    @property
    def tracks(self):
        return self._tracks

    def show(self):
        """
        Display self's elements on screen.
        """
        if self.current_range is None:
            raise ValueError("Frame range is uninitialized, "
                             "please `Frame.goto` method to initialize it.")

        return self.plot(self.current_range.chrom,
                         self.current_range.start,
                         self.current_range.end)

    def add_track(self, track, pos='tail'):
        """
        Add `Track` object to self.

        Parameters
        ----------
        track : Track
            The track need to be added to self.tracks

        pos : {'tail', 'head'}
            Add track to tail or head. default 'tail'

        Examples
        --------
        >>> from coolbox.core.track import XAxis, BigWig
        >>> frame = Frame()
        >>> frame.add_track(XAxis())
        >>> frame.add_track(BigWig("tests/test_data/bigwig_chrx_2e6_5e6.bw"))
        >>> len(frame.tracks)
        2
        """
        from ..track.base import Track

        assert isinstance(track, Track), "track must a Track object."
        if pos == 'tail':
            self.tracks.update({track.name: track})
        else:
            self.tracks.update({track.name: track})
            self.tracks.move_to_end(track.name, last=False)

    def get_tracks_height(self, default_height=3):
        """
        Get heights of all tracks.

        Return
        ------
        heights : list of float
            heights of all tracks.
        """
        heights = []
        for track in self.tracks.values():
            if hasattr(track, 'get_track_height'):
                frame_width = self.properties['width'] * self.properties['width_ratios'][1]
                height = track.get_track_height(frame_width)
                heights.append(height)
            elif 'height' in track.properties:
                heights.append(track.properties['height'])
            else:
                heights.append(default_height)
        return heights

    def plot(self, *args):
        """
        Plot all tracks.

        >>> from coolbox.api import *
        >>> frame = XAxis() + XAxis()
        >>> frame.plot("chr1", 100000, 200000)
        >>> frame.plot("chr1:100000-200000")
        """
        if len(args) == 3:
            gr = GenomeRange(*args[:3])
            gr2 = None
        elif len(args) == 1:
            gr = GenomeRange(args[0])
            gr2 = None
        elif len(args) == 6:
            gr = GenomeRange(args[:3])
            gr2 = GenomeRange(args[3:])
        elif len(args) == 2:
            gr = GenomeRange(args[0])
            gr2 = GenomeRange(args[1])
        elif self.current_range:
            gr = self.current_range
            gr2 = None
        else:
            raise ValueError("Please specify a genomic range in uscs format. For example: 'chr1:100000-200000'")
        # cache for the previous GenomeRange
        self.goto(gr, gr2)

        tracks_height = self.get_tracks_height()
        self.properties['height'] = sum(tracks_height)

        fig = plt.figure(figsize=cm2inch(self.properties['width'],
                                         self.properties['height']))
        if 'title' in self.properties:
            fig.suptitle(self.properties['title'])

        grids = matplotlib.gridspec.GridSpec(
            len(tracks_height), 3,
            height_ratios=tracks_height,
            width_ratios=self.properties['width_ratios'],
            wspace=0.01)

        axis_list = []
        for idx, track in enumerate(self.tracks.values()):
            y_ax = plt.subplot(grids[idx, 0])
            y_ax.set_axis_off()

            ax = axisartist.Subplot(fig, grids[idx, 1])
            fig.add_subplot(ax)
            ax.axis[:].set_visible(False)
            ax.patch.set_visible(False)

            label_ax = plt.subplot(grids[idx, 2])
            label_ax.set_axis_off()

            track.label_ax = label_ax
            track.y_ax = y_ax
            try:
                # Attention, copy is necessary, otherwise GenomeRange may change due to call of gr.change_chrom_names
                track.plot(ax, copy(gr), gr2=copy(gr2))

            except Exception as e:
                import sys, os
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = exc_tb.tb_frame.f_code.co_filename
                log.error("Error occured when plot track:\n"
                          "\ttrack name: {}\n\ttrack type:{}\n"
                          "\tError: {} {}\n"
                          "\toccurred in \"{}\", line {}".format(
                    track.name, type(track),
                    type(e), str(e), fname, exc_tb.tb_lineno)
                )
                log.exception(e)
            track.plot_coverages(ax, gr, gr2)

            axis_list.append(ax)

        margins = self.properties['margins']
        fig.subplots_adjust(wspace=0, hspace=0.0,
                            left=margins['left'],
                            right=margins['right'],
                            bottom=margins['bottom'],
                            top=margins['top'])

        plt.close()

        return fig


if __name__ == "__main__":
    import doctest

    doctest.testmod()
