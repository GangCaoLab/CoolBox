from collections import OrderedDict
from copy import copy

import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist

plt.rcParams['svg.fonttype'] = 'none'

from coolbox.utilities import (
    cm2inch,
    GenomeRange,
    op_err_msg,
    get_logger,
)

log = get_logger(__name__)


class Frame(object):
    """
    Frame for arrange and group plots.

    Parameters
    ----------
    genome_range : {GenomeRange, str}, optional
        The init genome range.

    width : {int, float}, optional
        Width of frame, default Frame.DEFAULT_WIDTH.

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
    >>> frame_2 = Frame(genome_range="chr1:1000-2000")
    >>> str(frame_2.current_range)
    'chr1:1000-2000'
    >>> frame_3 = Frame(genome_range=GenomeRange("chr1", 1000, 2000))
    >>> str(frame_3.current_range)
    'chr1:1000-2000'
    """

    DEFAULT_WIDTH = 40
    DEFAULT_WIDTH_RATIOS = (0.01, 0.93, 0.06)
    DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0, 'top': 1}

    def __init__(self, *args, **kwargs):
        # init range
        if 'genome_range' in kwargs:
            range_ = kwargs['genome_range']
            if isinstance(range_, GenomeRange):
                self.current_range = range_
            else:
                # init from genome range string
                # e.g. `frame = Frame(genome_range="chr1:1000-2000")`
                self.current_range = GenomeRange(range_)
        else:
            self.current_range = None

        tracks = OrderedDict()
        self.tracks = tracks
        properties_dict = {
            "width": Frame.DEFAULT_WIDTH,
            "width_ratios": Frame.DEFAULT_WIDTH_RATIOS,
            "margins": Frame.DEFAULT_MARGINS,
            "title": ""
        }
        properties_dict.update(kwargs)
        self.properties = properties_dict

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

    def goto(self, genome_range):
        """
        Go to the range on the genome.

        Parameters
        ----------
        genome_range : {str, GenomeRange}
            The range string,
            like "chr1:1000000-2000000", or GenomeRange object.

        Examples
        --------
        >>> frame = Frame()
        >>> frame.goto("chrX:3000000-5000000")
        >>> str(frame.current_range)
        'chrX:3000000-5000000'
        >>> frame.goto(GenomeRange("chr1", 1000, 2000))
        >>> str(frame.current_range)
        'chr1:1000-2000'
        """
        if isinstance(genome_range, GenomeRange):
            self.current_range = genome_range
        else:
            self.current_range = GenomeRange(genome_range)

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
        from .track.base import Track

        assert isinstance(track, Track), "track must a Track object."
        if pos == 'tail':
            self.tracks.update({track.name: track})
        else:
            self.tracks.update({track.name: track})
            self.tracks.move_to_end(track.name, last=False)

    def add_feature_to_tracks(self, feature):
        """
        Add feature to all tracks in this frame.

        Parameters
        ----------
        feature : Feature
            Feature object to be added to Frame's tracks.

        Examples
        --------
        >>> from coolbox.core.feature import Color
        >>> frame = Frame()
        >>> frame.add_feature_to_tracks(Color('#66ccff'))
        >>> assert all([track.properties['color'] == '#66ccff' for track in frame.tracks.values()])
        """
        from .feature import Feature
        assert isinstance(feature, Feature), "feature must a Feature object."
        for track in self.tracks.values():
            track.properties[feature.key] = feature.value

    def set_tracks_min_max(self, min_, max_, name=None):
        """
        set all bigwig and bedgraph tracks's min and max value.

        Examples
        --------
        >>> frame = Frame()
        >>> frame.set_tracks_min_max(-10, 10)
        >>> assert all([track.properties['min_value'] == -10 for track in frame.tracks.values()])
        >>> assert all([track.properties['min_value'] ==  10 for track in frame.tracks.values()])
        """
        from .track import BedGraph, BigWig
        if name is None:  # set all setable tracks
            for track in self.tracks.values():
                if isinstance(track, BedGraph) or \
                        isinstance(track, BigWig):
                    track.properties['min_value'] = min_
                    track.properties['max_value'] = max_
        else:  # set specified track
            if name not in self.tracks:
                log.warning("Track {name} not in frame")
                return
            track = self.tracks[name]
            if 'min_value' in track.properties:
                track.properties['min_value'] = min_
                track.properties['max_value'] = max_

    def add_cov_to_tracks(self, cov):
        """
        Add coverage to all tracks in this frame.

        Parameters
        ----------
        cov : Coverage
            Coverage object to be added to Frame's tracks.

        Examples
        --------
        >>> from coolbox.core.track import XAxis, BigWig
        >>> from coolbox.core.coverage import HighLights
        >>> frame = Frame()
        >>> frame.add_track(XAxis())
        >>> frame.add_track(BigWig("tests/test_data/bigwig_chrx_2e6_5e6.bw"))
        >>> highlights = HighLights([("chr1", 1000, 2000), ("chr2", 3000, 4000)])
        >>> frame.add_cov_to_tracks(highlights)
        >>> assert [track.coverages[0] is highlights for track in frame.tracks.values()]
        """
        for track in self.tracks.values():
            track.append_coverage(cov)

    def __add__(self, other):
        """
        Examples
        --------
        >>> from coolbox.api import *
        >>> frame1 = Frame()
        >>> frame2 = Frame()

        Rule: Frame + Track == Frame
        >>> frame3 = (frame1 + XAxis())
        >>> isinstance(frame3, Frame)
        True
        >>> len(frame3.tracks)
        1

        Rule: Frame + Frame == Frame
        >>> isinstance(frame1 + frame2, Frame)
        True

        Rule: Frame + Feature == Frame
        >>> f = XAxis() + BigWig("tests/test_data/bigwig_chrx_2e6_5e6.bw")
        >>> f = f + Color('#66ccff')
        >>> isinstance(f, Frame)
        True
        >>> 'color' not in list(f.tracks.values())[0].properties
        True
        >>> list(f.tracks.values())[1].properties['color']
        '#66ccff'

        Rule: Frame + Coverage == Frame
        >>> f = XAxis() + BigWig("tests/test_data/bigwig_chrx_2e6_5e6.bw")
        >>> highlight = HighLights([("chr1", 1000, 2000), ("chr2", 3000, 4000)])
        >>> f = f + highlight
        >>> isinstance(f, Frame)
        True
        >>> len(list(f.tracks.values())[0].coverages)
        0
        >>> list(f.tracks.values())[1].coverages[0] is highlight
        True

        Rule: Frame + CoverageStack == Frame
        >>> cov_stack = HighLights([("chr1", 1000, 2000)]) + HighLights([("chr1", 4000, 6000)])
        >>> f = XAxis() + XAxis()
        >>> f = f + cov_stack
        >>> isinstance(f, Frame)
        True
        >>> len(list(f.tracks.values())[0].coverages)
        0
        >>> len(list(f.tracks.values())[1].coverages)
        2

        Rule: Frame + WidgetsPanel == Browser
        >>> f = XAxis() + XAxis()
        >>> bsr = f + WidgetsPanel()
        >>> isinstance(bsr, Browser)
        True

        >>> f + 1 # error operation
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand type(s) for +: '<class 'coolbox.core.frame.Frame'>' and '<class 'int'>'


        """
        from .track.base import Track
        from .feature import Feature, FrameFeature
        from .coverage.base import Coverage, CoverageStack
        from .browser import WidgetsPanel, Browser

        result = copy(self)

        if isinstance(other, Track):
            result.add_track(other)
            return result
        elif isinstance(other, Frame):
            for track in other.tracks.values():
                result.add_track(track)
            result.properties.update(other.properties)
            return result
        elif isinstance(other, FrameFeature):
            result.properties[other.key] = other.value
            return result
        elif isinstance(other, Feature):
            if len(result.tracks) != 0:
                last = list(result.tracks.values())[-1]
                last.properties[other.key] = other.value
            return result
        elif isinstance(other, Coverage):
            if len(result.tracks) != 0:
                last = list(result.tracks.values())[-1]
                last.append_coverage(other)
            return result
        elif isinstance(other, CoverageStack):
            if len(result.tracks) != 0:
                last = list(result.tracks.values())[-1]
                last.pile_coverages(other.coverages, pos='top')
            return result
        elif isinstance(other, WidgetsPanel):
            return Browser(self, reference_genome=other.ref, widgets_box=other.type)
        else:
            raise TypeError(op_err_msg(self, other))

    def __mul__(self, other):
        """
        Examples
        --------
        >>> from coolbox.core.track import XAxis, BigWig
        >>> from coolbox.core.coverage import HighLights
        >>> from coolbox.core.feature import Color

        Rule: Frame * Feature == Frame
        >>> f = XAxis() + BigWig("tests/test_data/bigwig_chrx_2e6_5e6.bw")
        >>> f = f * Color("#ff9c9c")
        >>> assert all([track.properties['color'] == '#ff9c9c' for track in f.tracks.values()])

        Rule: Frame * Coverage == Frame
        >>> cov = HighLights([('chr1', 1000, 2000)])
        >>> f = f * cov
        >>> assert all([track.coverages[0] is cov for track in f.tracks.values()])

        """
        from .coverage.base import Coverage
        from .feature import Feature

        if isinstance(other, Coverage):
            result = copy(self)
            result.add_cov_to_tracks(other)
            return result
        elif isinstance(other, Feature):
            result = copy(self)
            result.add_feature_to_tracks(other)
            return result
        else:
            raise TypeError(op_err_msg(self, other, op='*'))

    def fetch_data(self, genome_range=None):
        if genome_range is None:
            genome_range = self.current_range

        if isinstance(genome_range, str):
            genome_range = GenomeRange(genome_range)

        tracks_data = OrderedDict()
        for name, track in self.tracks.items():
            if hasattr(track, 'fetch_data'):
                data = track.fetch_data(genome_range)
            else:
                data = []
            tracks_data.update([(name, data)])

        return tracks_data

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
        if len(args) >= 3:
            chrom, start, end = args[:3]
        elif len(args) >= 1:
            region_str = args[0]
            gr = GenomeRange(region_str)
            chrom, start, end = (gr.chrom, gr.start, gr.end)
        else:
            raise ValueError("Please specify a genomic range in uscs format. For example: 'chr1:100000-200000'")

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
                track.plot(ax, chrom, start, end)

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
            # plot coverages
            if hasattr(track, 'coverages'):
                for cov_idx, cov in enumerate(track.coverages):
                    cov.track = track
                    try:
                        cov.plot(ax, chrom, start, end)
                    except Exception as e:
                        log.error("Error occured when plot track's coverage:\n"
                                  "\ttrack name: {}\n\ttrack type:{}\n\tcoverage name: {}\n\tcov type: {}\n"
                                  "\tError: {} {}".format(
                            track.name, type(track), cov.name, type(cov),
                            type(e), str(e)))
                        log.exception(e)

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
