import abc
from typing import Union
from collections import OrderedDict
from copy import copy

from coolbox.utilities import (
    GenomeRange,
    op_err_msg,
    get_logger,
)

log = get_logger(__name__)

GR = Union[str, GenomeRange]


class FrameBase(abc.ABC):

    def __init__(self, properties_dict, *args, **kwargs):
        # init range
        if 'gr' in kwargs:
            range_ = kwargs['gr']
            if isinstance(range_, GenomeRange):
                self.current_range = range_
            else:
                # init from genome range string
                # e.g. `frame = Frame(gr="chr1:1000-2000")`
                self.current_range = GenomeRange(range_)
        else:
            self.current_range = None
        if 'gr2' in kwargs:
            range_ = kwargs['gr2']
            if isinstance(range_, GenomeRange):
                self.current_range = range_
            else:
                # init from genome range string
                # e.g. `frame = Frame(gr="chr1:1000-2000")`
                self.current_range = GenomeRange(range_)
        else:
            self.current_range2 = None

        self.properties = properties_dict

    @abc.abstractmethod
    def plot(self):
        pass

    @abc.abstractmethod
    def show(self):
        pass

    @property
    @abc.abstractmethod
    def tracks(self):
        pass

    @abc.abstractmethod
    def add_track(self, track):
        pass

    def goto(self, gr: GR = None, gr2: GR = None):
        """
        Go to the range on the genome.

        Parameters
        ----------
        gr2
        gr : {str, GenomeRange}
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
        if gr:
            self.current_range = gr if isinstance(gr, GenomeRange) else GenomeRange(gr)
        if gr2:
            self.current_range2 = gr if isinstance(gr2, GenomeRange) else GenomeRange(gr2)

    def fetch_data(self, gr: GR = None, gr2: GR = None):
        self.goto(gr, gr2)
        gr, gr2 = self.current_range, self.current_range2
        if gr is None:
            raise ValueError("No GenomeRange history found.")

        tracks_data = OrderedDict()
        for name, track in self.tracks.items():
            if hasattr(track, 'fetch_data'):
                data = track.fetch_data(gr, gr2=gr2)
            else:
                data = []
            tracks_data.update([(name, data)])

        return tracks_data

    def fetch_plot_data(self, gr: GR = None, gr2: GR = None):
        self.goto(gr, gr2)
        gr, gr2 = self.current_range, self.current_range2
        if gr is None:
            raise ValueError("No GenomeRange history found.")

        tracks_data = OrderedDict()
        for name, track in self.tracks.items():
            if hasattr(track, 'fetch_plot_data'):
                data = track.fetch_plot_data(gr, gr2=gr2)
            else:
                data = []
            tracks_data.update([(name, data)])

        return tracks_data

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
        from ..feature import Feature
        assert isinstance(feature, Feature), "feature must a Feature object."
        for track in self.tracks.values():
            track.properties.update(feature.properties)

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
        from ..track import BedGraph, BigWig
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
        from ..track.base import Track
        from ..feature import Feature, FrameFeature
        from ..coverage.base import Coverage, CoverageStack
        from ..browser import WidgetsPanel, Browser

        result = copy(self)

        if isinstance(other, Track):
            result.add_track(other)
            return result
        elif isinstance(other, FrameBase):
            for track in other.tracks.values():
                result.add_track(track)
            result.properties.update(other.properties)
            return result
        elif isinstance(other, FrameFeature):
            result.properties.update(other.properties)
            return result
        elif isinstance(other, Feature):
            if len(result.tracks) != 0:
                last = list(result.tracks.values())[-1]
                last.properties.update(other.properties)
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
        from ..coverage.base import Coverage
        from ..feature import Feature

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


if __name__ == "__main__":
    import doctest

    doctest.testmod()
