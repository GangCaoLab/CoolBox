"""
CoolBox API

example:

``` Python
from coolbox.api import *

bsr = Browser()

# open files
bsr.open("GM12878-rna-1.bw")
bsr.open("GM12878-chip-h3k27ac-1.bw")
bsr.open("GM12878-Hind3-1.cool")

bsr.show()

```

Or, using another style(ggplot like) API.

``` Python
from coolbox.api import *

frame = XAxis() + \
        BigWig("GM12878-rna-1.bw) + Color("#66ccff") + \
        BigWig("GM12878-chip-h3k27ac-1.bw") + Color("#ff9c9c") + \
        Cool("GM12878-Hind3-1.cool") + ColorMap("#0000ff", "#ff0000")

bsr = Browser(frame)

bsr.show()

```


Element types:

* Track
* Coverage
* Frame
* Feature
    - FrameFeature
* Browser


The rule of element composition:
```
    Track + Track = Frame
    Track + Feature = Track
    Track + Coverage = Track
    Frame + Track = Frame
    Frame + Coverage = Frame
    Frame + Feature = Frame
    Frame + FrameFeature = Frame
    Frame + Frame = Frame
    Frame + WidgetsPanel = Browser
    Coverage + Feature = Coverage

    Frame * Feature = Frame
    Frame * Coverage = Frame
```
"""

import warnings
warnings.filterwarnings('ignore')

from copy import copy
from collections import OrderedDict, deque

from .plots import *
from .fetchdata import *

from .interact import BrowserBase

from .utilities import GenomeRange, op_err_msg


import logging
log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)

logging.getLogger('coolbox.utilities').setLevel(logging.WARNING)
logging.getLogger('coolbox.plots').setLevel(logging.WARNING)
logging.getLogger('coolbox.fetchdata').setLevel(logging.WARNING)



__all__ = [
    "Frame", "WidgetsPanel", "Browser",
    "Feature", "Color", "ColorMap", "TrackHeight", "Inverted",
    "Title", "MaxValue", "MinValue", "ShowDataRange",
    "FrameFeature", "FrameTitle",
    "Track", "Spacer", "XAxis", "Bed", "BigWig", "ABCompartment", "BedGraph",
    "Boundaries", "Arcs", "TADs", "Cool",
    "Coverage", "CoverageStack", "VlinesFromFile", "Vlines",
    "HighLightsFromFile", "HighLights", "HiCPeaks",
]


FEATURES_STACK_NAME = "__COOLBOX_FEATURE_STACK__"
COVERAGE_STACK_NAME = "__COOLBOX_COVERAGE_STACK__"
global_scope = globals()
global_scope[FEATURES_STACK_NAME] = deque()
global_scope[COVERAGE_STACK_NAME] = deque()


class Frame(PlotFrame, FetchFrame):
    """
    Frame for arrange and group plots.

    Attributes:
        properties (:obj:`dict`, (`str`, `any type`)): frame properties dict.
        tracks (:obj:`OrderedDict`, key `str`, value :obj:`Track`): container of all tracks.
        current_range (:obj:`GenomeRange`, optional): current frame range.

    """

    DEFAULT_WIDTH = 40
    DEFAULT_WIDTH_RATIOS = (0.93, 0.07)
    DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0, 'top': 1}

    def __init__(self, *args, **kwargs):
        """
        Args:
            genome_range (:obj:`GenomeRange` or `str`, optional):
                the init genome range.
            width (`int` or `float`): width of frame. [Frame.DEFAULT_WIDTH]
            width_ratios (:obj: tuple): width ratios of track and track title,
                like `(0.9, 0.1)`. [Frame.DEFAULT_WIDTH_RATIOS]
            margins (:dict:): margins of frame. [Frame.DEFAULT_MARGINS]
            title (str): the title of this frame.

        >>> frame_1 = Frame()
        >>> frame_2 = Frame(genome_range="chr1:1000-2000")
        >>> str(frame_2.current_range)
        'chr1:1000-2000'
        >>> frame_3 = Frame(genome_range=GenomeRange("chr1", 1000, 2000))
        >>> str(frame_3.current_range)
        'chr1:1000-2000'
        """

        super().__init__({}, OrderedDict())

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

        # set properties
        if 'width' in kwargs:
            self.properties['width'] = kwargs['width']
        else:
            self.properties['width'] = Frame.DEFAULT_WIDTH

        if 'width_ratios' in kwargs:
            self.properties['width_ratios'] = kwargs['width_ratios']
        else:
            self.properties['width_ratios'] = Frame.DEFAULT_WIDTH_RATIOS

        if 'margins' in kwargs:
            self.properties['margins'] = kwargs['margins']
        else:
            self.properties['margins'] = Frame.DEFAULT_MARGINS

        if 'title' in kwargs:
            self.properties['title'] = kwargs['title']

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

        Args:
            genome_range (`str` or :obj:`GenomeRange`): the range string,
                like "chr1:1000000-2000000", or GenomeRange object.

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

        Args:
            track (:obj:`Track`): the track need to be added to self.tracks
            pos (`str`, 'tail' or 'head'): add track to tail or head. ['tail']

        >>> frame = Frame()
        >>> frame.add_track(XAxis())
        >>> frame.add_track(BigWig("tests/test_data/bigwig_chrx_2e6_5e6.bw"))
        >>> len(frame.tracks)
        2
        """
        assert isinstance(track, Track), "track must a Track object."
        if pos == 'tail':
            self.tracks.update({track.name: track})
        else:
            self.tracks.update({track.name: track})
            self.tracks.move_to_end(track.name, last=False)
    
    def add_feature_to_tracks(self, feature):
        """
        Add feature to all tracks in this frame.

        Args:
            feature (:obj:`Feature`): Feature object to be added to Frame's tracks.

        >>> frame = Frame()
        >>> frame.add_feature_to_tracks(Color('#66ccff'))
        >>> assert all([track.properties['color'] == '#66ccff' for track in frame.tracks.values()])
        """
        assert isinstance(feature, Feature), "feature must a Feature object."
        for track in self.tracks.values():
            track.properties[feature.key] = feature.value

    def set_tracks_min_max(self, min_, max_):
        """
        set all bigwig and bedgraph tracks's min and max value.

        >>> frame = Frame()
        >>> frame.set_tracks_min_max(-10, 10)
        >>> assert all([track.properties['min_value'] == -10 for track in frame.tracks.values()])
        >>> assert all([track.properties['min_value'] ==  10 for track in frame.tracks.values()])
        """
        for track in self.tracks.values():
            if isinstance(track, BedGraph) or \
               isinstance(track, BigWig):
                track.properties['min_value'] = min_
                track.properties['max_value'] = max_

    def add_cov_to_tracks(self, cov):
        """
        Add coverage to all tracks in this frame.

        Args:
            cov (:obj:`Coverage`): Coverage object to be added to Frame's tracks.

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
        TypeError: unsupported operand type(s) for +: '<class '__main__.Frame'>' and '<class 'int'>'


        """
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
            return Browser(self)
        else:
            raise TypeError(op_err_msg(self, other))

    def __mul__(self, other):
        """

        Rule: Frame * Feature == Frame
        >>> f = XAxis() + BigWig("tests/test_data/bigwig_chrx_2e6_5e6.bw")
        >>> f = f * Color("#ff9c9c")
        >>> assert all([track.properties['color'] == '#ff9c9c' for track in f.tracks.values()])

        Rule: Frame * Coverage == Frame
        >>> cov = HighLights([('chr1', 1000, 2000)])
        >>> f = f * cov
        >>> assert all([track.coverages[0] is cov for track in f.tracks.values()])

        """

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


class WidgetsPanel(object):
    """
    Widgets container.
    """

    def __init__(self, reference_genome="hg19"):
        self.ref = reference_genome
    
    def __add__(self, other):
        """

        >>>

        """
        if isinstance(other, Frame):
            return Browser(other, reference_genome=self.ref)
        else:
            raise TypeError(op_err_msg(self, other, op="+"))


class Browser(BrowserBase):
    """
    Genoeme browser.
    include:
        * Frame: for show plots
        * widgetsPanel: for control genome region showed in Frame.
    """

    @property
    def tracks(self):
        return self.frame.tracks

    def fetch_data(self, genome_range=None):
        return self.frame.fetch_data(genome_range)


class Feature(object):
    """
    Feature base class.
    """

    def __init__(self, key, value):
        self.key = key
        self.value = value

    def __add__(self, other):
        if isinstance(other, Track):
            result = copy(other)
            result.properties[self.key] = self.value
            return result
        elif isinstance(other, Frame):
            result = copy(other)
            if len(result.tracks) != 0:
                first = list(result.tracks.values())[0]
                first.properties[self.key] = self.value
            return result
        elif isinstance(other, Coverage):
            result = copy(other)
            result.properties[self.key] = self.value
            return result
        elif isinstance(other, CoverageStack):
            result = copy(other)
            if len(result.coverages) != 0:
                first = result.coverages[0]
                first[self.key] = self.value
            return result
        else:
            raise TypeError(op_err_msg(self, other))

    def __mul__(self, other):
        if isinstance(other, Frame):
            result = copy(other)
            result.add_feature_to_tracks(self)
            return result
        else:
            raise TypeError(op_err_msg(self, other, op='*'))

    def __enter__(self):
        scope = globals()
        stack = scope[FEATURES_STACK_NAME]
        stack.append(self)
        return self

    def __exit__(self, type, value, traceback):
        scope = globals()
        stack = scope[FEATURES_STACK_NAME]
        stack.pop()


class Color(Feature):
    """
    Track color.
    """
    def __init__(self, value):
        super().__init__('color', value)


class ColorMap(Feature):
    """
    Track color map.
    """
    def __init__(self, value):
        super().__init__('color', value)


class TrackHeight(Feature):
    """
    Track height.
    """
    def __init__(self, value):
        super().__init__('height', value)


class Inverted(Feature):
    """
    Invert the orientation of track.
    """
    def __init__(self):
        super().__init__('orientation', 'inverted')


class Title(Feature):
    """
    Track title.
    """
    def __init__(self, value):
        super().__init__('title', value)


class MaxValue(Feature):
    """
    Max value of track.
    """
    def __init__(self, value):
        super().__init__('max_value', value)


class MinValue(Feature):
    """
    Min value of track.
    """
    def __init__(self, value):
        super().__init__('min_value', value)


class ShowDataRange(Feature):
    """
    Show data range or not
    """
    def __init__(self, show=True):
        if show:
            super().__init__('show_data_range', 'yes')
        else:
            super().__init__('show_data_range', 'no')


class FrameFeature(Feature):
    """
    FrameFeature base class.
    """

    def __add__(self, other):
        if isinstance(other, Frame):
            result = copy(other)
            result.properties[self.key] = self.value
            return result
        else:
            raise TypeError(op_err_msg(self, other))


class FrameTitle(FrameFeature):
    """
    Frame title.
    """
    def __init__(self, value):
        super().__init__('title', value)


class Track(object):
    """
    Track base class.

    Attributes:
        name (str): The name of this Track.
        coverages (:obj:`list` of Coverage): Coverages on this Track.
    """

    DEFAULT_HEIGHT = 3
    DEFAULT_COLOR = "#808080"

    def __new__(cls, *args, **kwargs):
        if hasattr(cls, "_counts"):
            cls._counts += 1
        else:
            cls._counts = 1
        return super().__new__(cls)

    def __init__(self, properties_dict, name=None):
        self.properties = properties_dict
        if name is not None:
            assert isinstance(name, str), "Track name must be a `str`."
        else:
            name = self.__class__.__name__ + ".{}".format(self.__class__._counts)
        self.properties['name'] = name
        super().__init__()
        self.coverages = []

        scope = globals()
        features_stack = scope[FEATURES_STACK_NAME]
        for features in features_stack:
            self.properties[features.key] = features.value

        coverage_stack = scope[COVERAGE_STACK_NAME]
        for coverage in coverage_stack:
            self.coverages.append(coverage)

    def __del__(self):
        self.__class__._counts -= 1

    @property
    def name(self):
        return self.properties['name']
    
    @name.setter
    def name(self, value):
        self.properties['name'] = value

    def __add__(self, other):
        if isinstance(other, Track):
            result = Frame()
            result.add_track(self)
            result.add_track(other)
            return result
        elif isinstance(other, Frame):
            result = copy(other)
            result.add_track(self, pos='head')
            return result
        elif isinstance(other, Feature):
            result = copy(self)
            return other.__add__(result)
        elif isinstance(other, Coverage):
            result = copy(self)
            result.append_coverage(other)
            return result
        elif isinstance(other, CoverageStack):
            result = copy(self)
            result.pile_coverages(other.coverages, pos='top')
            return result
        else:
            raise TypeError(op_err_msg(self, other))
    
    def append_coverage(self, coverage, pos='top'):
        """
        Append coverage to this track.

        Args:
            coverage (:obj:`Coverage`): coverage object to be added.
            pos (`str`, 'top' or 'bottom'): add coverage to top layer or bottom layer. ['top']
        """
        if pos == 'top':
            self.coverages.append(coverage)
        else:
            self.coverages.insert(0, coverage)

    def pile_coverages(self, coverages, pos='top'):
        """
        Pile a stack of coverages with self's coverages

        Args:
            coverages (:obj:`list` of `Coverage`): coverage objects to be piled.
            pos (`str`, 'top' or 'bottom'): add coverages to top or bottom. ['top']
        """
        if not hasattr(self, 'coverages'):
            self.coverages = coverages
        else:
            if pos == 'top':
                self.coverages.extend(coverages)
            else:
                self.coverages = coverages + self.coverages


class Spacer(Track, PlotSpacer):
    """
    The spacer track,
    not have any real content, just used to split two tracks.
    """

    DEFAULT_HEIGHT = 1

    def __init__(self, height=None, name=None):
        """
        Args:
            height (float, optional): the height of Spacer track. default: Spacer.DEFAULT_HEIGHT
            name (str, optional): Track's name
        """

        properties_dict = {}

        if height is None:
            height = Spacer.DEFAULT_HEIGHT

        properties_dict['height'] = height

        super().__init__(properties_dict, name)


class XAxis(Track, PlotXAxis):
    """
    The x axis track.
    """

    DEFAULT_FONTSIZE = 15
    DEFAULT_HEIGHT = 2

    def __init__(self, height=None, fontsize=None, where='buttom', name=None):
        """
        Args:
            height (float, optional): height of Spacer track. [XAxis.DEFAULT_HEIGHT]
            fontsize (int, optional): font size of XAxis. [XAxis.DEFAULT_FONTSIZE]
            where (str, optional): the position of tick labels relative to the axis.
                selections ('buttom', 'top'), ['buttom']
            name (str, optional): Track's name
        """

        properties_dict = {}

        if height is None:
            height = XAxis.DEFAULT_HEIGHT
        if fontsize is None:
            fontsize = XAxis.DEFAULT_FONTSIZE

        properties_dict['height'] = height
        properties_dict['fontsize'] = fontsize
        properties_dict['where'] = where

        super().__init__(properties_dict, name)


class Bed(Track, PlotBed, FetchBed):
    """
    Bed track.
    """

    DEFAULT_FONTSIZE = 12

    def __init__(self, file_, height=None, color='bed_rgb', border_color='black',
                 fontsize=None, title='', labels='auto', style='flybase', display='stacked',
                 interval_height=None, global_max_row=None, gene_rows=None, max_value=None, min_value=None,
                 name=None):
        """
        Parameters
        ----------
            file_ (str): path to bed file.
            height (float, optional): height of track. [Bed.DEFAULT_HEIGHT]
            fontsize (int, optional): font size. [Bed.DEFAULT_FONTSIZE]
            color (str, optional): track color,
                'bed_rgb' for auto specify color according to record. ['bed_rgb']
            border_color (str, optional): border_color of gene. ['black']
            title (str, optional): label text. ['']
            labels (str, optional): draw bed name or not.
                'on' or 'off' or 'auto'. ['auto']
            display (str, optional): display mode,
                options ('stacked', 'interlaced', 'collapsed'). ['stacked']
            interval_height (int, optional)
            global_max_row (int, optional)
            gene_row (int, optional)
            max_value (float, optional)
            min_value (float, optional)
            name (str, optional): Track's name
        """

        properties_dict = {}

        if height is None:
            height = Bed.DEFAULT_HEIGHT
        if fontsize is None:
            fontsize = Bed.DEFAULT_FONTSIZE

        properties_dict['file'] = file_
        properties_dict['height'] = height
        properties_dict['color'] = color
        properties_dict['border_color'] = border_color
        properties_dict['fontsize'] = fontsize
        properties_dict['title'] = title
        properties_dict['labels'] = labels
        properties_dict['style'] = style
        properties_dict['display'] = display
        if interval_height is not None:
            properties_dict['interval_height'] = interval_height
        if global_max_row is not None:
            properties_dict['global_max_row'] = global_max_row
        if gene_rows is not None:
            properties_dict['gene_rows'] = gene_rows
        if max_value is not None:
            properties_dict['max_value'] = max_value
        if min_value is not None:
            properties_dict['min_value'] = min_value

        super().__init__(properties_dict, name)


class BigWig(Track, PlotBigWig, FetchBigWig):
    """
    BigWig track.
    """

    DEFAULT_COLOR = "#dfccde"

    def __init__(self, file_, height=None, color=None,
                 number_of_bins=700, type_=None, orientation=None, show_data_range='yes',
                 title='', max_value='auto', min_value='auto', name=None):
        """
        Args:
            file_ (str): path to bigwig file.
            height (float, optional): height of track. [BigWig.DEFAULT_HEIGHT]
            fontsize (int, optional): font size. [BigWig.DEFAULT_FONTSIZE]
            color (str, optional): track color, [BigWig.DEFAULT_COLOR]
            number_of_bins (int, optional): number_of_bins in current range. [700]
            type_ (str, optional): track graph type, format 'type:size', like 'line:2', 'points:0.5' 
            orientation (str, optional): track orientation, use 'inverted' for inverted track plot.
            show_data_range (str, optional): show_data_range or not. 'yes' or 'no' ['yes']
            title (str, optional): label text. ['']
            max_value (`float` or `'auto'`, optional): max value of track. ['auto']
            min_value (`float` or `'auto'`, optional): min value of track. ['auto']
            name (str, optional): Track's name
        """

        properties_dict = {}

        if height is None:
            height = BigWig.DEFAULT_HEIGHT
        if color is None:
            color = BigWig.DEFAULT_COLOR

        properties_dict['file'] = file_
        properties_dict['height'] = height
        properties_dict['color'] = color
        properties_dict['number_of_bins'] = number_of_bins
        if type_ is not None:
            properties_dict['type'] = type_
        if orientation is not None:
            properties_dict['orientation'] = orientation
        properties_dict['show_data_range'] = show_data_range
        properties_dict['title'] = title
        properties_dict['max_value'] = max_value
        properties_dict['min_value'] = min_value

        super().__init__(properties_dict, name)


class ABCompartment(BigWig):
    """
    A/B Comapartment BigWig track.
    """

    DEFAULT_POSITIVE_COLOR = "#ff9c9c"
    DEFAULT_NEGATIVE_COLOR = "#66ccff"

    def __init__(self, file_, positive_color=None, negative_color=None, **kwargs):
        """
        Args:
            file_ (str): path to bigwig file.
            height (float, optional): height of track. [BigWig.DEFAULT_HEIGHT]
            fontsize (int, optional): font size. [BigWig.DEFAULT_FONTSIZE]
            positive_color (str, optional): track's positive value color, [ABCompartment.DEFAULT_POSITIVE_COLOR]
            negative_color (str, optional): track's negative value color, [ABCompartment.DEFAULT_NEGATIVE_COLOR]
            number_of_bins (int, optional): number_of_bins in current range. [700]
            type_ (str, optional): track graph type, format 'type:size', like 'line:2', 'points:0.5' 
            orientation (str, optional): track orientation, use 'inverted' for inverted track plot.
            show_data_range (str, optional): show_data_range or not. 'yes' or 'no' ['yes']
            title (str, optional): label text. ['']
            max_value (`float` or `'auto'`, optional): max value of track. ['auto']
            min_value (`float` or `'auto'`, optional): min value of track. ['auto']
            name (str, optional): Track's name
        """
        super().__init__(file_, *kwargs)
        if positive_color is None:
            self.properties['positive_color'] = ABCompartment.DEFAULT_POSITIVE_COLOR
        if negative_color is None:
            self.properties['negative_color'] = ABCompartment.DEFAULT_NEGATIVE_COLOR


class BedGraph(Track, PlotBedGraph, FetchBedGraph):
    """
    BedGraph track.
    """

    DEFAULT_COLOR = '#a6cee3'

    def __init__(self, file_, height=None, color=None,
                 extra=None, show_data_range='yes', title='',
                 max_value='auto', min_value='auto', name=None):
        """
        Args:
            file_ (str): path to bedgraph file.
            height (float, optional): height of track. [BigWig.DEFAULT_HEIGHT]
            color (str, optional): track color, [BigWig.DEFAULT_COLOR]
            type_ (str, optional): track graph type, format 'type:size', like 'line:2', 'points:0.5' 
            extra (optional):
            show_data_range (str, optional): show_data_range or not. 'yes' or 'no' ['yes']
            title (str, optional): label text. ['']
            max_value (`float` or `'auto'`, optional): max value of track. ['auto']
            min_value (`float` or `'auto'`, optional): min value of track. ['auto']
            name (str, optional): Track's name
        """

        properties_dict = {}

        if height is None:
            height = BedGraph.DEFAULT_HEIGHT
        if color is None:
            color = BedGraph.DEFAULT_COLOR

        properties_dict['file'] = file_
        properties_dict['height'] = height
        properties_dict['color'] = color
        if extra:
            properties_dict['extra'] = extra
        properties_dict['show_data_range'] = show_data_range
        properties_dict['title'] = title
        properties_dict['max_value'] = max_value
        properties_dict['min_value'] = min_value

        super().__init__(properties_dict, name)


class Boundaries(Track, PlotBoundaries, FetchBed):
    """
    Boundaries track.
    """

    def __init__(self, file_, height=None, name=None):
        """
        Args:
            file_ (str): path to bedgraph file.
            height (float): height of track. [Boundaries.DEFAULT_HEIGHT]
            name (str, optional): Track's name
        """

        properties_dict = {}

        if height is None:
            height = Boundaries.DEFAULT_HEIGHT

        properties_dict['file'] = file_
        properties_dict['height'] = height

        super().__init__(properties_dict, name)


class Arcs(Track, PlotArcs, FetchArcs):
    """
    Arcs(link) track.
    """

    DEFAULT_COLOR = '#3297dc'
    DEFAULT_ALPHA = 0.8

    def __init__(self, file_, height=None, color=None, alpha=0.8,
                 line_width=None, orientation=None, title='', name=None):
        """
        Args:
            file_ (str): path to bedgraph file.
            height (float): height of track. [Boundaries.DEFAULT_HEIGHT]
            color (str, optional): track color, [BigWig.DEFAULT_COLOR]
            alpha (float, optional): alpha value of track. [0.8]
            line_width (float, optional): width of arc line.
            orientation (str, optional): track orientation, use 'inverted' for inverted track plot.
            title (str, optional): label text. ['']
            name (str, optional): Track's name
        """

        properties_dict = {}

        if height is None:
            height = Arcs.DEFAULT_HEIGHT
        if color is None:
            color = Arcs.DEFAULT_COLOR
        if alpha is None:
            alpha = Arcs.DEFAULT_ALPHA

        properties_dict['file'] = file_
        properties_dict['height'] = height
        properties_dict['color'] = color
        if line_width:
            properties_dict['line_width'] = line_width
        if orientation:
            properties_dict['orientation'] = orientation
        properties_dict['title'] = title

        super().__init__(properties_dict, name)


class TADs(Track, PlotTADs, FetchBed):
    """
    TADs track.
    """

    def __init__(self, file_, height=None, color='bed_rgb', border_color='black',
                 orientation=None, title='', name=None):
        """
        Args:
            file_ (str): path to bed file.
            height (float, optional): height of track. [TADs.DEFAULT_HEIGHT]
            fontsize (int, optional): font size. [TADs.DEFAULT_FONTSIZE]
            color (str, optional): track color, use 'bed_rgb' for specify color according to file. ['bed_rgb']
            border_color (str, optional): border_color of gene. ['black']
            orientation (str, optional): track orientation, use 'inverted' for inverted track plot.
            title (str, optional): label text. ['']
            name (str, optional): Track's name
        """

        properties_dict = {}

        if height is None:
            height = TADs.DEFAULT_HEIGHT
        if color is None:
            color = TADs.DEFAULT_COLOR

        properties_dict['file'] = file_
        properties_dict['height'] = height
        properties_dict['color'] = color
        properties_dict['border_color'] = border_color
        if orientation:
            properties_dict['orientation'] = orientation
        properties_dict['title'] = title

        super().__init__(properties_dict, name)


class Cool(Track, PlotCool, FetchCool):
    """
    Cool Hi-C matrix (or triangular matrix) track.
    """

    DEFAULT_COLOR = 'YlOrRd'

    def __init__(self, file_, cmap=None, triangular=True, balance=True,
                 depth_ratio='full', color_bar=True, transform='no',
                 orientation="normal",
                 max_value='auto', min_value='auto', title='',
                 name=None):
        """
        Args:
            file_ (str): path to bed file.
            cmap (str, optional): color map of hic matrix. [Cool.DEFAULT_COLOR]
            triangular (bool, optional): show traiangular form matrix or not. [True]
            balance (bool, optional): show balanced matrix or not. [True]
            depth_ratio (float, optional): depth ratio of triangular matrix, use 'full' for full depth. ['full']
            color_bar (bool, optional): show color_bar or not. [True]
            transform (str, optional): transform for matrix, like 'log2', 'log10', use 'no' for not transform. ['no']
            orientation (str, optional): track orientation, use 'inverted' for inverted track plot.
            max_value (`float` or `'auto'`, optional): max value of hic matrix. ['auto']
            min_value (`float` or `'auto'`, optional): min value of hic matrix. ['auto']
            title (str, optional): label text. ['']
            name (str, optional): Track's name
        """

        properties_dict = {}

        if cmap is None:
            cmap = Cool.DEFAULT_COLOR

        properties_dict['file'] = file_
        properties_dict['color'] = cmap
        properties_dict['triangular'] = 'yes' if triangular else 'no'
        properties_dict['balance'] = 'yes' if balance else 'no'
        properties_dict['depth_ratio'] = depth_ratio
        properties_dict['color_bar'] = 'yes' if color_bar else 'no'
        properties_dict['transform'] = transform
        properties_dict['orientation'] = orientation
        properties_dict['max_value'] = max_value
        properties_dict['min_value'] = min_value
        properties_dict['title'] = title

        super().__init__(properties_dict, name)


class Coverage(object):
    """
    Coverage base class.

    `Coverage` is the plots at the upper layer of Track.

    >>> c1 = Coverage({})
    >>> c1.properties
    {}
    >>> c2 = Coverage({})
    >>> t1 = Track({})
    >>> t2 = c1 + c2 + t1
    >>> len(t2.coverages)
    2
    >>> assert type(c1 + c2) is CoverageStack
    >>> t3 = Track({})
    >>> frame = t2 + t3
    >>> frame = frame + Coverage({})
    >>> len(list(frame.tracks.values())[-1].coverages)
    1
    """

    def __init__(self, properties_dict):
        self.properties = properties_dict
        super().__init__()

        scope = globals()
        stack = scope[FEATURES_STACK_NAME]
        for feature in stack:
            self.properties[feature.key] = feature.value

    def __add__(self, other):
        if isinstance(other, Track):
            result = copy(other)
            result.append_coverage(self)
            return result
        elif isinstance(other, Frame):
            result = copy(other)
            if len(result.tracks) > 1:
                first = list(result.tracks.values())[0]
                first.append_coverage(self, pos='bottom')
            return result
        elif isinstance(other, Feature):
            result = copy(self)
            result[other.key] = other.value
            return result
        elif isinstance(other, Coverage):
            stack = CoverageStack([self, other])
            return stack
        elif isinstance(other, CoverageStack):
            result = copy(other)
            result.to_bottom(self)
            return result
        else:
            raise TypeError(op_err_msg(self, other))

    def __mul__(self, other):
        if isinstance(other, Frame):
            result = copy(other)
            result.add_cov_to_tracks(self)
            return result
        else:
            raise TypeError(op_err_msg(self, other, op='*'))

    def __enter__(self):
        scope = globals()
        stack = scope[COVERAGE_STACK_NAME]
        stack.append(self)
        return self

    def __exit__(self, type, value, traceback):
        scope = globals()
        stack = scope[COVERAGE_STACK_NAME]
        stack.pop()


class CoverageStack(object):
    """
    Denote a stack of Coverage.

    ps: this "Stack" is actually a "Deque",
        name it "Stack" is just for imaging it vertically. 
    """

    def __init__(self, coverages):
        """
        Args:
            coverages (:obj:`list` of :obj:`Coverage`): coverages list.
        """
        self.coverages = coverages

    def to_top(self, cov):
        self.coverages.append(cov)

    def to_bottom(self, cov):
        self.coverages.insert(0, cov)

    def __add__(self, other):
        if isinstance(other, Coverage):
            result = copy(self)
            result.to_top(other)
            return result
        elif isinstance(other, Track):
            result = copy(other)
            result.pile_coverages(self.coverages, pos='bottom')
            return result
        elif isinstance(other, Frame):
            result = copy(other)
            if len(result.tracks) != 0:
                first = list(result.tracks.values())[0]
                first.pile_coverages(self.coverages, pos='bottom')
            return result
        elif isinstance(other, Feature):
            result = copy(self)
            if len(result.coverages) != 0:
                last = result.coverages[-1]
                last.properties[other.key] = other.value
            return result
        else:
            raise TypeError(op_err_msg(self, other))


class VlinesFromFile(Coverage, PlotVlines):
 
    def __init__(self, file_, color='#1e1e1e', alpha=0.8,
                 line_style='dashed', line_width=1):
        """
        Args:
            file_ (str): 
            color (str, optional): ['#1e1e1e']
            alpha (float, optional): [0.8]
            line_style (str, optional): ['dashed']
            line_width (float, optional): [0.5]
        """

        properties_dict = dict()

        properties_dict['file'] = file_
        properties_dict['color'] = color
        properties_dict['alpha'] = 0.8
        properties_dict['line_style'] = line_style
        properties_dict['line_width'] = line_width

        super().__init__(properties_dict)


class Vlines(Coverage, PlotVlines):

    def __init__(self, vlines, chr=None, color='#1e1e1e', alpha=0.8,
                 line_style='dashed', line_width=1):
        """
        Args:
            vlines (:obj:`list` of `int`): A list of vline positions.
            chr (str, optional): chromosome of vline, if not specify will plot in all chromosome.
            color (str, optional): ['#1e1e1e']
            alpha (float, optional): [0.8]
            line_style (str, optional): ['dashed']
            line_width (float, optional): [0.5]
        """

        properties_dict = dict()

        properties_dict['vlines_list'] = vlines
        properties_dict['color'] = color
        properties_dict['alpha'] = 0.8
        properties_dict['line_style'] = line_style
        properties_dict['line_width'] = line_width
        properties_dict['chr'] = chr

        super().__init__(properties_dict)


class HighLightsFromFile(Coverage, PlotHighLightRegions):

    def __init__(self, file_, color='bed_rgb', alpha=0.6,
                 border_line='yes', border_line_style='dashed',
                 border_line_width=1.0, border_line_color='#000000',
                 border_line_alpha=0.8):
        """
        Args:
            file_ (str):
            color (str, optional): ['bed_rgb']
            alpha (float, optional): [0.6]
            border_line (str, optional): plot border line or not. ['yes']
            border_line_style (str, optional): ['dashed']
            border_line_width (float, optional): [1.0]
            border_line_color (str, optional): ['#000000']
            border_line_alpha (float, optional): [0.8]
        """

        properties_dict = dict()

        properties_dict['file'] = file_
        properties_dict['color'] = color
        properties_dict['alpha'] = alpha
        properties_dict['border_line'] = border_line
        properties_dict['border_line_style'] = border_line_style
        properties_dict['border_line_width'] = border_line_width
        properties_dict['border_line_color'] = border_line_color
        properties_dict['border_line_alpha'] = border_line_alpha

        super().__init__(properties_dict)


class HighLights(Coverage, PlotHighLightRegions):

    DEFAULT_COLOR = "#ff9c9c"

    def __init__(self, highlight_regions, chr=None, color=None, alpha=0.6, border_line='yes',
                 border_line_style='dashed', border_line_width=1.0,
                 border_line_color='#000000', border_line_alpha=0.8):
        """
        Args:
            highlight_regions (:obj:`list` of :obj:`tuple`): A list of regions for highlights,
                region tuple format: `(start, end)` like, [(100000, 120000), (130000, 150000)].
            chr (str, optional): chromosome of highlight regions, if not specify will plot in all chromosome.
            color (str, optional): [HighLights.DEFAULT_COLOR]
            alpha (float, optional): [0.6]
            border_line (str, optional): plot border line or not. ['yes']
            border_line_style (str, optional): ['dashed']
            border_line_width (float, optional): [1.0]
            border_line_color (str, optional): ['#000000']
            border_line_alpha (float, optional): [0.8]
        """

        if color is None:
            color = HighLights.DEFAULT_COLOR

        properties_dict = dict()

        properties_dict['highlight_regions'] = highlight_regions
        properties_dict['color'] = color
        properties_dict['alpha'] = alpha
        properties_dict['border_line'] = border_line
        properties_dict['border_line_style'] = border_line_style
        properties_dict['border_line_width'] = border_line_width
        properties_dict['border_line_color'] = border_line_color
        properties_dict['border_line_alpha'] = border_line_alpha
        properties_dict['chr'] = chr

        super().__init__(properties_dict)


class HiCPeaks(Coverage, PlotHiCPeaks):
    """
    Hi-C Peaks(Loops) Coverge is a special kind of Coverage.
    Used to show the peaks on the Hi-C interaction map.
    """

    def __init__(self, file_, color='bed_rgb', alpha=0.8,
                 line_width=1.5, line_style='solid',
                 fill='no', fill_color='bed_rgb'):
        """
        Args:
            file_ (str): path to the loop file, loop file is a tab splited text file, fields:
                chr1, x1, x2, chr2, y1, y2, [color], ... (other optional fields)
            color (str, optional): ['bed_rgb']
            alpha (float, optional): [0.8]
            line_width (float, optional): [1.0]
            line_style (str, optional): ['solid']
            fill (str, optional): ['no']
            fill_color (str, optional): ['bed_rgb']
        """

        properties_dict = {}

        properties_dict['file'] = file_
        properties_dict['color'] = color
        properties_dict['alpha'] = alpha
        properties_dict['line_width'] = line_width
        properties_dict['line_style'] = line_style
        properties_dict['fill'] = fill
        properties_dict['fill_color'] = fill_color

        super().__init__(properties_dict)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
