from collections import deque
from copy import copy

from coolbox.plots.track import *
from coolbox.fetchdata import (
    FetchBed, FetchBedGraph, FetchBigWig,
    FetchArcs, FetchCool
)
from coolbox.utilities import op_err_msg

import coolbox.api.frame
import coolbox.api.feature
import coolbox.api.coverage

Frame = coolbox.api.frame.Frame
Feature = coolbox.api.feature.Feature
FrameFeature = coolbox.api.feature.FrameFeature
Coverage = coolbox.api.coverage.Coverage
CoverageStack = coolbox.api.coverage.CoverageStack


FEATURES_STACK_NAME = "__COOLBOX_FEATURE_STACK__"
COVERAGE_STACK_NAME = "__COOLBOX_COVERAGE_STACK__"
global_scope = globals()
global_scope[FEATURES_STACK_NAME] = deque()
global_scope[COVERAGE_STACK_NAME] = deque()


__all__ = [
    "Spacer", "XAxis", "Bed", "TADs",
    "BigWig", "ABCompartment", "BedGraph",
    "Arcs", "Cool"
]


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


if __name__ == "__main__":
    import doctest
    doctest.testmod()
