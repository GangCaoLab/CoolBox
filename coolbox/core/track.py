from collections import deque
from copy import copy

from coolbox.plots.track import *
from coolbox.fetchdata import (
    FetchBed, FetchBedGraph, FetchBigWig,
    FetchArcs, FetchCool
)
from coolbox.utilities import op_err_msg


FEATURES_STACK_NAME = "__COOLBOX_FEATURE_STACK__"
COVERAGE_STACK_NAME = "__COOLBOX_COVERAGE_STACK__"
global_scope = globals()
global_scope[FEATURES_STACK_NAME] = deque()
global_scope[COVERAGE_STACK_NAME] = deque()


__all__ = [
    "Spacer", "XAxis", "Bed", "TADs",
    "BigWig", "ABCompartment", "BedGraph",
    "Arcs", "Cool", "HicCompare"
]


class Track(object):
    """
    Track base class.

    Parameters
    ----------
    properties_dict : dict
        The properties(features) of this track. For example 'height', 'color'...

    name : str, optional
        The name of Track.
        (Default: "{self.__class__.__name__}.{self.__class__._counts}")


    Attributes
    ----------
    properties : dict
        The properties(features) of this track. For example 'height', 'color'...

    name : str
        The name of Track.

    coverages : list of `coolbox.api.coverage.Coverage`
        Coverages on this Track.
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

        # load features from global feature stack
        scope = globals()
        features_stack = scope[FEATURES_STACK_NAME]
        for features in features_stack:
            self.properties[features.key] = features.value

        # load coverages from global coverages stack
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
        from .frame import Frame
        from .coverage import Coverage
        from .coverage import CoverageStack
        from .feature import Feature

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

        Parameters
        ----------
        coverage : `coolbox.api.coverage.Coverage`
            Coverage object to be piled.

        pos : {'top', 'bottom'}, optional
            Add coverages to top or bottom. (Default: 'top')
        """
        if pos == 'top':
            self.coverages.append(coverage)
        else:
            self.coverages.insert(0, coverage)

    def pile_coverages(self, coverages, pos='top'):
        """
        Pile a stack of coverages with self's coverages

        Parameters
        ----------
        coverages : list of `coolbox.api.coverage.Coverage` or `coolbox.api.coverage.CoverageStack`
            Coverage objects to be piled.

        pos : {'top', 'bottom'}, optional
            Add coverages to top or bottom. (Default: 'top')
        """
        from .coverage import CoverageStack

        if isinstance(coverages, CoverageStack):
            coverages = coverages.coverages
        elif isinstance(coverages, list):
            pass
        else:
            raise TypeError("coverages must a list of Coverage or CoverageStack")

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

    Parameters
    ----------
    height : float, optional
        The height of Spacer track. (Default: Spacer.DEFAULT_HEIGHT)

    name : str, optional
        Track's name. (Default: "Spacer.{Spacer._counts}")
    """

    DEFAULT_HEIGHT = 1

    def __init__(self, height=None, name=None):
        properties_dict = {}

        if height is None:
            height = Spacer.DEFAULT_HEIGHT

        properties_dict['height'] = height

        super().__init__(properties_dict, name)


class XAxis(Track, PlotXAxis):
    """
    The x axis track.

    Parameters
    ----------
    height : float, optional
        Height of Spacer track. (Default: XAxis.DEFAULT_HEIGHT)

    fontsize : int, optional
        Font size of XAxis. (Default: XAxis.DEFAULT_FONTSIZE)

    where : {'top', 'bottom'}, optional
        The position of tick labels relative to the axis.
        (Default: 'bottom')

    name (str, optional):
        Track's name. (Default: "XAxis.{XAxis._counts}")
    """

    DEFAULT_FONTSIZE = 15
    DEFAULT_HEIGHT = 2

    def __init__(self, height=None, fontsize=None, where='bottom', name=None):
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

    Parameters
    ----------
    file_ : str
        Path to bed file.

    height : float, optional
        Height of track. (Default: Bed.DEFAULT_HEIGHT)

    fontsize : int, optional
        Font size. (Default: Bed.DEFAULT_FONTSIZE)

    color : str, optional
        Track color, 'bed_rgb' for auto specify color according to bed record.
        (Default: 'bed_rgb')

    border_color : str, optional
        Border_color of gene. (Default: 'black')

    title : str, optional
        Label text. (Default: '')

    labels : {True, False, 'auto'}, optional
        Draw bed name or not. 'auto' for automate decision according to density.
        'on' or 'off' or 'auto'. (Default: 'auto')

    display : {'stacked', 'interlaced', 'collapsed'}, optional
        Display mode. (Default: 'stacked')

    interval_height : int, optional
        The height of the interval. (Default: 100)

    global_max_row : bool, optional
        If set to True, will process the whole bed regions
        at the given figure length and font size to
        determine the maximum number of rows required. (Default: False)

    gene_row : int, optional
        Set the max interval rows. (Default: unlimited interval rows)

    max_value : float, optional
        Max score. (Default: inf)

    min_value : float, optional
        Min score. (Default: -inf)

    name : str, optional
        Track's name (Default: "Bed.{Bed._counts}"

    """

    DEFAULT_FONTSIZE = 12

    def __init__(self, file_, height=None, color='bed_rgb', border_color='black',
                 fontsize=None, title='', labels='auto', style='flybase', display='stacked',
                 interval_height=None, global_max_row=None, gene_rows=None, max_value=None, min_value=None,
                 name=None):
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
        properties_dict['style'] = style
        properties_dict['display'] = display

        if labels == 'auto':
            properties_dict['labels'] = 'auto'
        elif labels is True:
            properties_dict['labels'] = 'on'
        else:
            properties_dict['labels'] = 'off'

        if interval_height is not None:
            properties_dict['interval_height'] = interval_height
        if global_max_row is not None:
            properties_dict['global_max_row'] = "yes" if global_max_row else "no"
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

    Parameters
    ----------
    file_ : str
        Path to bigwig file.

    height : float, optional
        Height of track, default BigWig.DEFAULT_HEIGHT.

    fontsize : int, optional
        Font size, default BigWig.DEFAULT_FONTSIZE.

    color : str, optional
        Track color, default BigWig.DEFAULT_COLOR.

    number_of_bins : int, optional
        Number_of_bins in current range, default 700.

    type_ : str, optional
        Track graph type, format 'type:size', like 'line:2', 'points:0.5'

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    show_data_range : bool, optional
        Show_data_range or not, default True.

    title : str, optional
        Label text, default ''

    max_value : {float, 'auto'}, optional
        Max value of track. 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of track. 'auto' for specify max value automatically, default 'auto'.

    name : str, optional
        Track's name.
    """

    DEFAULT_COLOR = "#dfccde"

    def __init__(self, file_, height=None, color=None,
                 number_of_bins=700, type_=None, orientation=None, show_data_range='yes',
                 title='', max_value='auto', min_value='auto', name=None):

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
        properties_dict['show_data_range'] = 'yes' if show_data_range else 'no'
        properties_dict['title'] = title
        properties_dict['max_value'] = max_value
        properties_dict['min_value'] = min_value

        super().__init__(properties_dict, name)


class ABCompartment(BigWig):
    """
    A/B Comapartment BigWig track.

    Parameters
    ----------
    file_ : str
        Path to bigwig file.

    height : float, optional
        Height of track, default BigWig.DEFAULT_HEIGHT.

    fontsize : int, optional
        Font size, default BigWig.DEFAULT_FONTSIZE

    positive_color : str, optional
        Track's positive value color, default ABCompartment.DEFAULT_POSITIVE_COLOR

    negative_color : str, optional
        Track's negative value color, default ABCompartment.DEFAULT_NEGATIVE_COLOR

    number_of_bins : int, optional
        Number_of_bins in current range, default 700

    type_ : str, optional
        Track graph type, format 'type:size', like 'line:2', 'points:0.5'

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    show_data_range : bool, optional
        Show_data_range or not, default True.

    title : str, optional
        Label text. default ''

    max_value : {float, 'auto'}, optional
        Max value of track, use 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of track, use 'auto' for specify min value automatically, default 'auto'.

    name : str, optional
        Track's name
    """

    DEFAULT_POSITIVE_COLOR = "#ff9c9c"
    DEFAULT_NEGATIVE_COLOR = "#66ccff"

    def __init__(self, file_, positive_color=None, negative_color=None, **kwargs):
        super().__init__(file_, *kwargs)
        if positive_color is None:
            self.properties['positive_color'] = ABCompartment.DEFAULT_POSITIVE_COLOR
        if negative_color is None:
            self.properties['negative_color'] = ABCompartment.DEFAULT_NEGATIVE_COLOR


class BedGraph(Track, PlotBedGraph, FetchBedGraph):
    """
    BedGraph track.

    Parameters
    ----------
    file_ : str
        Path to bedgraph file.

    height : float, optional
        Height of track, default BigWig.DEFAULT_HEIGHT

    color : str, optional
        Track color, default BigWig.DEFAULT_COLOR

    type_ : str, optional
        Track graph type, format 'type:size', like 'line:2', 'points:0.5'

    extra : optional

    show_data_range : bool, optional
        Show_data_range or not, default True.

    title : str, optional
        Label text, default ''.

    max_value : {float, 'auto'}, optional
        Max value of track, use 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of track, use 'auto' for specify min value automatically, default 'auto'.

    name : str, optional
        Track's name.
    """

    DEFAULT_COLOR = '#a6cee3'

    def __init__(self, file_, height=None, color=None,
                 extra=None, show_data_range='yes', title='',
                 max_value='auto', min_value='auto', name=None):

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
        properties_dict['show_data_range'] = 'yes' if show_data_range else 'no'
        properties_dict['title'] = title
        properties_dict['max_value'] = max_value
        properties_dict['min_value'] = min_value

        super().__init__(properties_dict, name)


class Arcs(Track, PlotArcs, FetchArcs):
    """
    Arcs(link) track.

    Parameters
    ----------
    file_ : str
        Path to bedgraph file.

    height : float
        Height of track, default Boundaries.DEFAULT_HEIGHT.

    color : str, optional
        Track color, default BigWig.DEFAULT_COLOR.

    alpha : float, optional
        Alpha value of track, default 0.8.

    line_width : float, optional
        Width of arc line.

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    title : str, optional
        Label text. default ''

    name : str, optional
        Track's name.
    """

    DEFAULT_COLOR = '#3297dc'
    DEFAULT_ALPHA = 0.8

    def __init__(self, file_, height=None, color=None, alpha=0.8,
                 line_width=None, orientation=None, title='', name=None):

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

    def __init__(self, file_, height=None, color='bed_rgb', border_color='black',
                 orientation=None, title='', name=None):

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

    Parameters
    ----------
    file_ : str
        Path to bed file.

    cmap : str, optional
        Color map of hic matrix, default Cool.DEFAULT_COLOR.

    triangular : bool, optional
        Show traiangular form matrix or not, default True.

    balance : bool, optional
        Show balanced matrix or not, default True

    depth_ratio : float, optional
        Depth ratio of triangular matrix, use 'full' for full depth. default 'full'.

    color_bar : bool, optional
        Show color_bar or not, default True.

    transform : str, optional
        Transform for matrix, like 'log2', 'log10', use 'no' for not transform, default 'no'.

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    title : str, optional
        Label text, default ''.

    max_value : {float, 'auto'}, optional
        Max value of hic matrix, use 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of hic matrix, use 'auto' for specify min value automatically, default 'auto'.

    name : str, optional
        Track's name.
    """

    DEFAULT_COLOR = 'YlOrRd'

    def __init__(self, file_, cmap=None, triangular=True, balance=True,
                 depth_ratio='full', color_bar=True, transform='no',
                 orientation="normal",
                 max_value='auto', min_value='auto', title='',
                 name=None):

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


class HicCompare(Track, PlotHicCompare):
    """
    Track for express the comparison between two HiC Track.

    Parameters
    ----------
    hic1 : coolbox.api.track.Cool
        First HiC Track.

    hic2 : coolbox.api.track.Cool
        Second HiC Track.

    cmap : {str, matplotlib.colors.Colormap}, optional
        A diverging colormap, left color represent the first HiC file,
        and right represent the second HiC file.

    color_bar : bool, optional
        Show color bar or not.

    name : str, optional
        Track's name (Default: "Bed.{Bed._counts}"

    """

    DEFAULT_COLOR = 'bwr'

    def __init__(self, hic1, hic2,
                 cmap=None,
                 color_bar=True,
                 title='', name=None):
        properties_dict = {}

        if cmap is None:
            properties_dict['color'] = HicCompare.DEFAULT_COLOR

        properties_dict['hic1'] = hic1
        properties_dict['hic2'] = hic2
        if color_bar:
            properties_dict['color_bar'] = 'yes'
        else:
            properties_dict['color_bar'] = 'no'
        properties_dict['title'] = title

        super().__init__(properties_dict, name)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
