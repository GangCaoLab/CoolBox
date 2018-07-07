from collections import deque
from copy import copy

from coolbox.plots.coverage import *
from coolbox.utilities import op_err_msg


FEATURES_STACK_NAME = "__COOLBOX_FEATURE_STACK__"
COVERAGE_STACK_NAME = "__COOLBOX_COVERAGE_STACK__"
global_scope = globals()
global_scope[FEATURES_STACK_NAME] = deque()
global_scope[COVERAGE_STACK_NAME] = deque()


__all__ = [
    "CoverageStack",
    "Vlines", "VlinesFromFile",
    "HighLights", "HighLightsFromFile",
    "HiCPeaks"
]


class Coverage(object):
    """
    Coverage base class.

    `Coverage` is the plots at the upper layer of Track.

    >>> from coolbox.core.track import Track
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
        from .track import Track
        from .frame import Frame
        from .feature import Feature

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
        from .frame import Frame
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

    Parameters
    ----------
    coverages : list of coolbox.core.coverage.Coverage
        coverages list.
    """

    def __init__(self, coverages):
        self.coverages = coverages

    def to_top(self, cov):
        self.coverages.append(cov)

    def to_bottom(self, cov):
        self.coverages.insert(0, cov)

    def __add__(self, other):
        from .track import Track
        from .frame import Frame
        from .feature import Feature

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

    """
    Vertical lines from the file.

    Parameters
    ----------
    file_ : str
        Path to file.

    color : str, optional
        Line color, default '#1e1e1e'.

    alpha : float, optional
        Line alpha value, default 0.8.

    line_style : str, optional
        Line style, default 'dashed'.

    line_width : (float, optional)
        Line width, default 0.5.
    """

    def __init__(self, file_, color='#1e1e1e', alpha=0.8,
                 line_style='dashed', line_width=1):
        properties_dict = dict()

        properties_dict['file'] = file_
        properties_dict['color'] = color
        properties_dict['alpha'] = 0.8
        properties_dict['line_style'] = line_style
        properties_dict['line_width'] = line_width

        super().__init__(properties_dict)


class Vlines(Coverage, PlotVlines):

    """
    Vertical lines.

    Parameters
    ----------
    vlines : list of int
        A list of vline positions.

    chr : str, optional
        Chromosome of vline, if not specify will plot in all chromosome.

    color : str, optional
        Line color, default '#1e1e1e'.

    alpha : float, optional
        Line alpha value, default 0.8.

    line_style : str, optional
        Line style, default 'dashed'.

    line_width : float, optional
        Line width, default 0.5.
    """

    def __init__(self, vlines, chr=None, color='#1e1e1e', alpha=0.8,
                 line_style='dashed', line_width=1):
        properties_dict = dict()

        properties_dict['vlines_list'] = vlines
        properties_dict['color'] = color
        properties_dict['alpha'] = 0.8
        properties_dict['line_style'] = line_style
        properties_dict['line_width'] = line_width
        properties_dict['chr'] = chr

        super().__init__(properties_dict)


class HighLightsFromFile(Coverage, PlotHighLightRegions):

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
        High light region alpha value, default 0.6.

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
    """

    def __init__(self, file_, color='bed_rgb', alpha=0.6,
                 border_line=True, border_line_style='dashed',
                 border_line_width=1.0, border_line_color='#000000',
                 border_line_alpha=0.8):
        properties_dict = dict()

        if border_line:
            properties_dict['border_line'] = 'yes'
        else:
            properties_dict['border_line'] = 'no'

        properties_dict['file'] = file_
        properties_dict['color'] = color
        properties_dict['alpha'] = alpha
        properties_dict['border_line_style'] = border_line_style
        properties_dict['border_line_width'] = border_line_width
        properties_dict['border_line_color'] = border_line_color
        properties_dict['border_line_alpha'] = border_line_alpha

        super().__init__(properties_dict)


class HighLights(Coverage, PlotHighLightRegions):

    """
    High light region.

    Parameters
    ----------
    highlight_regions : list of tuple
        A list of regions for highlights, region tuple format: `(start, end)`,
        like, [(100000, 120000), (130000, 150000)].

    chr : str, optional
        Chromosome of highlight regions, if not specify will plot in all chromosome.

    color : str, optional
        High light region color, default HighLights.DEFAULT_COLOR.

    alpha : float, optional
        High light region alpha value, default 0.6

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
    """

    DEFAULT_COLOR = "#ff9c9c"

    def __init__(self, highlight_regions, chr=None, color=None, alpha=0.6, border_line='yes',
                 border_line_style='dashed', border_line_width=1.0,
                 border_line_color='#000000', border_line_alpha=0.8):

        if color is None:
            color = HighLights.DEFAULT_COLOR

        properties_dict = dict()

        if border_line:
            properties_dict['border_line'] = 'yes'
        else:
            properties_dict['border_line'] = 'no'

        properties_dict['highlight_regions'] = highlight_regions
        properties_dict['color'] = color
        properties_dict['alpha'] = alpha
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

    Parameters
    ----------
    file_ : str
        Path to the loop file, loop file is a tab splited text file, fields:
        chr1, x1, x2, chr2, y1, y2, [color], ... (other optional fields)

    color : str, optional
        Peak color, use 'bed_rgb' for specify color in file,
        default 'bed_rgb'.

    alpha : float, optional
        Peak alpha value, default 0.8.

    line_width : float, optional
        Peak border line width, default 1.0

    line_style : str, optional
        Border line style, default 'solid'

    fill : bool, optional
        Fill center or not, default False.

    fill_color : str, optional
        Fill color, use 'bed_rgb' for specify color in file,
        default 'bed_rgb'.
    """

    def __init__(self, file_, color='bed_rgb', alpha=0.8,
                 line_width=1.5, line_style='solid',
                 fill='no', fill_color='bed_rgb'):
        properties_dict = {}

        if fill:
            properties_dict['fill'] = 'yes'
        else:
            properties_dict['fill'] = 'no'

        properties_dict['file'] = file_
        properties_dict['color'] = color
        properties_dict['alpha'] = alpha
        properties_dict['line_width'] = line_width
        properties_dict['line_style'] = line_style
        properties_dict['fill_color'] = fill_color

        super().__init__(properties_dict)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
