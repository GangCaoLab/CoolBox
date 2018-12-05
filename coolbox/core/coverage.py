from copy import copy

from coolbox.plots.coverage import *
from coolbox.utilities import op_err_msg, get_coverage_stack, get_feature_stack


__all__ = [
    "CoverageStack",
    "Vlines", "VlinesFromFile",
    "HighLights", "HighLightsFromFile",
    "HiCPeaks", "TADCoverage"
]


class Coverage(object):
    """
    Coverage base class.

    `Coverage` is the plots at the upper layer of Track.

    >>> from coolbox.core.track import Track
    >>> c1 = Coverage({})
    >>> c1.properties
    {'name': 'Coverage.1'}
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

    def __new__(cls, *args, **kwargs):
        if hasattr(cls, "_counts"):
            cls._counts += 1
        else:
            cls._counts = 1
        return super().__new__(cls)

    def __init__(self, properties_dict):
        self.properties = properties_dict
        self.__bool2str()
        name = self.properties.get("name")
        if name is not None:
            assert isinstance(name, str), "Coverage name must be a `str`."
        else:
            name = self.__class__.__name__ + ".{}".format(self.__class__._counts)
        self.properties['name'] = name

        super().__init__()

        stack = get_feature_stack()
        for feature in stack:
            self.properties[feature.key] = feature.value

    def __bool2str(self):
        """
        Conver bool value to 'yes' or 'no', for compatible with pyGenomeTracks
        """
        for key, value in self.properties.items():
            if isinstance(value, bool):
                if value:
                    self.properties[key] = 'yes'
                else:
                    self.properties[key] = 'no'

    @property
    def name(self):
        return self.properties['name']

    @name.setter
    def name(self, value):
        self.properties['name'] = value

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
            result.properties[other.key] = other.value
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
        stack = get_coverage_stack()
        stack.append(self)
        return self

    def __exit__(self, type, value, traceback):
        stack = get_coverage_stack()
        stack.pop()

    def check_track_type(self, allow):
        valid = any([isinstance(self.track, type_) for type_ in allow])
        if not valid:
            msg = "{} coverage's track must be a instance of {}".format(self.track.__class__.__name__,
                                                                        [type_.__name__ for type_ in allow])
            raise ValueError(msg)


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

    line_width : float, optional
        Line width, default 0.5.

    name : str, optional
        The name of thr Coverage.
    """

    def __init__(self, file_, **kwargs):
        properties_dict = {
            "file": file_,
            "color": "#1e1e1e",
            "alpha": 0.8,
            "line_style": "dashed",
            "line_width": 1,
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)


class Vlines(Coverage, PlotVlines):

    """
    Vertical lines.

    Parameters
    ----------
    vlines : list of {int, str}
        A list of vline positions. position can be expressed as a tuple like:
        [('chr1', 10000), ('chr2', 50000)]
        or a genome range string like:
        ['chr1:10000-10000', 'chr2:50000-50000']

    color : str, optional
        Line color, default '#1e1e1e'.

    alpha : float, optional
        Line alpha value, default 0.8.

    line_style : str, optional
        Line style, default 'dashed'.

    line_width : float, optional
        Line width, default 0.5.

    name : str, optional
        The name of thr Coverage.
    """

    def __init__(self, vlines, **kwargs):
        properties_dict = {
            "vlines_list": vlines,
            "color": "#1e1e1e",
            "alpha": 0.8,
            "line_style": "dashed",
            "line_width": 1,
        }
        properties_dict.update(kwargs)
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
        High light region alpha value, default 0.5.

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

    name : str, optional
        The name of thr Coverage.
    """

    def __init__(self, file_, **kwargs):
        properties_dict = {
            "file": file_,
            "color": "bed_rgb",
            "alpha": 0.5,
            "border_line": True,
            "border_line_style": "dashed",
            "border_line_width": 1.0,
            "border_line_color": "#000000",
            "border_line_alpha": 0.8,
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)


class HighLights(Coverage, PlotHighLightRegions):

    """
    High light region.

    Parameters
    ----------
    highlight_regions : list of {str, tuple}
        A list of regions for highlights, region can be expressed as a tuple or string.
        region tuple like:
        [('chr1', 100000, 120000), ('chr2', 130000, 150000)]
        region string format: `chr:start-end` like:
        ['chr1:100000-120000', 'chr2:130000-150000'].

    color : str, optional
        High light region color, default HighLights.DEFAULT_COLOR.

    alpha : float, optional
        High light region alpha value, default 0.5

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

    name : str, optional
        The name of thr Coverage.
    """

    DEFAULT_COLOR = "#ff9c9c"

    def __init__(self, highlight_regions, **kwargs):
        properties_dict = {
            "highlight_regions": highlight_regions,
            "color": HighLights.DEFAULT_COLOR,
            "alpha": 0.5,
            "border_line": True,
            "border_line_style": "dashed",
            "border_line_width": 1.0,
            "border_line_color": "#000000",
            "border_line_alpha": 0.8,
        }
        properties_dict.update(kwargs)
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
        Peak alpha value, default 0.6.

    line_width : float, optional
        Peak border line width, default 1.0

    line_style : str, optional
        Border line style, default 'solid'

    fill : bool, optional
        Fill center or not, default False.

    fill_color : str, optional
        Fill color, use 'bed_rgb' for specify color in file,
        default 'bed_rgb'.

    side : {'upper', 'lower', 'both'}
        Plot peak in which side of the matrix.
        NOTE: This parameters is useful only if the Cool track in matrix format.
    """

    def __init__(self, file_, **kwargs):
        properties_dict = {
            "file": file_,
            "color": "bed_rgb",
            "alpha": 0.6,
            "line_width": 1.5,
            "line_style": "solid",
            "fill": False,
            "fill_color": "bed_rgb",
            "side": "both",
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)


class TADCoverage(Coverage, PlotTADCoverage):
    """
    TAD Coverage is used for plot TAD on upper layer of a track.

    Parameters
    ----------
    file_ : str
        Path to the loop file, loop file is a tab splited text file, fields:
        chr1, x1, x2, chr2, y1, y2, [color], ... (other optional fields)

    show_score : bool
        Show bed score or not.
        default False.

    score_font_size : {'auto', int}
        Score text font size.
        default 'auto'

    score_font_color : str
        Score text color.
        default '#000000'

    score_height_ratio : float
        (text tag height) / (TAD height). used for adjust the position of Score text.
        default 0.5

    color : str, optional
        Peak color, use 'bed_rgb' for specify color in file,
        default 'bed_rgb'.

    alpha : float, optional
        Peak alpha value, default 0.3.

    line_color : str, optional
        Border line color, default '#000000'.

    line_width : float, optional
        Border line width, default 1.0.

    line_style : str, optional
        Border line style, default 'solid'.

    fill : bool, optional
        Fill center or not, default True.

    """

    def __init__(self, file_, **kwargs):
        properties_dict = {
            "file": file_,
            "show_score": False,
            "score_font_size": 'auto',
            "score_font_color": "#000000",
            "score_height_ratio": 0.4,
            "color": "bed_rgb",
            "alpha": 0.3,
            "border_color": "#000000",
            "border_width": 1.0,
            "border_style": "solid",
            "fill": True,
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
