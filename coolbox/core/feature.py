from copy import copy

from coolbox.utilities import op_err_msg, get_feature_stack

__all__ = [
    "Feature", "Color", "ColorMap",
    "TrackHeight", "Inverted", "Title",
    "MaxValue", "MinValue", "HistStyle", "ShowDataRange",
    "ShowColorBar", "DepthRatio", "CoolStyle",
    "FrameFeature", "FrameTitle"
]


class Feature(object):
    """
    Feature base class.
    """

    def __init__(self, **kwargs):
        self.properties = kwargs

    def __add__(self, other):
        from .track.base import Track
        from .frame.base import FrameBase
        from .coverage.base import Coverage
        from .coverage.base import CoverageStack

        if isinstance(other, Track):
            result = copy(other)
            result.properties.update(self.properties)
            return result
        elif isinstance(other, FrameBase):
            result = copy(other)
            if len(result.tracks) != 0:
                first = list(result.tracks.values())[0]
                first.properties.update(self.properties)
            return result
        elif isinstance(other, Coverage):
            result = copy(other)
            result.properties.update(self.properties)
            return result
        elif isinstance(other, CoverageStack):
            result = copy(other)
            if len(result.coverages) != 0:
                first = result.coverages[0]
                first.properties.update(self.properties)
            return result
        else:
            raise TypeError(op_err_msg(self, other))

    def __mul__(self, other):
        from .frame.base import FrameBase

        if isinstance(other, FrameBase):
            result = copy(other)
            result.add_feature_to_tracks(self)
            return result
        else:
            raise TypeError(op_err_msg(self, other, op='*'))

    def __enter__(self):
        stack = get_feature_stack()
        stack.append(self)
        return self

    def __exit__(self, type, value, traceback):
        stack = get_feature_stack()
        stack.pop()


class Color(Feature):
    """
    Track color.
    """

    def __init__(self, value):
        super().__init__(color=value)


class ColorMap(Feature):
    """
    Track color map.
    """

    def __init__(self, value):
        # TODO which one?
        super().__init__(color=value, cmap=value)


class TrackHeight(Feature):
    """
    Track height.
    """

    def __init__(self, value):
        super().__init__(height=value)


class Inverted(Feature):
    """
    Invert the orientation of track.
    """

    def __init__(self):
        super().__init__(orientation='inverted')


class Title(Feature):
    """
    Track title.
    """

    def __init__(self, value):
        super().__init__(title=value)


class MaxValue(Feature):
    """
    Max value of track.
    """

    def __init__(self, value):
        super().__init__(max_value=value)


class MinValue(Feature):
    """
    Min value of track.
    """

    def __init__(self, value):
        super().__init__(min_value=value)


class HistStyle(Feature):
    """
    Style of BigWig or BedGraph.
    """

    def __init__(self, style='fill', fmt="-", size=10, line_width=2.0):
        from coolbox.core.track import HistBase
        assert style in HistBase.STYLES, f"The style should be one of {HistBase.STYLES}"
        super().__init__(style=style, fmt=fmt, size=size, line_width=line_width)


class ShowDataRange(Feature):
    """
    Show data range or not.
    """

    def __init__(self, data_range_style="y-axis"):
        assert data_range_style in ("no", False, "text", "y-axis"), "Should be one of ['no', 'text', 'y-axis']"
        super().__init__(data_range_style=data_range_style)


class ShowColorBar(Feature):
    """
    Show color bar or not.
    """

    def __init__(self, show=True):
        super().__init__(color_bar='yes' if bool(show) else 'no')


class DepthRatio(Feature):
    """
    Control Cool track's depth ratio.
    """

    def __init__(self, depth_ratio):
        super().__init__(depth_ratio=depth_ratio)


class CoolStyle(Feature):
    """
    Control Cool track's style.
    """

    def __init__(self, style='triangular'):
        super().__init__(style=style)


class FrameFeature(Feature):
    """
    FrameFeature base class.
    """

    def __add__(self, other):
        from .frame.base import FrameBase

        if isinstance(other, FrameBase):
            properties = self.properties.copy()
            properties.update(other.properties)
            result = copy(other)
            result.properties = properties
            return result
        else:
            raise TypeError(op_err_msg(self, other))


class FrameTitle(FrameFeature):
    """
    Frame title.
    """

    def __init__(self, value):
        super().__init__(title=value)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
