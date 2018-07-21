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

    def __init__(self, key, value):
        self.key = key
        self.value = value

    def __add__(self, other):
        from .track import Track
        from .frame import Frame
        from .coverage import Coverage
        from .coverage import CoverageStack

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
        from .frame import Frame

        if isinstance(other, Frame):
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


class HistStyle(Feature):
    """
    Style of BigWig or BedGraph.
    """
    def __init__(self, type='fill', size=0.5):
        if type != 'fill':
            value = type + ":" + str(size)
        else:
            value = 'fill'
        super().__init__('type', value)


class ShowDataRange(Feature):
    """
    Show data range or not.
    """
    def __init__(self, show=True):
        if show:
            super().__init__('show_data_range', 'yes')
        else:
            super().__init__('show_data_range', 'no')


class ShowColorBar(Feature):
    """
    Show color bar or not.
    """
    def __init__(self, show=True):
        if show:
            super().__init__('color_bar', 'yes')
        else:
            super().__init__('color_bar', 'no')


class DepthRatio(Feature):
    """
    Control Cool track's depth ratio.
    """
    def __init__(self, depth_ratio):
        super().__init__('depth_ratio', depth_ratio)


class CoolStyle(Feature):
    """
    Control Cool track's style.
    """
    def __init__(self, style='triangular'):
        super().__init__('style', style)


class FrameFeature(Feature):
    """
    FrameFeature base class.
    """

    def __add__(self, other):
        from .frame import Frame

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


if __name__ == "__main__":
    import doctest
    doctest.testmod()
