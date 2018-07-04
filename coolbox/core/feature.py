from collections import deque
from copy import copy

from coolbox.utilities import op_err_msg


FEATURES_STACK_NAME = "__COOLBOX_FEATURE_STACK__"
global_scope = globals()
global_scope[FEATURES_STACK_NAME] = deque()


__all__ = [
    "Feature", "Color", "ColorMap",
    "TrackHeight", "Inverted", "Title",
    "MaxValue", "MinValue", "ShowDataRange",
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
