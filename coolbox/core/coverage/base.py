from copy import copy

from coolbox.utilities import (
    op_err_msg, get_coverage_stack, get_feature_stack,
    format_properties, GenomeRange
)


class Coverage(object):
    """
    Coverage base class.

    `Coverage` is the plots at the upper layer of Track.

    >>> from coolbox.core.track.base import Track
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
        self.properties = format_properties(properties_dict)
        name = self.properties.get("name")
        if name is not None:
            assert isinstance(name, str), "Coverage name must be a `str`."
        else:
            name = self.__class__.__name__ + ".{}".format(self.__class__._counts)
        self.properties['name'] = name

        stack = get_feature_stack()
        for feature in stack:
            self.properties.update(feature.properties)

        self.track = None

    @property
    def name(self):
        return self.properties['name']

    @name.setter
    def name(self, value):
        self.properties['name'] = value

    def __add__(self, other):
        from ..track.base import Track
        from ..frame.base import FrameBase
        from ..feature import Feature

        if isinstance(other, Track):
            result = copy(other)
            result.append_coverage(self)
            return result
        elif isinstance(other, FrameBase):
            result = copy(other)
            if len(result.tracks) > 1:
                first = list(result.tracks.values())[0]
                first.append_coverage(self, pos='bottom')
            return result
        elif isinstance(other, Feature):
            result = copy(self)
            result.properties.update(other.properties)
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
        from ..frame.base import FrameBase
        if isinstance(other, FrameBase):
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
            msg = "{} coverage's track must be a instance of {}".format(
                self.track.__class__.__name__,
                [type_.__name__ for type_ in allow]
            )
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
        from ..track.base import Track
        from ..frame.base import FrameBase
        from ..feature import Feature

        if isinstance(other, Coverage):
            result = copy(self)
            result.to_top(other)
            return result
        elif isinstance(other, Track):
            result = copy(other)
            result.pile_coverages(self.coverages, pos='bottom')
            return result
        elif isinstance(other, FrameBase):
            result = copy(other)
            if len(result.tracks) != 0:
                first = list(result.tracks.values())[0]
                first.pile_coverages(self.coverages, pos='bottom')
            return result
        elif isinstance(other, Feature):
            result = copy(self)
            if len(result.coverages) != 0:
                last = result.coverages[-1]
                last.properties.update(other.properties)
            return result
        else:
            raise TypeError(op_err_msg(self, other))


def has_data_range(track_class):
    from coolbox.core.track import HistBase, Hist, ArcsBase, Arcs
    if isinstance(track_class, type):
        return HistBase in track_class.__mro__ or ArcsBase in track_class.__mro__
    else:
        return track_class in (Hist, Arcs)


def track_to_coverage(track_class):
    def init(self, *args, **kwargs):
        if has_data_range(track_class):
            kwargs.update({
                "show_data_range": False,
            })
        self.track_instance = track_class(*args, **kwargs)
        self.properties = self.track_instance.properties

    def plot(self, ax, gr: GenomeRange, **kwargs):
        if hasattr(self, 'track'):
            if has_data_range(track_class):
                # update height when plot
                self.track_instance.properties['height'] = self.track.properties['height']
        self.track_instance.plot(ax, gr, **kwargs)

    # TODO , other track methods
    cov_class = type(
        track_class.__name__ + "Coverage",
        (Coverage,),
        {
            "__init__": init,
            "__doc__": track_class.__doc__,
            "plot": plot
        }
    )

    return cov_class
