from copy import copy

from coolbox.utilities import (
    op_err_msg, get_feature_stack,
    get_coverage_stack, format_properties,
    get_logger
)
from coolbox.utilities.genome import GenomeRange
from coolbox.utilities.doctool import NumpyDocInheritor


log = get_logger(__name__)


class Track(object, metaclass=NumpyDocInheritor):
    """
    Track base class.

    Parameters
    ----------
    properties_dict : dict
        The properties(features) of this track. For example 'height', 'color'...

    name : str, optional
        The name of Track.
        (Default: "{self.__class__.__name__}.{self.__class__._counts}")

    title : str, optional
        Title of ax

    height : int, optional
        Height of ax

    color : str, optional;
        Color of ax

    Attributes
    ----------
    properties : dict
        The properties(features) of this track. For example 'height', 'color'...

    name : str
        The name of Track.

    coverages : list of `coolbox.api.coverage.Coverage`
        Coverages on this Track.
    """

    DEFAULT_PROPERTIES = {
        'height': 2,
        "color": "#808080",
        "name": "",
        "title": "",
    }

    counts = 0
    tracks = {}

    def __new__(cls, *args, **kwargs):
        cls.counts += 1
        # An environment to record all exist tracks
        return super().__new__(cls)

    def __init__(self, properties_dict={}, **kwargs):
        # TODO replace properties_dict with **kwargs ? Then the initialization in all classes are unified
        properties = Track.DEFAULT_PROPERTIES.copy()
        properties.update(properties_dict)
        properties.update(kwargs)
        self.properties = format_properties(properties)

        name = self.properties.get('name')
        if name:
            assert isinstance(name, str), "Track name must be a `str`."
        else:
            name = self.__class__.__name__ + ".{}".format(self.__class__.counts)
        Track.tracks[name] = self
        self.properties['name'] = name
        self.name = name

        # disable call mixin class
        # super().__init__()
        self.coverages = []

        # load features from global feature stack
        features_stack = get_feature_stack()
        for features in features_stack:
            self.properties.update(features.properties)

        # load coverages from global coverages stack
        coverage_stack = get_coverage_stack()
        for coverage in coverage_stack:
            self.coverages.append(coverage)

        self.ax = None

    def __del__(self):
        name = self.properties['name']
        if name in Track.tracks:
            del Track.tracks[name]

    def fetch_data(self, gr: GenomeRange, **kwargs):
        """
        Fetch the raw data for the given GenomeRange.

        Parameters
        ----------
        gr: GenomeRange
        kwargs: dict
            Other optional parameters including `gr2`, `resolution` for some tracks.

        Returns
        -------

        """
        raise NotImplementedError

    def fetch_plot_data(self, gr: GenomeRange, **kwargs):
        """
        Fetch the plot data for the given GenomeRange.
        The returned data would be used in `plot` function, the default implementation directly return `fetch_data`.

        Parameters
        ----------
        gr: GenomeRange
            GenomeRange is an object with properties of [chrom, start, end].

        kwargs: dict
            Other optional parameters including `gr2`, `resolution` for some tracks.

        Returns
        -------

        """
        return self.fetch_data(gr, **kwargs)

    def plot(self, ax, gr: GenomeRange, **kwargs):
        """
        Plot figure based on data returned by `fetch_plot_data`.

        Parameters
        ----------
        gr: GenomeRange
        kwargs: dict
            Other optional parameters including `gr2`, `resolution` for some tracks.

        Returns
        -------

        """
        raise NotImplementedError

    @property
    def name(self):
        return self.properties['name']

    @name.setter
    def name(self, value):
        self.properties['name'] = value

    def __add__(self, other):
        from ..frame.base import FrameBase
        from ..frame.frame import Frame
        from ..coverage.base import Coverage
        from ..coverage.base import CoverageStack
        from ..feature import Feature

        if isinstance(other, Track):
            result = Frame()
            result.add_track(self)
            result.add_track(other)
            return result
        elif isinstance(other, FrameBase):
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
        from ..coverage.base import CoverageStack

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

    def has_prop(self, p):
        return (p in self.properties) and (self.properties[p] is not None)

    def plot_label(self):
        if hasattr(self, 'label_ax') and self.label_ax is not None:
            self.label_ax.text(0.15, 0.5, self.properties['title'],
                               horizontalalignment='left', size='large',
                               verticalalignment='center')

    def plot_coverages(self, ax, gr: GenomeRange, gr2: GenomeRange):
        """
        Plot all coverages on given axes.

        Parameters
        ----------
        gr : GenomeRange
            First genome range.

        gr2 : GenomeRange
            Second genome range.
        """
        if not hasattr(self, 'coverages'):
            return
        for cov_idx, cov in enumerate(self.coverages):
            cov.track = self
            if hasattr(cov, 'track_instance'):
                cov.track_instance.track = self
            try:
                cov.plot(ax, copy(gr), gr2=copy(gr2))
            except Exception as e:
                log.error("Error occured when plot track's coverage:\n"
                          "\ttrack name: {}\n\ttrack type:{}\n\tcoverage name: {}\n\tcov type: {}\n"
                          "\tError: {} {}".format(
                            self.name, type(self), cov.name, type(cov),
                            type(e), str(e)))
                log.exception(e)
