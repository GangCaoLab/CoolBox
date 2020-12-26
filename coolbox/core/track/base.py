from copy import copy

from coolbox.utilities import (
    op_err_msg, get_feature_stack, get_coverage_stack,
    split_genome_range, bool2str
)


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

    def __init__(self, properties_dict):
        self.properties = properties_dict
        self.properties = bool2str(self.properties)
        name = self.properties.get('name')
        if name is not None:
            assert isinstance(name, str), "Track name must be a `str`."
        else:
            name = self.__class__.__name__ + ".{}".format(self.__class__._counts)
        self.properties['name'] = name
        super().__init__()
        self.coverages = []

        # load features from global feature stack
        features_stack = get_feature_stack()
        for features in features_stack:
            self.properties[features.key] = features.value

        # load coverages from global coverages stack
        coverage_stack = get_coverage_stack()
        for coverage in coverage_stack:
            self.coverages.append(coverage)

        self.ax = None

    @property
    def name(self):
        return self.properties['name']

    @name.setter
    def name(self, value):
        self.properties['name'] = value

    def __add__(self, other):
        from ..frame import Frame
        from ..coverage.base import Coverage
        from ..coverage.base import CoverageStack
        from ..feature import Feature

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

    def plot_genome_range(self, ax, genome_range):
        """
        Plot the track within a genome range.

        Parameters
        ----------
        ax: matplotlib.axes.Axes
            Axis to use to plot the scale.

        genome_range : {str, GenomeRange}
            Genome range to plot.
        """
        chrom, start, end = split_genome_range(genome_range)
        self.plot(ax, chrom, start, end)

    def plot_y_axis(self, plot_axis, y_ax):
        """
        Plot the scale of the y axis with respect to the plot_axis

        plot something that looks like this:
        ymax ┐
             │
             │
        ymin ┘

        Parameters
        ----------
        plot_axis : matplotlib.axes.Axes
            Main plot axis.

        y_ax : matplotlib.axes.Axes
            Axis to use to plot the scale
        """

        if 'show_data_range' in self.properties and self.properties['show_data_range'] == 'no':
            return

        def value_to_str(value):
            if value % 1 == 0:
                str_value = str(int(value))
            else:
                if value < 0.01:
                    str_value = "{:.4f}".format(value)
                else:
                    str_value = "{:.2f}".format(value)
            return str_value

        ymin, ymax = plot_axis.get_ylim()

        ymax_str = value_to_str(ymax)
        ymin_str = value_to_str(ymin)
        x_pos = [0, 0.5, 0.5, 0]
        y_pos = [0.01, 0.01, 0.99, 0.99]
        y_ax.plot(x_pos, y_pos, color='black', linewidth=1, transform=y_ax.transAxes)
        y_ax.text(-0.2, -0.01, ymin_str, verticalalignment='bottom', horizontalalignment='right',
                  transform=y_ax.transAxes)
        y_ax.text(-0.2, 1, ymax_str, verticalalignment='top', horizontalalignment='right', transform=y_ax.transAxes)
        y_ax.patch.set_visible(False)

    def plot_label(self):
        if hasattr(self, 'label_ax') and self.label_ax is not None:
            self.label_ax.text(0.15, 0.5, self.properties['title'],
                               horizontalalignment='left', size='large',
                               verticalalignment='center')
