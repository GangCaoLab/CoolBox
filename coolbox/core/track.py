from copy import copy

from coolbox.plots.track import *
from coolbox.fetchdata import *
from coolbox.utilities import op_err_msg, get_feature_stack, get_coverage_stack


__all__ = [
    "Spacer", "HLine", "XAxis", "Bed", "TADs",
    "BigWig", "ABCompartment", "BedGraph",
    "Arcs", "Cool", "DotHiC", "HicCompare", "Virtual4C",
    "Ideogram",
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

    def __init__(self, properties_dict):
        self.properties = properties_dict
        self.__bool2str()
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
        Track's name.
    """

    DEFAULT_HEIGHT = 1

    def __init__(self, *args, **kwargs):
        properties_dict = {
            'height': Spacer.DEFAULT_HEIGHT,
        }
        if args:
            height = args[0]
            properties_dict['height'] = height
        properties_dict.update(kwargs)

        super().__init__(properties_dict)


class HLine(Track, PlotHLine):
    """
    Horizontal line track.
    Used for add a horizontal line between two tracks.

    Parameters
    ----------
    line_width : float, optional
        (Default: HLine.DEFAULT_LINE_WIDTH)

    line_style : str, optional
        (Default: HLine.DEFAULT_LINE_STYLE)

    color : str, optional
        (Default: HLine.DEFAULT_COLOR)

    alpha : float, optional
        (Default: HLine.DEFAULT_ALPHA)

    height : float, optional
        The height of Spacer track. (Default: HLine.DEFAULT_HEIGHT)

    name : str, optional
        Track's name.
    """

    DEFAULT_LINE_WIDTH = 1.0
    DEFAULT_HEIGHT = max(DEFAULT_LINE_WIDTH / 50, 0.05)  # this just a empiric value
    DEFAULT_LINE_STYLE = '--'
    DEFAULT_COLOR = '#000000'
    DEFAULT_ALPHA = 0.75

    def __init__(self, **kwargs):
        properties_dict = {
            'height': HLine.DEFAULT_HEIGHT,
            'line_width': HLine.DEFAULT_LINE_WIDTH,
            'line_style': HLine.DEFAULT_LINE_STYLE,
            'color': HLine.DEFAULT_COLOR,
            'alpha': HLine.DEFAULT_ALPHA,
        }
        properties_dict.update(**kwargs)
        super().__init__(properties_dict)


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
        Track's name.
    """

    DEFAULT_FONTSIZE = 15
    DEFAULT_HEIGHT = 2

    def __init__(self, **kwargs):
        properties_dict = {
            'height': XAxis.DEFAULT_HEIGHT,
            'fontsize': XAxis.DEFAULT_FONTSIZE,
            'where': 'bottom',
        }
        properties_dict.update(kwargs)

        super().__init__(properties_dict)


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
        (Default: 'auto')

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
        Track's name.

    """

    DEFAULT_FONTSIZE = 12

    def __init__(self, file_, **kwargs):
        properties_dict = {
            'file': file_,
            'height': Bed.DEFAULT_HEIGHT,
            'color': 'bed_rgb',
            'border_color': 'black',
            'fontsize': Bed.DEFAULT_FONTSIZE,
            'title': '',
            'labels': 'auto',
            'style': 'flybase',
            'display': 'stacked',
            'global_max_row': False,
        }
        properties_dict.update(kwargs)

        labels = properties_dict.get('labels')
        if labels == 'auto':
            properties_dict['labels'] = 'auto'
        elif labels is True:
            properties_dict['labels'] = 'on'
        else:
            properties_dict['labels'] = 'off'

        super().__init__(properties_dict)


class BigWig(Track, PlotBigWig, FetchBigWig):
    """
    BigWig track.

    Parameters
    ----------
    file_ : str
        Path to bigwig file.

    height : float, optional
        Height of track, default BigWig.DEFAULT_HEIGHT.

    color : str, optional
        Track color, default BigWig.DEFAULT_COLOR.

    number_of_bins : int, optional
        Number_of_bins in current range, default 700.

    style : str, optional
        Track graph type, format {'fill', 'line:`size`', 'points:`size`'},
        example: 'line:2', 'points:0.5'. default: 'fill'

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    show_data_range : bool, optional
        Show_data_range or not, default True.

    data_range_style : {'text', 'y-axis'}
        The style of the data range. default: 'y-axis'

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

    def __init__(self, file_, **kwargs):
        properties_dict = {
            'file': file_,
            'height': BigWig.DEFAULT_HEIGHT,
            'color': BigWig.DEFAULT_COLOR,
            'number_of_bins': 700,
            'style': 'fill',
            'show_data_range': True,
            'data_range_style': 'y-axis',
            'title': '',
            'max_value': 'auto',
            'min_value': 'auto',
        }
        properties_dict.update(kwargs)
        properties_dict['type'] = properties_dict['style']  # change key word

        super().__init__(properties_dict)


class ABCompartment(BigWig):
    """
    A/B Comapartment BigWig track.

    Parameters
    ----------
    file_ : str
        Path to bigwig file.

    height : float, optional
        Height of track, default BigWig.DEFAULT_HEIGHT.

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

    def __init__(self, file_, **kwargs):
        properties_dict = {
            'positive_color': ABCompartment.DEFAULT_POSITIVE_COLOR,
            'negative_color': ABCompartment.DEFAULT_NEGATIVE_COLOR,
        }
        properties_dict.update(kwargs)
        super().__init__(file_, **properties_dict)


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

    style : str, optional
        Track graph type, format {'fill', 'line:`size`', 'points:`size`'},
        example: 'line:2', 'points:0.5'. default: 'fill'

    extra : optional

    show_data_range : bool, optional
        Show_data_range or not, default True.

    data_range_style : {'text', 'y-axis'}
        The style of the data range. default: 'y-axis'

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

    def __init__(self, file_, **kwargs):

        properties_dict = {
            'file': file_,
            'height': BedGraph.DEFAULT_HEIGHT,
            'color': BedGraph.DEFAULT_COLOR,
            'style': 'fill',
            'show_data_range': True,
            'data_range_style': 'y-axis',
            'title': '',
            'max_value': 'auto',
            'min_value': 'auto',
        }
        properties_dict.update(kwargs)
        properties_dict['type'] = properties_dict['style']  # change key word

        super().__init__(properties_dict)


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

    def __init__(self, file_, **kwargs):

        properties_dict = {
            'file': file_,
            'height': Arcs.DEFAULT_HEIGHT,
            'color': Arcs.DEFAULT_COLOR,
            'alpha': Arcs.DEFAULT_ALPHA,
            'title': '',
        }
        properties_dict.update(kwargs)

        super().__init__(properties_dict)


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

    def __init__(self, file_, **kwargs):

        properties_dict = {
            "file": file_,
            "height": TADs.DEFAULT_HEIGHT,
            "color": TADs.DEFAULT_COLOR,
            "border_color": 'black',
            "title": '',
        }
        properties_dict.update(kwargs)

        super().__init__(properties_dict)


class Cool(Track, PlotCool, FetchCool):
    """
    Cool Hi-C matrix (or triangular matrix) track.

    Parameters
    ----------
    file_ : str
        Path to bed file.

    cmap : str, optional
        Color map of hic matrix, default Cool.DEFAULT_COLOR.

    style : {'triangular', 'window', 'matrix'}, optional
        Matrix style, default 'triangular'.

    balance : bool, optional
        Show balanced matrix or not, default True

    depth_ratio : float, optional
        Depth ratio of triangular matrix, use 'full' for full depth. default 'full'.

    color_bar : bool, optional
        Show color_bar or not, default True.

    transform : {str, bool}, optional
        Transform for matrix, like 'log2', 'log10', default False.

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

    def __init__(self, file_, **kwargs):

        properties_dict = {
            "file": file_,
            "cmap": Cool.DEFAULT_COLOR,
            "style": 'triangular',
            "balance": True,
            "depth_ratio": "full",
            "color_bar": True,
            "transform": False,
            "norm": 'log',
            "max_value": "auto",
            "min_value": "auto",
            "title": '',
        }
        properties_dict.update(kwargs)
        properties_dict['color'] = properties_dict['cmap']

        super().__init__(properties_dict)


class DotHiC(Track, PlotDotHiC, FetchDotHiC):

    """
    .hic Hi-C matrix (or triangular matrix) track.

    Parameters
    ----------
    file_ : str
        Path to bed file.

    cmap : str, optional
        Color map of hic matrix, default Cool.DEFAULT_COLOR.

    style : {'triangular', 'window', 'matrix'}, optional
        Matrix style,
        default 'triangular'.

    balance : {bool, 'KR', 'VC', 'VC_SQRT'}, optional
        Matrix balance method,
        default True('KR' balance)

    depth_ratio : float, optional
        Depth ratio of triangular matrix, use 'full' for full depth. default 'full'.

    color_bar : bool, optional
        Show color_bar or not, default True.

    transform : {str, bool}, optional
        Transform for matrix, like 'log2', 'log10', default False.

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
    DEFAULT_COLOR = 'Reds'

    def __init__(self, file_, **kwargs):

        properties_dict = {
            "file": file_,
            "cmap": Cool.DEFAULT_COLOR,
            "style": 'triangular',
            "balance": True,
            "depth_ratio": "full",
            "color_bar": True,
            "transform": False,
            "norm": 'log',
            "max_value": "auto",
            "min_value": "auto",
            "title": '',
        }
        properties_dict.update(kwargs)
        properties_dict['color'] = properties_dict['cmap']

        super().__init__(properties_dict)


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

    title : str, optional
        Label text, default ''.

    name : str, optional
        Track's name

    """

    DEFAULT_COLOR = 'bwr'

    def __init__(self, hic1, hic2, **kwargs):
        properties_dict = {
            "hic1": hic1,
            "hic2": hic2,
            "cmap": HicCompare.DEFAULT_COLOR,
            "color_bar": True,
            "title": '',
        }
        properties_dict.update(kwargs)
        properties_dict['color'] = properties_dict['cmap'] # change key word

        super().__init__(properties_dict)


class Virtual4C(Track, PlotVirtual4C, FetchVirtual4C):
    """
    Track for view virtual 4C related to a certain genome position,
    and a HiC Track (include `Cool` and `DotHiC`).

    Parameters
    ----------
    hic_track : {`Cool`, `DotHiC`}
        related hic track.

    genome_position : str
        related genome position, like: 'chr1:2000000-2000000'

    bin_width : int, optional
        How many bin used for calculate the mean value.
        default 3

    color : str, optional
        Track color.

    height : int, optional
        Track height

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    max_value : {float, 'auto'}, optional
        Max value of track, use 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of track, use 'auto' for specify min value automatically, default 'auto'.

    show_data_range : bool, optional
        Show_data_range or not, default True.

    data_range_style : {'text', 'y-axis'}
        The style of the data range. default: 'y-axis'

    style : str, optional
        Track graph type, format {'fill', 'line:`size`', 'points:`size`'},
        example: 'line:2', 'points:0.5'. default: 'line:2'

    title : str, optional
        Label text, default ''.

    name : str, optional
        Track's name.

    """

    DEFAULT_COLOR = '#2855d8'

    def __init__(self, hic_track, genome_position, **kwargs):
        properties_dict = {
            'hic': hic_track,
            'color': Virtual4C.DEFAULT_COLOR,
            'height': Virtual4C.DEFAULT_HEIGHT,
            'genome_position': genome_position,
            'bin_width': 3,
            'max_value': 'auto',
            'min_value': 'auto',
            'show_data_range': True,
            'data_range_style': 'y-axis',
            'style': 'line:1',
            'title': '',
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)


class Ideogram(Track, PlotIdeogram):
    """
    The chromosome ideograme track.

    Parameters
    ----------
    file_ : str
        Path to chromosome ideogram txt file,
        ideogram file is download from the UCSC Table Browser CytoBandIdeo table (in "all table" group).
        see: http://genome.ucsc.edu/cgi-bin/hgTables?hgta_group=allTables&hgta_table=cytoBandIdeo
    color_scheme : dict, optional
        Color scheme of ideogram, default: Ideogram.DEFAULT_COLOR_SCHEME
    show_name : bool, optional
        Show band name or not. default True.
    font_size : int, optional
        Band name font size.
    border_color : str, optional
        Border color. default: '#000000'
    border_width : float, optional
        Border width. default: 1.2
    height : int, optional
        Track height.
    title : str, optional
        Label text, default ''.
    name : str, optional
        Track's name.
    """

    DEFAULT_HEIGHT = 1.2
    DEFAULT_COLOR_SCHEME = {
        'gneg':    '#ffffff',
        'gpos25':  '#999999',
        'gpos50':  '#666666',
        'gpos75':  '#333333',
        'gpos100': '#000000',
        'acen':    '#cc6666',
        'gvar':    '#cccccc',
        'stalk':   '#e5e5e5',
    }
    DEFAULT_FONT_SIZE = 12
    DEFAULT_BORDER_WIDTH = 1.2

    def __init__(self, file_, **kwargs):
        properties_dict = {
            'file': file_,
            'color_scheme': Ideogram.DEFAULT_COLOR_SCHEME,
            'show_name': True,
            'font_size': Ideogram.DEFAULT_FONT_SIZE,
            'border_color': '#000000',
            'border_width': Ideogram.DEFAULT_BORDER_WIDTH,
            'height': Ideogram.DEFAULT_HEIGHT,
            'title': '',
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
