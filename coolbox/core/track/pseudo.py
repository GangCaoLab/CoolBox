from coolbox.utilities import GenomeRange
from .base import Track


# TODO change properties to new mode
class Spacer(Track):
    """
    The spacer track,
    not have any real content, just used to split two tracks.

    Parameters
    ----------
    height : float, optional
        The height of Spacer track. (Default: Spacer.DEFAULT_HEIGHT)

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

    def fetch_data(self, gr: GenomeRange, **kwargs):
        pass

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax
        ax.set_xlim(gr.start, gr.end)


class HLine(Track):
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

    def fetch_data(self, gr: GenomeRange, **kwargs):
        pass

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax

        ax.set_xlim(gr.start, gr.end)
        ax.hlines(0, gr.start, gr.end,
                  linestyles=self.properties['line_style'],
                  linewidth=self.properties['line_width'],
                  colors=self.properties['color'],
                  alpha=self.properties['alpha'])


class XAxis(Track):
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

    def fetch_data(self, gr: GenomeRange, **kwargs):
        pass

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax

        ax.set_xlim(gr.start, gr.end)
        ticks = ax.get_xticks()
        if ticks[-1] - ticks[1] <= 1e5:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"

        elif 1e5 < ticks[-1] - ticks[1] < 4e6:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"
        else:
            labels = ["{:,.1f} ".format((x / 1e6))
                      for x in ticks]
            labels[-2] += " Mbp"

        ax.axis["x"] = ax.new_floating_axis(0, 0.5)

        ax.axis["x"].axis.set_ticklabels(labels)
        ax.axis['x'].axis.set_tick_params(which='minor', bottom='on')

        ax.axis["x"].major_ticklabels.set(size=int(self.properties['fontsize']))

        if 'where' in self.properties and self.properties['where'] == 'top':
            ax.axis["x"].set_axis_direction("top")


class ChromName(Track):
    """
    Track for show chromosome name.

    Parameters
    ----------
    fontsize : float
        Font name to show.

    offset : float
        Offset ratio to the start position.
    """

    def __init__(self, fontsize=50, offset=0.45):
        super().__init__({
            "fontsize": fontsize,
            "offset": offset,
        })

    def fetch_data(self, gr: GenomeRange, **kwargs):
        return gr.chrom  # return chromosome name

    def plot(self, ax, gr: GenomeRange, **kwargs):
        x = gr.start + self.properties['offset'] * (gr.end - gr.start)
        ax.text(x, 0, gr.chrom, fontsize=self.properties['fontsize'])
        ax.set_xlim(gr.start, gr.end)
