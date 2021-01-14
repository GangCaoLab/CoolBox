from typing import Union

from .base import Coverage
from coolbox.utilities.genome import GenomeRange


class HLines(Coverage):
    """
    Horizontal line coverage, for show threshold.

    Parameters
    ----------
    values : float
        Y position. A list of y value or a single value.

    percent_mode : bool
        If in percent mode, y position will equal to ymin + (ymax - ymin) * val.

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

    def __init__(self, values, **kwargs):
        if not isinstance(values, list):
            values = [values]
        properties_dict = {
            "values": values,
            "percent_mode": False,
            "color": "#1e1e1e",
            "alpha": 0.8,
            "line_style": "dashed",
            "line_width": 1,
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)

    def plot(self, ax, gr: GenomeRange, **kwargs):
        gr = GenomeRange(gr)
        ymin, ymax = ax.get_ylim()
        if self.properties['percent_mode'] != 'no':
            hlines_list = [ymin + val * ymax for val in self.properties['values']]
        else:
            hlines_list = self.properties['values']

        ax.hlines(hlines_list, gr.start, gr.end,
                  linestyle=self.properties['line_style'],
                  linewidth=self.properties['line_width'],
                  color=self.properties['color'],
                  alpha=self.properties['alpha'])
