from coolbox.utilities import (
    file_to_intervaltree, GenomeRange,
    to_gr,
)
from .base import Coverage


class VlinesBase(object):
    def fetch_data(self, gr: GenomeRange):
        vlines_list = []

        if gr.chrom not in list(self.vlines_intval_tree):
            gr.change_chrom_names()

        for region in sorted(self.vlines_intval_tree[gr.chrom][gr.start - 1:gr.end + 1]):
            vlines_list.append(region.begin)
            if region.end != region.begin:
                vlines_list.append(region.end)

        return vlines_list

    def plot(self, ax, gr: GenomeRange, **kwargs):
        gr = GenomeRange(gr)
        vlines_list = self.fetch_data(gr)

        ymin, ymax = ax.get_ylim()

        ax.vlines(vlines_list, ymin, ymax,
                  linestyle=self.properties['line_style'],
                  linewidth=self.properties['line_width'],
                  color=self.properties['color'],
                  alpha=self.properties['alpha'])


class VlinesFromFile(Coverage, VlinesBase):
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

    DEFAULT_LINE_WIDTH = 0.5
    DEFAULT_LINE_STYLE = 'dashed'
    DEFAULT_COLOR = '#1e1e1e'
    DEFAULT_ALPHA = 0.8

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
        self.vlines_intval_tree, _, _ = file_to_intervaltree(self.properties['file'])


class Vlines(Coverage, VlinesBase):
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
        if not isinstance(vlines, list):
            vlines = [vlines]
        properties_dict = {
            "vlines_list": vlines,
            "color": "#1e1e1e",
            "alpha": 0.8,
            "line_style": "dashed",
            "line_width": 1,
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)
        self.vlines_intval_tree = self.__intervaltree_from_list(self.properties['vlines_list'])

    def __intervaltree_from_list(self, vlines_list):
        from intervaltree import IntervalTree
        itree = {}
        for v in vlines_list:
            if isinstance(v, str):
                grange = GenomeRange(v)
            elif isinstance(v, tuple):
                grange = GenomeRange(v[0], v[1], v[1])
            elif isinstance(v, GenomeRange):
                grange = v
            else:
                raise ValueError("position must be a tuple or string.")
            chr_ = grange.chrom
            itree.setdefault(chr_, IntervalTree())
            itree[chr_][grange.start:grange.end + 1] = grange
        return itree
