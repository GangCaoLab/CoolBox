from coolbox.utilities import (
    file_to_intervaltree, change_chrom_names, GenomeRange
)

from coolbox.plots.coverage.base import CoveragePlot


class PlotVlines(CoveragePlot):

    DEFAULT_LINE_WIDTH = 0.5
    DEFAULT_LINE_STYLE = 'dashed'
    DEFAULT_COLOR = '#1e1e1e'
    DEFAULT_ALPHA = 0.8

    def __init__(self, *args, **kwargs):

        CoveragePlot.__init__(self, *args, **kwargs)
        if 'line_width' not in self.properties:
            self.properties['line_width'] = PlotVlines.DEFAULT_LINE_WIDTH
        if 'line_style' not in self.properties:
            self.properties['line_style'] = PlotVlines.DEFAULT_LINE_STYLE
        if 'color' not in self.properties:
            self.properties['color'] = PlotVlines.DEFAULT_COLOR
        if 'alpha' not in self.properties:
            self.properties['alpha'] = PlotVlines.DEFAULT_ALPHA

        if 'file' in self.properties:
            # plot vlines from file
            self.vlines_intval_tree, _, _ = file_to_intervaltree(self.properties['file'])
        else:
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
            itree[chr_][grange.start:grange.end+1] = grange
        return itree

    def __get_vlines(self, chrom, start, end):
        vlines_list = []

        if chrom not in list(self.vlines_intval_tree):
            chrom = change_chrom_names(chrom)

        for region in sorted(self.vlines_intval_tree[chrom][start-1:end+1]):
            vlines_list.append(region.begin)

        return vlines_list

    def plot(self, ax, chrom_region, start_region, end_region):
        vlines_list = self.__get_vlines(chrom_region, start_region, end_region)

        ymin, ymax = ax.get_ylim()

        ax.vlines(vlines_list, ymin, ymax,
                  linestyle=self.properties['line_style'],
                  linewidth=self.properties['line_width'],
                  color=self.properties['color'],
                  alpha=self.properties['alpha'])
