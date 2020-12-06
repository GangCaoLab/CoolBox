import abc

import numpy as np
import pandas as pd
from matplotlib.patches import Arc

from coolbox.utilities import (
    to_gr, GenomeRange
)
from coolbox.utilities.bed import (
    pairix_query, process_bedpe, process_pairs
)
from .base import Track


class _Arcs(Track, abc.ABC):
    """
    Arcs(link) track.

    Parameters
    ----------
    file_ : str
        Path to file.

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

    point_at : str, optional
        Link anchor point method: 'start', 'end', or 'mid', default 'mid'

    score_to_width : str, optional
        Mapping function of score to width, default: '0.5 + sqrt(score)'

    diameter_to_height : str, optional
        Mapping function of arc diameter(interval end - start) to height.
        You can specify to 'max_diameter' let all arcs has same height.
        default 'max_height * diameter / max_diameter'.

    title : str, optional
        Label text. default ''

    name : str, optional
        Track's name.
    """

    DEFAULT_HEIGHT = 2.0
    DEFAULT_COLOR = '#3297dc'
    DEFAULT_ALPHA = 0.8

    def __init__(self, file_, **kwargs):

        properties_dict = {
            'file': file_,
            'height': self.DEFAULT_HEIGHT,
            'color': self.DEFAULT_COLOR,
            'alpha': self.DEFAULT_ALPHA,
            'title': '',
            'point_at': 'mid',
            'score_to_width': '0.5 + sqrt(score)',
            'diameter_to_height': 'max_height * diameter / max_diameter',
        }
        properties_dict.update(kwargs)

        super().__init__(properties_dict)

    def fetch_data(self, genome_range):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        Return
        ------
        intervals : pandas.core.frame.DataFrame
            Arcs interval table.
        """
        return self.fetch_intervals(genome_range)

    def fetch_intervals(self, genome_range, open_query=True):
        gr = to_gr(genome_range)
        rows = self.load(gr, open_query)
        fields = self.fields
        if len(rows) == 0:
            gr.change_chrom_names()
            rows = self.load(gr, open_query)
        if len(rows) == 0:
            return pd.DataFrame([], columns=fields)

        _row = rows[0]
        _diff = len(_row) - len(fields)
        if _diff > 0:
            fields += [f"extra_{i}" for i in range(_diff)]
        df = pd.DataFrame(rows, columns=fields)
        df = self.convert_type(df, _row)
        return df

    def load(self, gr, open_):
        rows = []
        second = gr.chrom if open_ else None
        g = pairix_query(self.bgz_file, str(gr), second=second, split=True)
        for row_item in g:
            rows.append(row_item)
        return rows

    @abc.abstractmethod
    def convert_type(self, df, _row):
        pass

    def plot_arcs(self, ax, genome_range, intervals):
        """
        Parameters
        ----------
        intervals : List[Tuple(int, int, float)]
            List of intervals (start, end, score).
        """
        gr = to_gr(genome_range)
        max_height = self.properties.get('height', 1.0)
        alpha = self.properties.get('alpha', 1.0)
        self.__adjust_yaxis(ax, max_height)
        ax.set_xlim(gr.start, gr.end)
        color = self.properties['color']

        if len(intervals) == 0:
            return

        max_itv = max(intervals, key=lambda t: t[1] - t[0])
        max_diameter = max_itv[1] - max_itv[0]

        for itv in intervals:
            start, end, score = itv
            line_width = self.__get_linewidth(score)
            diameter = (end - start)
            height = 2 * self.__get_height(max_height, max_diameter, diameter)
            center = (start + end) / 2
            ax.plot([center], [diameter])
            arc = Arc(
                (center, 0), diameter,
                height, 0, 0, 180,
                color=color,
                alpha=alpha,
                lw=line_width,
            )
            ax.add_patch(arc)

    def __get_height(self, max_height, max_diameter, diameter):
        max_height = 0.97 * max_height
        key_ = 'diameter_to_height'
        if key_ in self.properties:
            try:
                return eval(self.properties[key_])
            except Exception:
                pass
        return max_height * diameter / max_diameter

    def __get_linewidth(self, score):
        if 'line_width' in self.properties:
            return float(self.properties['line_width'])
        elif 'score_to_width' in self.properties:
            try:
                return eval(self.properties['score_to_width'])
            except Exception:
                pass
        line_width = 0.5 * np.sqrt(score)
        return line_width

    def __adjust_yaxis(self, ax, height):
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(height, 0.001)
        else:
            ax.set_ylim(-0.001, height)


class BEDPE(_Arcs):
    __doc__ = _Arcs.__doc__
    DEFAULT_COLOR = '#3297dc'
    DEFAULT_ALPHA = 0.8
    fields = ["chrom1", "start1", "end1", "chrom2", "start2", "end2",
              "name", "score", "strand1", "strand2"]

    def __init__(self, file_, **kwargs):
        super().__init__(file_, **kwargs)
        self.bgz_file = process_bedpe(file_)

    def convert_type(self, df, _row):
        df['start1'] = df['start1'].astype(int)
        df['end1'] = df['end1'].astype(int)
        df['start2'] = df['start2'].astype(int)
        df['end2'] = df['start2'].astype(int)
        try:
            float(_row[7])
            df['score'] = df['score'].astype(float)
        except ValueError:
            pass
        return df

    def plot(self, ax, region_chrom, region_start, region_end):
        """
        Makes and arc connecting two points on a linear scale representing
        interactions between bins.
        """
        self.ax = ax
        gr = GenomeRange(region_chrom, region_start, region_end)
        itv_df = self.fetch_intervals(gr)
        point_at = self.properties.get('point_at', 'mid')
        intervals = []
        for _, row in itv_df.iterrows():
            if point_at == 'start':
                s, e = row['start1'], row['start2']
            elif point_at == 'end':
                s, e = row['end1'], row['end2']
            else:
                s = (row['start1'] + row['end1']) // 2
                e = (row['start2'] + row['end2']) // 2
            s, e = sorted([s, e])
            intervals.append((s, e, row['score']))
        self.plot_arcs(ax, gr, intervals)
        self.plot_label()


class Pairs(_Arcs):
    __doc__ = _Arcs.__doc__
    DEFAULT_COLOR = '#dc9732'
    DEFAULT_ALPHA = 0.8
    fields = ["name", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]

    def __init__(self, file_, **kwargs):
        super().__init__(file_, **kwargs)
        self.bgz_file = process_pairs(file_)

    def convert_type(self, df, *args):
        df['pos1'] = df['pos1'].astype(int)
        df['pos2'] = df['pos2'].astype(int)
        return df

    def plot(self, ax, region_chrom, region_start, region_end):
        """
        Makes and arc connecting two points on a linear scale representing
        interactions between bins.
        """
        self.ax = ax
        gr = GenomeRange(region_chrom, region_start, region_end)
        itv_df = self.fetch_intervals(gr)
        intervals = []
        for _, row in itv_df.iterrows():
            s, e = row['pos1'], row['pos2']
            s, e = sorted([s, e])
            intervals.append((s, e, 1))
        self.plot_arcs(ax, gr, intervals)
        self.plot_label()


def Arcs(file_, *args, **kwargs):
    """Compose Arcs track(.bedpe, .pairs), determine type by file extension."""
    if file_.endswith(".bedpe") or file_.endswith('.bedpe.bgz'):
        return BEDPE(file_, *args, **kwargs)
    elif file_.endswith(".pairs") or file_.endswith('.pairs.bgz'):
        return Pairs(file_, *args, **kwargs)
    else:
        raise NotImplementedError("Arcs track only support .bedpe or .pairs input format.")
