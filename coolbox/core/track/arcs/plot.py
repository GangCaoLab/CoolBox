import numpy as np
import pandas as pd
from matplotlib.patches import Arc, Rectangle, Ellipse
from matplotlib import cm

from coolbox.utilities.genome import GenomeRange, get_logger


log = get_logger(__name__)


class PlotContacts(object):

    def plot_contacts(self, ax, gr: GenomeRange, gr2: GenomeRange, intervals: pd.DataFrame):
        style = self.properties['style']
        columns = intervals.columns
        if style == self.STYLE_ARCS:
            assert len(intervals) == 0 or all(c in columns for c in ['pos1', 'pos2']), \
                "The 'arcs' style's input DataFrame should have columns ['pos1', 'pos2']"
            self.plot_arcs(ax, gr, gr2, intervals)
        elif style == self.STYLE_HICPEAKS:
            assert len(intervals) == 0 or all(c in columns for c in ['start1', 'end1', 'start2', 'end2']), \
                "The 'hicpeaks' style's input DataFrame should have columns:" \
                " ['start1', 'end1', 'start2', 'end2'] or ['pos1', 'pos2']"
            self.plot_hicpeaks(ax, gr, gr2, intervals)
        else:
            raise ValueError("Style should be one of ['arcs', 'peaks']")

    # only support for 1 gr
    def plot_arcs(self, ax, gr: GenomeRange, gr2: GenomeRange, intervals: pd.DataFrame):
        """
        Parameters
        ----------
        ax
        gr
        gr2
        intervals : pd.DataFrame
            Should be with columns ['pos1', 'pos2', 'score']. The `score` column is optional
        """
        properties = self.properties

        def get_height(diameter):
            key = 'diameter_to_height'
            if key in properties:
                try:
                    return eval(properties[key])
                except Exception:
                    pass
            return 0.97 * max_height * diameter / max_diameter

        def get_linewidth(score):
            if ('line_width' in properties) and properties['line_width'] is not None:
                return float(properties['line_width'])
            elif 'score_to_width' in properties:
                try:
                    import numpy
                    import math
                    return eval(properties['score_to_width'])
                except Exception as e:
                    log.warning(str(e))
                    log.warning("Fail to evaluate width according to score")
            return 0.5 * np.sqrt(score)

        def adjust_yaxis(height):
            if properties['orientation'] == 'inverted':
                ax.set_ylim(height, 0.001)
            else:
                ax.set_ylim(-0.001, height)

        def get_color(score, min_score, max_score):
            if properties['cmap'] is not None:
                cmap = cm.get_cmap(properties['cmap'])
                return cmap((max(score, min_score)-min_score)/(max_score-min_score))
            else:
                return properties['color']

        alpha = properties['alpha']
        max_height = properties['height']
        adjust_yaxis(max_height)
        ax.set_xlim(gr.start, gr.end)
        if len(intervals) == 0:
            return

        if properties['cmap'] is not None:
            if properties['vmax'] is not None:
                max_score = properties['vmax']
            elif 'score' in intervals.columns:
                max_score = intervals['score'].max()
            else:
                max_score = 1
        else:
            max_score = 1
        min_score = properties['vmin'] or 0

        max_diameter = (intervals['pos2'] - intervals['pos1']).max()
        for row in intervals.itertuples():
            start, end = row.pos1, row.pos2
            score = row.score if hasattr(row, 'score') else 1
            line_width = get_linewidth(score)
            color = get_color(score, min_score, max_score)
            diameter = (end - start)
            height = 2 * get_height(diameter)
            center = (start + end) / 2
            ax.plot([center], [diameter])
            arc = Arc(
                (center, 0), diameter,
                height, 0, 0, 180,
                color=color,
                alpha=alpha,
                lw=line_width,
                linestyle=properties['line_style']
            )
            if properties['fill'] == 'yes':
                fill = Ellipse((center, 0), diameter, height,
                               alpha=properties['fill_alpha'],
                               color=properties['fill_color'],
                               fill=True)
                ax.add_patch(fill)
            ax.add_patch(arc)

    def plot_hicpeaks(self, ax, gr: GenomeRange, gr2: GenomeRange, intervals: pd.DataFrame):
        """ plot hic peaks uppon a HicMatBase track.

        Parameters
        ----------
        ax
        gr
        gr2
        intervals : pd.DatatFrame
            Should be with columns ['start1', 'end1', 'start2', 'end2']. The `score` column is optional
        """
        import math
        from coolbox.core.track.hicmat import HicMatBase
        # coverage only
        assert 'track' in self.__dict__ and isinstance(self.track, HicMatBase), \
            "No parent HicMatBase Track found. The hicpeaks style " \
            "must be plotted on a HicMatBase Track as a coverage track"
        hicmat_tri_style = (HicMatBase.STYLE_WINDOW, HicMatBase.STYLE_TRIANGULAR)
        hicmat_ma_style = (HicMatBase.STYLE_MATRIX,)

        if len(intervals) == 0:
            return

        # Calculate the position and size of the box from loop's start and end postion.
        def pos_triangular(st1, ed1, st2, ed2):
            m1 = (st1 + ed1) / 2
            m2 = (st2 + ed2) / 2
            x = (m1 + m2) / 2
            y = x - m1
            w = ((ed2 - st2) / 2) / math.cos(math.pi / 4)
            h = ((ed1 - st1) / 2) / math.cos(math.pi / 4)
            return x, y, (w, h)

        def pos_matrix_upper(st1, ed1, st2, ed2):
            x = (st2 + ed2) / 2
            y = (st1 + ed1) / 2
            w = ed2 - st2
            h = ed1 - st1
            return x, y, (w, h)

        def pos_matrix_lower(st1, ed1, st2, ed2):
            x = (st1 + ed1) / 2
            y = (st2 + ed2) / 2
            w = ed1 - st1
            h = ed2 - st2
            return x, y, (w, h)

        hictrack = self.track
        hicmat_style = hictrack.properties['style']
        hicmat_depth_ratio = hictrack.properties['depth_ratio']
        hicmat_depth_ratio = hicmat_depth_ratio if hicmat_depth_ratio != HicMatBase.DEPTH_FULL else 1
        depth_limit = hicmat_depth_ratio * (gr.end - gr.start) / 2

        properties = self.properties
        plot_kwargs = {
            'edgecolor': properties['color'],
            'facecolor': properties['fill_color'],
            'linewidth': properties['line_width'],
            'linestyle': properties['line_style'],
            'fill': True if properties['fill'] == 'yes' else False
        }
        side = properties['side']

        def set_alpha(rec):
            fc = list(rec.get_facecolor())
            fc[-1] = properties['fill_alpha']
            rec.set_facecolor(fc)
            ec = list(rec.get_edgecolor())
            ec[-1] = properties['alpha']
            rec.set_edgecolor(ec)
            return rec

        for loop in intervals.itertuples():
            (st1, ed1, st2, ed2) = loop.start1, loop.end1, loop.start2, loop.end2
            if hicmat_style in hicmat_tri_style:
                x, y, (w, h) = pos_triangular(st1, ed1, st2, ed2)
                if y >= depth_limit:
                    continue
                rec = Rectangle((x, y), w, h, angle=45, **plot_kwargs)
                rec = set_alpha(rec)
                ax.add_patch(rec)

            elif hicmat_style in hicmat_ma_style:
                if side in ('upper', 'both'):
                    # plot upper rectangle
                    x, y, (w, h) = pos_matrix_upper(st1, ed1, st2, ed2)
                    rec = Rectangle((x, y), w, h, **plot_kwargs)
                    rec = set_alpha(rec)
                    ax.add_patch(rec)

                if side in ('lower', 'both'):
                    # plot lower rectangle
                    x, y, (w, h) = pos_matrix_lower(st1, ed1, st2, ed2)
                    rec = Rectangle((x, y), w, h, **plot_kwargs)
                    rec = set_alpha(rec)
                    ax.add_patch(rec)
