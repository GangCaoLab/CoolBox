import abc
from collections import OrderedDict

import svgutils.compose as sc
from IPython.display import SVG

from coolbox.utilities.filetool import get_uniq_tmp_file
from coolbox.utilities import GenomeRange
from ..base import FrameBase


class SuperFrame(FrameBase, abc.ABC):
    """Bigger frame composed by normal frames,
    compose figure using svgutils,
    this allow compose big figures which the matplotlib can not do.
    """

    def __init__(self, properties_dict, **kwargs):
        assert "sub_frames" in properties_dict
        super().__init__(properties_dict, **kwargs)
        self.current_range = [None, None]

    def plot_frames(self, frame2grange):
        """Plot each frame by a given GenomeRange object."""
        res = OrderedDict()
        for k, f in self.properties['sub_frames'].items():
            gr, gr2 = frame2grange[k]
            path = get_uniq_tmp_file(prefix="frame_", suffix=".svg")
            fig = f.plot(gr, gr2)
            fig.subplots_adjust(wspace=0, hspace=0.0, left=0, right=1, bottom=0, top=1)
            fig.savefig(path)
            svg = sc.SVG(path)
            res[k] = svg
        return res

    @property
    def tracks(self):
        sub_frames = self.properties['sub_frames']
        tracks = OrderedDict()
        for f in sub_frames.values():
            tracks.update(f.tracks)
        return tracks

    def goto(self, gr1=None, gr2=None):
        if gr1 is not None:
            gr1 = GenomeRange(gr1)
        if gr2 is not None:
            gr2 = GenomeRange(gr2)
        if gr1 is None:
            gr1 = self.current_range[0]
        if gr2 is None:
            gr2 = gr1

        if gr1 is None or gr2 is None:
            raise ValueError("No history gr found.")
        self.current_range = [gr1, gr2]

    def plot(self, gr1=None, gr2=None):
        """

        Parameters
        ----------
        gr1 : {str, GenomeRange}
            First genome range

        gr2 : {str, GenomeRange}, optional
            Second genome range
        """
        raise NotImplementedError

    def show(self):
        gr1, gr2 = self.current_range
        svg = self.plot(gr1, gr2=gr2)
        path = get_uniq_tmp_file(prefix="sub_frame", suffix=".svg")
        svg.save(path)
        return SVG(path)
