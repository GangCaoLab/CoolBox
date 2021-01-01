import abc
from collections import OrderedDict

import svgutils.compose as sc
from IPython.display import SVG

from coolbox.utilities.filetool import get_uniq_tmp_file
from ..base import FrameBase


class SuperFrame(FrameBase, abc.ABC):
    """Bigger frame composed by normal frames,
    compose figure using svgutils,
    this allow compose big figures which the matplotlib can not do.
    """

    def __init__(self, properties_dict, **kwargs):
        assert "sub_frames" in properties_dict
        super().__init__(properties_dict, **kwargs)

    def plot_frames(self, frame2grange):
        """Plot each frame by a given GenomeRange object."""
        res = OrderedDict()
        for k, f in self.properties['sub_frames'].items():
            gr = frame2grange[k]
            path = get_uniq_tmp_file(prefix="frame_", suffix=".svg")
            fig = f.plot(gr)
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

    def show(self):
        if isinstance(self.current_range, list):
            svg = self.plot(*self.current_range)
        else:
            svg = self.plot(self.current_range)
        path = get_uniq_tmp_file(prefix="sub_frame", suffix=".svg")
        svg.save(path)
        return SVG(path)
