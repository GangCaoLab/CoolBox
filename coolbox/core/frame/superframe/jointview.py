from typing import Union

from collections import OrderedDict

import svgutils.compose as sc
import matplotlib.pyplot as plt

from coolbox.utilities.figtools import cm2inch
from coolbox.utilities.filetool import get_uniq_tmp_file
from coolbox.utilities import GenomeRange
from coolbox.core.track import Track
from coolbox.core.frame import Frame
from .base import SuperFrame


class JointView(SuperFrame):
    """Compose two track and a center matrix.

    https://github.com/GangCaoLab/CoolBox/issues/12

    Parameters
    ----------
    center : Track
        Center track for show 'contact map-like' plot,
        should support for plotting images with 2d GenomeRanges.

    top : Frame, optional
        Frame plot in the top of the center track.

    right : Frame, optional

    bottom : Frame, optional

    left : Frame, optional

    center_width : float
        The width of center track, unit in cm. default 20.

    trbl : str
        sub-frames(top right bottom left) use which genome range(first or second),
        default '1212', which means: top -> 1, right -> 2, bottom -> 1, left -> 2

    space : float
        Space between frame and center, unit in cm. default 0.5
    """

    def __init__(self,
                 center: Track,
                 top: Union[Frame, Track] = None,
                 right: Union[Frame, Track] = None,
                 bottom: Union[Frame, Track] = None,
                 left: Union[Frame, Track] = None,
                 **kwargs,
                 ):
        sub_frames = OrderedDict({
            "top": top,
            "right": right,
            "bottom": bottom,
            "left": left,
        })
        sub_frames = {k: v for k, v in sub_frames.items() if v is not None}
        self.__check_sub_frames(center, sub_frames)
        # explicitely use matrix style
        center.properties['style'] = 'matrix'
        properties = {
            "sub_frames": sub_frames,
            "center_track": center,
            "center_width": 20,
            "trbl": "1212",
            "space": 1,
            "cm2px": 28.5,
            "padding_left": 1,
        }
        properties.update(**kwargs)

        super().__init__(properties)
        self.__adjust_sub_frames_width()

    def cm2px(self, vec):
        return [i * self.properties['cm2px'] for i in vec]

    @staticmethod
    def __check_sub_frames(center, sub_frames):
        from ..frame import Frame
        from ...track.base import Track
        sub_f_names = ", ".join(sub_frames.keys())
        if (not isinstance(center, Track)) and (not hasattr(center, "plot_joint")):
            raise TypeError("center should be a Track type instance with plot_joint method, "
                            "for example Cool, DotHiC, ...")
        for k, f in sub_frames.items():
            if not isinstance(f, Frame):
                if isinstance(f, Track):
                    sub_frames[k] = (Frame() + f)  # convert track to frame
                else:
                    raise TypeError(f"{sub_f_names} should be Frame object.")

    def __adjust_sub_frames_width(self):
        for k, f in self.properties['sub_frames'].items():
            width_ratios = f.properties['width_ratios']
            middle_ratio = width_ratios[1]
            new_width = self.properties['center_width'] / middle_ratio
            f.properties['width'] = new_width

    def plot_center(self, gr1: GenomeRange, gr2: GenomeRange):
        center_track = self.properties['center_track']
        size = cm2inch(self.properties['center_width'])
        fig, ax = plt.subplots(figsize=(size, size))
        center_track.plot(ax, gr2, gr2=gr1)
        center_track.plot_coverages(ax, gr2, gr1)
        ax.set_axis_off()
        path = get_uniq_tmp_file(prefix='center', suffix='.svg')
        fig.subplots_adjust(wspace=0, hspace=0.0, left=0, right=1, bottom=0, top=1)
        fig.savefig(path)
        plt.close()
        return sc.SVG(path)

    def frame_granges(self, gr1=None, gr2=None):
        self.goto(gr1, gr2)
        gr1, gr2 = self.current_range

        trbl = self.properties['trbl']

        orientations = ['top', 'right', 'bottom', 'left']
        frame2grange = {
            k: (gr1, gr2) if (trbl[orientations.index(k)] == '1') else (gr2, gr1)
            for k in orientations
        }

        return frame2grange

    def plot(self, gr1=None, gr2=None):
        """

        Parameters
        ----------
        gr1 : {str, GenomeRange}
            First genome range

        gr2 : {str, GenomeRange}, optional
            Second genome range
        """
        frame2grange = self.frame_granges(gr1, gr2)
        gr1, gr2 = self.current_range
        sub_frames = self.properties['sub_frames']

        frame_svgs = self.plot_frames(frame2grange)
        center_svg = self.plot_center(gr1, gr2)

        center_offsets = self.__get_center_offsets(sub_frames)

        center_svg.move(*self.cm2px(center_offsets))
        self.__transform_sub_svgs(frame_svgs, sub_frames, center_offsets)

        figsize = self.cm2px(self.__get_figsize(sub_frames))
        fig = sc.Figure(f"{figsize[0]}px", f"{figsize[1]}px",
                        sc.Panel(center_svg),
                        *[sc.Panel(svg) for svg in frame_svgs.values()])
        return fig

    def fetch_data(self, gr1=None, gr2=None) -> dict:
        """

        Parameters
        ----------
        gr1 : {str, GenomeRange}, optional
            First genome range

        gr2 : {str, GenomeRange}, optional
            Second genome range
        """
        tracks_data = {}
        frame2grange = self.frame_granges(gr1, gr2)
        gr1, gr2 = self.current_range
        sub_frames = self.properties['sub_frames']

        for pos, fr in sub_frames.items():
            tracks_data[pos] = sub_frames[pos].fetch_data(frame2grange[pos][0], gr2=frame2grange[pos][1])

        tracks_data['center'] = self.properties['center_track'].fetch_data(gr1, gr2=gr2)

        return tracks_data

    def __transform_sub_svgs(self, sub_svgs, sub_frames, center_offsets):
        c_width = self.properties['center_width']
        space = self.properties['space']
        pd_left = self.properties['padding_left']
        if 'top' in sub_svgs:
            s = sub_svgs['top']
            offsets = [pd_left, 0]
            if 'left' in sub_svgs:
                offsets[0] += sub_frames['left'].properties['height'] + space
            s.move(*self.cm2px(offsets))
        if 'right' in sub_svgs:
            f = sub_frames['right']
            s = sub_svgs['right']
            s.rotate(90)
            wr = f.properties['width_ratios']
            right_offsets = [
                center_offsets[0] + c_width + f.properties['height'] + space,
                center_offsets[1] - f.properties['width'] * wr[0]
            ]
            right_offsets = self.cm2px(right_offsets)
            s.move(*right_offsets)
        if 'bottom' in sub_svgs:
            s = sub_svgs['bottom']
            offsets = [pd_left, c_width]
            if 'left' in sub_svgs:
                offsets[0] += sub_frames['left'].properties['height'] + space
            if 'top' in sub_svgs:
                offsets[1] += sub_frames['top'].properties['height']
            offsets = self.cm2px(offsets)
            s.move(*offsets)
        if 'left' in sub_svgs:
            s = sub_svgs['left']
            f = sub_frames['left']
            s.rotate(90)
            offsets = [pd_left + f.properties['height'], 0]
            if 'top' in sub_svgs:
                offsets[1] += sub_frames['top'].properties['height']
            offsets = self.cm2px(offsets)
            s.move(*offsets)

    def __get_center_offsets(self, sub_frames):
        space = self.properties['space']
        center_offsets = [self.properties['padding_left'], 0]  # x, y (left, top)
        if 'top' in sub_frames:
            f = sub_frames['top']
            center_offsets[0] += f.properties['width'] * f.properties['width_ratios'][0]
            center_offsets[1] += space + f.properties['height']
        if 'bottom' in sub_frames:
            f = sub_frames['bottom']
            if 'top' not in sub_frames:
                center_offsets[0] += f.properties['width'] * f.properties['width_ratios'][0]
        if 'left' in sub_frames:
            center_offsets[0] += space + sub_frames['left'].properties['height']
        return center_offsets

    def __get_figsize(self, sub_frames):
        space = self.properties['space']
        center_width = self.properties['center_width']
        size = [center_width, center_width]  # width, height
        if 'top' in sub_frames:
            f = sub_frames['top']
            size[0] = f.properties['width']
            size[1] += f.properties['height'] + space
        if 'bottom' in sub_frames:
            f = sub_frames['bottom']
            size[0] = f.properties['width']
            size[1] += f.properties['height'] + space
        if 'left' in sub_frames:
            f = sub_frames['left']
            size[0] += f.properties['height'] + space
        if 'right' in sub_frames:
            f = sub_frames['right']
            size[0] += f.properties['height'] + space
        self.properties['width'] = size[0]
        self.properties['height'] = size[1]
        return size

    def add_track(self, track, pos=None):
        sub_frames = self.properties['sub_frames']
        if pos is None:
            pos = list(sub_frames.keys())[-1]
        frame = sub_frames[pos]
        frame.add_track(track)
