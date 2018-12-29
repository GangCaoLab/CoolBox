import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist

from coolbox.utilities import cm2inch, get_logger, GenomeRange
log = get_logger(__name__)


class PlotFrame(object):

    DEFAULT_WIDTH = 40
    DEFAULT_WIDTH_RATIOS = (0.01, 0.93, 0.06)
    DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0, 'top': 1}

    def __init__(self, properties_dict, tracks, *args, **kwargs):

        self.properties = properties_dict
        self.tracks = tracks

        if 'width' not in self.properties:
            self.properties['width'] = PlotFrame.DEFAULT_WIDTH

        if 'width_ratios' not in self.properties:
            self.properties['width_ratios'] = PlotFrame.DEFAULT_WIDTH_RATIOS

        if 'margins' not in self.properties:
            self.properties['margins'] = PlotFrame.DEFAULT_MARGINS

        super().__init__()

    def get_tracks_height(self, default_height=3):
        """
        Get heights of all tracks.

        Return
        ------
        heights : list of float
            heights of all tracks.
        """
        heights = []
        for track in self.tracks.values():
            if hasattr(track, 'get_track_height'):
                frame_width = self.properties['width'] * self.properties['width_ratios'][1]
                height = track.get_track_height(frame_width)
                heights.append(height)
            elif 'height' in track.properties:
                heights.append(track.properties['height'])
            else:
                heights.append(default_height)
        return heights

    def plot(self, *args):
        """
        Plot all tracks.

        >>> from coolbox.api import *
        >>> frame = XAxis() + XAxis()
        >>> frame.plot("chr1", 100000, 200000)
        >>> frame.plot("chr1:100000-200000)
        """
        if len(args) >= 3:
            chrom, start, end = args[:3]
        else:
            region_str = args[0]
            gr = GenomeRange(region_str)
            chrom, start, end = (gr.chrom, gr.start, gr.end)

        tracks_height = self.get_tracks_height()
        self.properties['height'] = sum(tracks_height)

        fig = plt.figure(figsize=cm2inch(self.properties['width'],
                                         self.properties['height']))
        if 'title' in self.properties:
            fig.suptitle(self.properties['title'])

        grids = matplotlib.gridspec.GridSpec(
            len(tracks_height), 3,
            height_ratios=tracks_height,
            width_ratios=self.properties['width_ratios'],
            wspace=0.01)

        axis_list = []
        for idx, track in enumerate(self.tracks.values()):
            y_ax = plt.subplot(grids[idx, 0])
            y_ax.set_axis_off()

            ax = axisartist.Subplot(fig, grids[idx, 1])
            fig.add_subplot(ax)
            ax.axis[:].set_visible(False)
            ax.patch.set_visible(False)

            label_ax = plt.subplot(grids[idx, 2])
            label_ax.set_axis_off()

            track.label_ax = label_ax
            track.y_ax = y_ax
            try:
                track.plot(ax, chrom, start, end)

            except Exception as e:
                import sys, os
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = exc_tb.tb_frame.f_code.co_filename
                log.error("Error occured when plot track:\n"
                          "\ttrack name: {}\n\ttrack type:{}\n"
                          "\tError: {} {}\n"
                          "\toccurred in \"{}\", line {}".format(
                    track.name, type(track),
                    type(e), str(e), fname, exc_tb.tb_lineno)
                )
                log.exception(e)
            # plot coverages
            if hasattr(track, 'coverages'):
                for cov_idx, cov in enumerate(track.coverages):
                    cov.track = track
                    try:
                        cov.plot(ax, chrom, start, end)
                    except Exception as e:
                        log.error("Error occured when plot track's coverage:\n"
                                  "\ttrack name: {}\n\ttrack type:{}\n\tcoverage name: {}\n\tcov type: {}\n"
                                  "\tError: {} {}".format(
                            track.name, type(track), cov.name, type(cov),
                            type(e), str(e)))
                        log.exception(e)

            axis_list.append(ax)

        margins = self.properties['margins']
        fig.subplots_adjust(wspace=0, hspace=0.0,
                            left=margins['left'],
                            right=margins['right'],
                            bottom=margins['bottom'],
                            top=margins['top'])

        plt.close()

        return fig
