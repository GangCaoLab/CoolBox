import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist

from coolbox.utilities import cm2inch, get_logger
log = get_logger(__name__)


class PlotFrame(object):

    DEFAULT_WIDTH = 40
    DEFAULT_WIDTH_RATIOS = (0.93, 0.07)
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
            if 'height' in track.properties:

                # auto specify height for Cool Track
                if track.properties['height'] == 'cool_auto':
                    cool_height = track.get_tracks_height(
                        self.properties['width'] * self.properties['width_ratios'][0])
                    heights.append(cool_height)
                else:
                    heights.append(track.properties['height'])

            else:
                heights.append(default_height)
        return heights

    def plot(self, chrom, start, end):
        """
        Plot all tracks.
        """

        tracks_height = self.get_tracks_height()
        self.properties['height'] = sum(tracks_height)

        fig = plt.figure(figsize=cm2inch(self.properties['width'],
                                         self.properties['height']))
        if 'title' in self.properties:
            fig.suptitle(self.properties['title'])

        grids = matplotlib.gridspec.GridSpec(
            len(tracks_height), 2,
            height_ratios=tracks_height,
            width_ratios=self.properties['width_ratios'])
        axis_list = []
        for idx, track in enumerate(self.tracks.values()):
            axis = axisartist.Subplot(fig, grids[idx, 0])
            fig.add_subplot(axis)
            axis.axis[:].set_visible(False)
            axis.patch.set_visible(False)
            label_axis = plt.subplot(grids[idx, 1])
            label_axis.set_axis_off()
            try:
                track.plot(axis, label_axis, chrom, start, end)

            except Exception as e:
                log.error("Error occured when plot track 'name: {}, type:{}', {}:{}".format(
                    track.name, type(track),
                    type(e), str(e)))
            # plot coverages
            if hasattr(track, 'coverages'):
                for cov_idx, cov in enumerate(track.coverages):
                    cov.track = track
                    try:
                        cov.plot(axis, chrom, start, end)
                    except Exception as e:
                        log.error("Error occured when plot track's coverage "
                                  "'track name: {}, track type:{}, cov idx: {}, cov type: {}', "
                                  "{}:{}".format(
                            track.name, type(track), cov_idx, type(cov),
                            type(e), str(e)))

            axis_list.append(axis)

        margins = self.properties['margins']
        fig.subplots_adjust(wspace=0, hspace=0.0,
                            left=margins['left'],
                            right=margins['right'],
                            bottom=margins['bottom'],
                            top=margins['top'])

        plt.close()

        return fig
