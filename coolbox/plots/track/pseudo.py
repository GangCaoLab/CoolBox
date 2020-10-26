from coolbox.plots.track.base import TrackPlot


class PlotSpacer(TrackPlot):

    def plot(self, ax, chrom_region, start_region, end_region):
        self.ax = ax
        ax.set_xlim(start_region, end_region)


class PlotHLine(TrackPlot):

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)
        if 'line_style' not in self.properties:
            self.properties['line_style'] = '--'
        if 'line_width' not in self.properties:
            self.properties['line_width'] = 0.5
        if 'color' not in self.properties:
            self.properties['color'] = '#000000'
        if 'alpha' not in self.properties:
            self.properties['alpha'] = 0.75

    def plot(self, ax, chrom_region, region_start, region_end):
        self.ax = ax

        ax.set_xlim(region_start, region_end)
        ax.hlines(0, region_start, region_end,
                  linestyles=self.properties['line_style'],
                  linewidth=self.properties['line_width'],
                  colors=self.properties['color'],
                  alpha=self.properties['alpha'])


class PlotXAxis(TrackPlot):

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)
        if 'fontsize' not in self.properties:
            self.properties['fontsize'] = 15

    def plot(self, ax, chrom_region, region_start, region_end):
        self.ax = ax

        ax.set_xlim(region_start, region_end)
        ticks = ax.get_xticks()
        if ticks[-1] - ticks[1] <= 1e5:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"

        elif 1e5 < ticks[-1] - ticks[1] < 4e6:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"
        else:
            labels = ["{:,.1f} ".format((x / 1e6))
                      for x in ticks]
            labels[-2] += " Mbp"

        ax.axis["x"] = ax.new_floating_axis(0, 0.5)

        ax.axis["x"].axis.set_ticklabels(labels)
        ax.axis['x'].axis.set_tick_params(which='minor', bottom='on')

        ax.axis["x"].major_ticklabels.set(size=int(self.properties['fontsize']))

        if 'where' in self.properties and self.properties['where'] == 'top':
            ax.axis["x"].set_axis_direction("top")

