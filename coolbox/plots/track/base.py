class TrackPlot(object):
    """
    The TrackPlot object is a holder for all tracks that are to be plotted.
    For example, to plot a bedgraph file a new class that extends TrackPlot
    should be created.

    It is expected that all TrackPlot objects have a plot method.

    """

    def __init__(self, *args, **kwargs):
        if not hasattr(self, 'properties'):
            self.properties = args[0]
        super().__init__()

    def plot_y_axis(self, y_ax):
        """
        Plot the scale of the y axis with respect to the plot_axis

        plot something that looks like this:
        ymax ┐
             │
             │
        ymin ┘

        Parameters
        ----------
        y_ax : matplotlib.axes.Axes
            Axis to use to plot the scale
        """
        plot_axis = self.ax

        if 'show_data_range' in self.properties and self.properties['show_data_range'] == 'no':
            return

        def value_to_str(value):
            if value % 1 == 0:
                str_value = str(int(value))
            else:
                if value < 0.01:
                    str_value = "{:.4f}".format(value)
                else:
                    str_value = "{:.2f}".format(value)
            return str_value

        ymin, ymax = plot_axis.get_ylim()

        ymax_str = value_to_str(ymax)
        ymin_str = value_to_str(ymin)
        x_pos = [0, 0.5, 0.5, 0]
        y_pos = [0.01, 0.01, 0.99, 0.99]
        y_ax.plot(x_pos, y_pos, color='black', linewidth=1, transform=y_ax.transAxes)
        y_ax.text(-0.2, -0.01, ymin_str, verticalalignment='bottom', horizontalalignment='right', transform=y_ax.transAxes)
        y_ax.text(-0.2, 1, ymax_str, verticalalignment='top', horizontalalignment='right', transform=y_ax.transAxes)
        y_ax.patch.set_visible(False)

    def plot_label(self):
        if hasattr(self, 'label_ax'):
            self.label_ax.text(0.15, 0.5, self.properties['title'],
                               horizontalalignment='left', size='large',
                               verticalalignment='center')
