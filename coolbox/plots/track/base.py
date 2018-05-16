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
