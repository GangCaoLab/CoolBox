class CoveragePlot(object):

    """
    Coverage plot holder, similar to TrackPlot.
    """

    def __init__(self, *args, **kwargs):
        if not hasattr(self, 'properties'):
            self.properties = args[0]
        super().__init__()
