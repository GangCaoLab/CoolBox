class FetchTrackData(object):
    """
    FetchTrackData base class.
    FetchTrackData object used for fetch data from specify data format.

    All `FetchTrackData` must have a `fetch_data` method.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)