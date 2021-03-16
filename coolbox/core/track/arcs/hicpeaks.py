from .bedpe import BEDPE
from .base import ArcsBase


class HiCPeaks(BEDPE):
    """
    Hi-C Peaks(Loops) from .bedpe file. Used to show the peaks on the Hi-C interaction map.
    """
    DEFAULT_PROPERTIES = {
        "style": ArcsBase.STYLE_HICPEAKS,
        'color': "#2255ff",
        "open_region": True,
        "alpha": 0.6,
        "line_width": 5,
    }

    def __init__(self, file, **kwargs):
        properties = HiCPeaks.DEFAULT_PROPERTIES.copy()
        properties.update(kwargs)
        super().__init__(file, **properties)
