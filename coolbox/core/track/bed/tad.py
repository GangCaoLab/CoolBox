from .bed import BED
from .base import BedBase


class TAD(BED):
    """
    Tad tack from bed file
    """

    DEFAULT_PROPERTIES = {
        'style': BedBase.STYLE_TAD,
        'alpha': 0.3,
    }

    def __init__(self, file, **kwargs):
        properties = TAD.DEFAULT_PROPERTIES.copy()
        properties.update(kwargs)
        super().__init__(file, **properties)
