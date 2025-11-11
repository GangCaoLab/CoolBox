from .base import ArcsBase
from coolbox.utilities.genome import GenomeRange


class Pairs(ArcsBase):
    """
    Arcs track from .pairs file.

    Parameters
    ----------
    file: str
        Path of .pairs file

    """
    DEFAULT_PROPERTIES = {
        'color': '#dc9732'
    }
    FIELDS = ["name", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]

    def __init__(self, file, **kwargs):
        properties = Pairs.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(**properties)