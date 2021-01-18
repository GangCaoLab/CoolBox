from coolbox.core.track.bed.fetch import FetchBed
from coolbox.utilities import (
    get_logger
)
from coolbox.utilities.genome import GenomeRange
from coolbox.utilities.bed import build_bed_index
from .base import BedBase

log = get_logger(__name__)


class BED(BedBase, FetchBed):
    """
    Bed Track for plotting 1d intervals data from .bed file.
    The input bed file can be bed3/bed6/bed9/bed12

    Parameters
    ----------
    file: str
        The file path of `.bed` file.


    """

    DEFAULT_PROPERTIES = {
        'labels': "on",
    }

    def __init__(self, file, **kwargs):
        properties = BED.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(**properties)
        self.bgz_file = build_bed_index(file)

    def fetch_data(self, gr: GenomeRange, **kwargs):
        intervals = self.fetch_intervals(self.bgz_file, gr)
        return intervals


