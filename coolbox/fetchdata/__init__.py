"""
Mixin classes,
API for fetch the data of tracks.
"""

from collections import OrderedDict

from coolbox.utilities import GenomeRange

from .bigwig import FetchBigWig
from .bedgraph import FetchBedGraph
from .bed import FetchBed
from .arcs import FetchBEDPE, FetchPairs
from .cool import FetchCool
from .dothic import FetchDotHiC
from .v4c import FetchVirtual4C
from .gtf import FetchGTF
from .bam import FetchBAM


__all__ = [
    "FetchBigWig", "FetchBedGraph", "FetchBed",
    "FetchBEDPE", "FetchPairs",
    "FetchCool", "FetchDotHiC",
    "FetchVirtual4C", "FetchGTF",
    "FetchFrame", "FetchBAM"
]


class FetchFrame(object):
    """
    FetchFrame used for fetch data from all tracks in the frame.
    """

    def fetch_data(self, genome_range=None):
        if genome_range is None:
            genome_range = self.current_range

        if isinstance(genome_range, str):
            genome_range = GenomeRange(genome_range)

        tracks_data = OrderedDict()
        for name, track in self.tracks.items():
            if hasattr(track, 'fetch_data'):
                data = track.fetch_data(genome_range)
            else:
                data = []
            tracks_data.update([(name, data)])

        return tracks_data


