from .base import track_to_coverage
from ..track.arcs import Arcs, Pairs, BEDPE


ArcsCoverage = track_to_coverage(Arcs)
PairsCoverage = track_to_coverage(Pairs)
BEDPECoverage = track_to_coverage(BEDPE)
