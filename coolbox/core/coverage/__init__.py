from .highlights import HighLights, HighLightsFromFile
from .vlines import Vlines, VlinesFromFile
from .hlines import HLines

from .base import track_to_coverage
from ..track.bed import TAD
from ..track.hist import Hist, BigWig, BedGraph
from ..track.arcs import Arcs, Pairs, BEDPE, HiCPeaks

TADCoverage = track_to_coverage(TAD)

ArcsCoverage = track_to_coverage(Arcs)
PairsCoverage = track_to_coverage(Pairs)
BEDPECoverage = track_to_coverage(BEDPE)
HiCPeaksCoverage = track_to_coverage(HiCPeaks)

HistCoverage = track_to_coverage(Hist)
BigWigCoverage = track_to_coverage(BigWig)
BedGraphCoverage = track_to_coverage(BedGraph)
