from .base import track_to_coverage
from ..track.hist import Hist, BigWig, BedGraph


HistCoverage = track_to_coverage(Hist)
BigWigCoverage = track_to_coverage(BigWig)
BedGraphCoverage = track_to_coverage(BedGraph)

