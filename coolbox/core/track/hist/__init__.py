from .bedgraph import BedGraph
from .bigwig import BigWig, ABCompartment
from .hicfeature import DiScore, InsuScore, Virtual4C
from .base import HistBase


def Hist(file, *args, **kwargs) -> HistBase:
    """Hist track(bigwig, bedgraph), determine track type by file extension."""
    if file.endswith(".bw") or \
            file.endswith(".bigwig"):
        return BigWig(file, *args, **kwargs)
    elif file.endswith('.bedgraph') or \
            file.endswith('.bg') or \
            file.endswith('.bedgraph.bgz') or \
            file.endswith('.bg.bgz'):
        return BedGraph(file, *args, **kwargs)
    else:
        raise NotImplementedError("Hist only support .bigwig or .bedgraph file now.")
