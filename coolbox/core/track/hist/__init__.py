from .bedgraph import BedGraph
from .bigwig import BigWig, ABCompartment


def Hist(file, **kwargs):
    """Hist track(bigwig, bedgraph), determine track type by file extension."""
    if file.endswith(".bw") or file.endswith(".bigwig"):
        return BigWig(file, **kwargs)
    elif file.endswith('.bedgraph') or file.endswith('.bg') or \
            file.endswith('.bedgraph.bgz') or file.endswith('.bg.bgz'):
        return BedGraph(file, **kwargs)
    else:
        raise NotImplementedError("Hist only support .bigwig or .bedgraph file now.")
