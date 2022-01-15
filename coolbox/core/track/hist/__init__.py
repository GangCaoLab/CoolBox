from .base import HistBase
from .bedgraph import BedGraph
from .bigwig import BigWig
from .hicfeature import DiScore, InsuScore, Virtual4C
from .bam import BAMCov
from .snp import SNP


def Hist(file, *args, **kwargs) -> HistBase:
    """Compose BigWig, BedGraph, SNP, BAMCov track automatically based on type of file
    extension.(.bw, .bigwig, .snp, .vcf, .bedgraph, .bam, .sam)
    """
    if file.endswith((".bw", ".bigwig", ".bigWig")):
        return BigWig(file, *args, **kwargs)
    elif file.endswith(('.bedgraph', ".bg", ".bedgraph.bgz", ".bg.bgz", ".bedGraph", "bedGraph.bgz")):
        return BedGraph(file, *args, **kwargs)
    elif file.endswith((".snp", ".vcf")):
        return SNP(file, *args, **kwargs)
    elif file.endswith((".bam", ".sam")):
        return BAMCov(file, *args, **kwargs)
    else:
        raise NotImplementedError("Hist only support .bigwig or .bedgraph file now.")


def ABCompartment(file, *args, **kwargs):
    kwargs['threshold'] = 0
    kwargs['color'] = '#0000ff'
    kwargs['threshold_color'] = '#ff0000'
    return Hist(file, *args, **kwargs)
