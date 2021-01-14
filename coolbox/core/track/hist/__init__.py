from .base import HistBase
from .bedgraph import BedGraph
from .bigwig import BigWig, ABCompartment
from .hicfeature import DiScore, InsuScore, Virtual4C
from .bam import BAMCov
from .snp import SNP


def Hist(file, *args, **kwargs) -> HistBase:
    """Compose BigWig, BedGraph, SNP, BAMCov track automatically based on type of file
    extension.(.bw, .bigwig, .snp, .vcf, .bedgraph, .bam, .sam)
    """
    if file.endswith((".bw", ".bigwig", ".bigWig")):
        return BigWig(file, *args, **kwargs)
    elif file.endswith(('.bedgraph', ".bg", ".bedgraph.bgz", ".bg.bgz")):
        return BedGraph(file, *args, **kwargs)
    elif file.endswith((".snp", ".vcf")):
        return SNP(file, *args, **kwargs)
    elif file.endswith(".bam", ".sam"):
        return BAMCov(file, *args, **kwargs)
    else:
        raise NotImplementedError("Hist only support .bigwig or .bedgraph file now.")
