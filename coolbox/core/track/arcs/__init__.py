from .base import ArcsBase
from .bedpe import BEDPE
from .hicpeaks import HiCPeaks
from .pairs import Pairs


def Arcs(file_, *args, **kwargs) -> ArcsBase:
    """Compose BEDPE, Pairs track automatically based on type of file extension.(.bedpe, .pairs)
    """
    if file_.endswith((".bedpe", ".bedpe.bgz")):
        return BEDPE(file_, *args, **kwargs)
    elif file_.endswith((".pairs", ".pairs.bgz")):
        return Pairs(file_, *args, **kwargs)
    else:
        raise NotImplementedError("Arcs track only support .bedpe or .pairs input format.")
