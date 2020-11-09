from .cool import Cool
from .dothic import DotHiC


def HiCMat(file_, *args, **kwargs):
    """Compose Hi-C track(.cool, .mcool, .hic), determine type by file extension."""
    if file_.endswith(".hic"):
        return DotHiC(file_, *args, **kwargs)
    elif file_.endswith(".cool") or file_.endswith(".mcool"):
        return Cool(file_, *args, **kwargs)
    else:
        raise NotImplementedError("Hi-C Matrix only support .hic or .cool input format.")
