from coolbox.utilities.hic.tools import hicmat_filetype
from .base import HicMatBase
from .cool import Cool
from .dothic import DotHiC


def HiCMat(file_, *args, **kwargs) -> HicMatBase:
    """Compose Hi-C track(.cool, .mcool, .hic), determine type by file extension."""
    ftype = hicmat_filetype(file_)
    if ftype == ".hic":
        return DotHiC(file_, *args, **kwargs)
    else:
        return Cool(file_, *args, **kwargs)
