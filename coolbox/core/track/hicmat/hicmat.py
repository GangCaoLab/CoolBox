from pathlib import Path
from typing import Union

from .base import HicMatBase
from .cool import Cool
from .dothic import DotHiC


def HiCMat(file_or_hicmat: Union[str, HicMatBase], *args, **kwargs) -> HicMatBase:
    """
    Compose DotHic or Cool track automatically based on tpye of file extension (.cool, .mcool, .hic)
    """
    if isinstance(file_or_hicmat, HicMatBase):
        return file_or_hicmat
    elif not Path(file_or_hicmat).is_file():
        raise ValueError("The file path does not exist.")

    p = file_or_hicmat
    if p.endswith(".hic"):
        return DotHiC(file_or_hicmat, *args, **kwargs)
    elif p.endswith((".cool", ".mcool")):
        return Cool(file_or_hicmat, *args, **kwargs)
    else:
        raise NotImplementedError(f"File type of {p} not supported for HicMat. "
                                  f"The file type should be one of .cool/.mcool/.hic")
