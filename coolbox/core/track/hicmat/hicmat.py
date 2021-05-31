from pathlib import Path
from typing import Union

from .base import HicMatBase
from .cool import Cool
from .dothic import DotHiC


def HiCMat(file: Union[str, HicMatBase], *args, **kwargs) -> HicMatBase:
    """
    Compose DotHic or Cool track automatically based on tpye of file extension (.cool, .mcool, .hic)
    """
    if isinstance(file, HicMatBase):
        return file
    elif not Path(file).is_file():
        raise ValueError("The file path does not exist.")

    if file.endswith(".hic"):
        return DotHiC(file, *args, **kwargs)
    elif file.endswith((".cool", ".mcool")):
        return Cool(file, *args, **kwargs)
    else:
        raise NotImplementedError(f"File type of {p} not supported for HicMat. "
                                  f"The file type should be one of .cool/.mcool/.hic")
