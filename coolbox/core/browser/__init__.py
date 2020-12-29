from .base import Browser
from coolbox.utilities import op_err_msg

__all__ = ["WidgetsPanel", "Browser"]


class WidgetsPanel(object):
    """
    Widgets container.

    Parameters
    ----------
    type : {'simple', 'full'}, optional
        Widgets type 'simple' or 'full',
        default 'simple'

    reference_genome : str, optional
        Reference genome,
        built-in references:('hg19', 'hg38', 'mm9', 'mm10')
        if you want use other genome, you can specify the "chromosome length file",
        that is a tab splited file, first column is the chromosomes,
        and second column is the length of correspond chromosome.
        default 'hg19'

    """

    def __init__(self, type="simple", reference_genome="hg19"):
        self.ref = reference_genome
        self.type = type

    def __add__(self, other):
        """

        >>>

        """
        from ..frame.base import FrameBase
        if isinstance(other, FrameBase):
            return Browser(other, reference_genome=self.ref, widgets_box=self.type)
        else:
            raise TypeError(op_err_msg(self, other, op="+"))


if __name__ == "__main__":
    import doctest

    doctest.testmod()
