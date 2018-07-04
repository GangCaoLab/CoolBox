from coolbox.utilities import op_err_msg
from coolbox.interact import BrowserBase


__all__ = ["WidgetsPanel", "Browser"]


class WidgetsPanel(object):
    """
    Widgets container.
    """

    def __init__(self, reference_genome="hg19"):
        self.ref = reference_genome

    def __add__(self, other):
        """

        >>>

        """
        from .frame import Frame
        if isinstance(other, Frame):
            return Browser(other, reference_genome=self.ref)
        else:
            raise TypeError(op_err_msg(self, other, op="+"))


class Browser(BrowserBase):
    """
    Genoeme browser.
    include:
        * Frame: for show plots
        * widgetsPanel: for control genome region showed in Frame.
    """

    @property
    def tracks(self):
        return self.frame.tracks

    def fetch_data(self, genome_range=None):
        return self.frame.fetch_data(genome_range)


if __name__ == "__main__":
    import doctest
    doctest.testmod()

