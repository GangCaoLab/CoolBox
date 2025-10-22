import numpy as np

from coolbox.utilities import (
    split_genome_range, change_chrom_names,
    GenomeRange, get_logger, to_gr
)
import oxbow as ox
from .base import HistBase

log = get_logger(__name__)


class BigWig(HistBase):
    """
    BigWig track

    Parameters
    -----------
    file : str
        File path of bigwig file.

    num_bins : int, optional
        Number of bins to plot the hist in current range, default 700.


    """

    DEFAULT_PROPERTIES = {
        "color": "#dfccde",
        "style": HistBase.STYLE_FILL,
        "num_bins": 700,
        "threshold": "inf"
    }

    def __init__(self, file, **kwargs):
        properties = BigWig.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(**properties)
        self.ds = ox.from_bigwig(self.properties['file'])

    def fetch_plot_data(self, gr: GenomeRange, **kwargs):
        self.check_chrom_name(gr)
        intervals = self.fetch_data(gr)
        values = intervals['value'].values
        return values

    def fetch_data(self, gr: GenomeRange, **kwargs):
        """
        Parameters
        ----------
        gr : GenomeRange

        Return
        ------
        intervals : pandas.core.frame.DataFrame
            BigWig interval table.
        """
        chrom, start, end = split_genome_range(gr)
        if chrom not in self.ds.chrom_names:
            chrom = change_chrom_names(chrom)

        intervals = self.ds.regions(f"{chrom}:{start}-{end}").pd()
        return intervals

    def check_chrom_name(self, genome_range):
        if genome_range.chrom not in self.ds.chrom_names:
            genome_range = genome_range.change_chrom_names()

        if genome_range.chrom not in self.ds.chrom_names:
            log.warning("Can not read region {} from bigwig file:\n\n"
                        "{}\n\nPlease check that the chromosome name is part of the bigwig file "
                        "and that the region is valid".format(str(genome_range), self.properties['file']))

