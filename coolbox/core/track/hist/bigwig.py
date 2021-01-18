import numpy as np
import pandas as pd

from coolbox.utilities import (
    split_genome_range, change_chrom_names,
    GenomeRange, get_logger, to_gr
)
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
        import pyBigWig
        self.bw = pyBigWig.open(self.properties['file'])

    def fetch_plot_data(self, gr: GenomeRange, **kwargs):
        num_bins = self.get_num_bins()
        self.check_chrom_name(gr)
        scores_per_bin = self.fetch_scores(gr, num_bins)
        return scores_per_bin

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
        if chrom not in self.bw.chroms():
            chrom = change_chrom_names(chrom)

        intervals = self.bw.intervals(chrom, start, end)

        col_chrom = [chrom] * len(intervals)
        col_start = []
        col_end = []
        col_score = []
        for s, e, v in intervals:
            col_start.append(s)
            col_end.append(e)
            col_score.append(v)

        intval_table = pd.DataFrame(
            {
                "chromsome": col_chrom,
                "start": col_start,
                "end": col_end,
                "score": col_score,
            },
            columns=['chromsome', 'start', 'end', 'score']
        )

        return intval_table

    def get_num_bins(self, default_num=700):
        num_bins = default_num
        if 'number_of_bins' in self.properties:
            try:
                num_bins = int(self.properties['number_of_bins'])
            except TypeError:
                num_bins = default_num
                log.warning("'number_of_bins' value: {} for bigwig file {} "
                            "is not valid. Using default value (700)".format(self.properties['number_of_bins'],
                                                                             self.properties['file']))
        return num_bins

    def fetch_scores(self, genome_range, num_bins, max_try_nums=5):
        """
        Fetch bins scores within input chromosome range.

        on rare occasions pyBigWig may throw an error, apparently caused by a corruption
        of the memory. This only occurs when calling trackPlot from different
        processors. Reloading the file solves the problem.
        """
        num_tries = 0
        scores_per_bin = np.zeros(num_bins)
        gr = to_gr(genome_range)
        if gr.chrom not in self.bw.chroms():
            gr.change_chrom_names()
        while num_tries < max_try_nums:
            num_tries += 1
            try:
                scores_per_bin = np.array(
                    self.bw.stats(gr.chrom, gr.start, gr.end, nBins=num_bins)
                ).astype(float)
            except Exception as e:
                import pyBigWig
                self.bw = pyBigWig.open(self.properties['file'])

                log.warning("error found while reading bigwig scores ({}).\nTrying again. Iter num: {}".
                            format(e, num_tries))
                pass
            else:
                if num_tries > 1:
                    log.warning("After {} the scores could be computed".format(num_tries))
                break
        return scores_per_bin

    def check_chrom_name(self, genome_range):
        if genome_range.chrom not in self.bw.chroms().keys():
            genome_range.change_chrom_names()

        if genome_range.chrom not in self.bw.chroms().keys():
            log.warning("Can not read region {} from bigwig file:\n\n"
                        "{}\n\nPlease check that the chromosome name is part of the bigwig file "
                        "and that the region is valid".format(str(genome_range), self.properties['file']))


class ABCompartment(BigWig):
    """
    A/B Compartment BigWig track.
    """

    DEFAULT_PROPERTIES = {
        "color": "#66ccff",
        "threshold_color": "#ff9c9c",
        "threshold": 0,
    }

    def __init__(self, file, **kwargs):
        properties = ABCompartment.DEFAULT_PROPERTIES.copy()
        properties.update(kwargs)
        super().__init__(file, **properties)
