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
    BigWig track.

    Parameters
    ----------
    file_ : str
        Path to bigwig file.

    height : float, optional
        Height of track, default BigWig.DEFAULT_HEIGHT.

    color : str, optional
        Track color, default BigWig.DEFAULT_COLOR.

    positive_color : str, optional
        Track's positive value color

    negative_color : str, optional
        Track's negative value color

    alpha : float, optional
        Alpha value of plot, default 1.0

    number_of_bins : int, optional
        Number_of_bins in current range, default 700.

    style : str, optional
        Track graph type, format {'fill', 'line:`size`', 'points:`size`'},
        example: 'line:2', 'points:0.5'. default: 'fill'

    orientation : str, optional
        Track orientation, use 'inverted' for inverted track plot.

    show_data_range : bool, optional
        Show_data_range or not, default True.

    data_range_style : {'text', 'y-axis'}
        The style of the data range. default: 'y-axis'

    title : str, optional
        Label text, default ''

    max_value : {float, 'auto'}, optional
        Max value of track. 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of track. 'auto' for specify max value automatically, default 'auto'.

    name : str, optional
        Track's name.
    """

    DEFAULT_COLOR = "#dfccde"

    def __init__(self, file_, **kwargs):
        properties_dict = {
            'file': file_,
            'color': self.DEFAULT_COLOR,
            'alpha': 1.0,
            'number_of_bins': 700,
            'style': 'fill',
        }
        properties_dict.update(kwargs)
        super().__init__(**properties_dict)
        import pyBigWig
        self.bw = pyBigWig.open(self.properties['file'])

    def fetch_data(self, genome_range):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        Return
        ------
        intervals : pandas.core.frame.DataFrame
            BigWig interval table.
        """
        chrom, start, end = split_genome_range(genome_range)
        if chrom not in self.bw.chroms():
            chrom_ = change_chrom_names(chrom)
        else:
            chrom_ = chrom

        intervals = self.bw.intervals(chrom_, start, end)

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

    def fetch_plot_data(self, genome_range: GenomeRange):
        num_bins = self.__get_bins_num()
        self.__check_chrom_name(genome_range)
        scores_per_bin = self.fetch_scores(genome_range, num_bins)
        return scores_per_bin

    def __get_bins_num(self, default_num=700):
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

    def __check_chrom_name(self, genome_range):
        if genome_range.chrom not in self.bw.chroms().keys():
            genome_range.change_chrom_names()

        if genome_range.chrom not in self.bw.chroms().keys():
            log.warning("Can not read region {} from bigwig file:\n\n"
                        "{}\n\nPlease check that the chromosome name is part of the bigwig file "
                        "and that the region is valid".format(str(genome_range), self.properties['file']))


class ABCompartment(BigWig):
    """
    A/B Comapartment BigWig track.
    """

    DEFAULT_POSITIVE_COLOR = "#ff9c9c"
    DEFAULT_NEGATIVE_COLOR = "#66ccff"

    def __init__(self, file_, **kwargs):
        properties_dict = {
            'positive_color': ABCompartment.DEFAULT_POSITIVE_COLOR,
            'negative_color': ABCompartment.DEFAULT_NEGATIVE_COLOR,
        }
        properties_dict.update(kwargs)
        super().__init__(file_, **properties_dict)
