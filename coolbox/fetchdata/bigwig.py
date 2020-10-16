import pandas as pd

from coolbox.fetchdata.base import FetchTrackData
from coolbox.utilities import split_genome_range, change_chrom_names


class FetchBigWig(FetchTrackData):

    DEFAULT_NUM_BINS = 700

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
        try:
            num_bins = int(self.properties['number_of_bins'])
        except KeyError:
            num_bins = FetchBigWig.DEFAULT_NUM_BINS
        scores = self.fetch_intervals(genome_range)
        return scores

    def fetch_scores(self, genome_range, num_bins=700):
        """
        Fetch bins scores within input chromosome range.
        """
        chrom, start, end = split_genome_range(genome_range)
        if chrom not in self.bw.chroms():
            chrom = change_chrom_names(chrom)
        scores = self.bw.stats(chrom, start, end, nBins=num_bins)
        return scores

    def fetch_intervals(self, genome_range):
        """
        Fetch BigWig intervals within input chromosome range.
        """
        chrom, start, end = split_genome_range(genome_range)
        if chrom not in self.bw.chroms():
            chrom_ = change_chrom_names(chrom)
        else:
            chrom_ = chrom

        intervals = self.bw.intervals(chrom_, start, end)

        col_chrom = [chrom] * len(intervals)
        col_start = []
        col_end   = []
        col_score = []
        for s, e, v in intervals:
            col_start.append(s)
            col_end.append(e)
            col_score.append(v)

        intval_table = pd.DataFrame(
            {
                "chromsome": col_chrom,
                "start": col_start,
                "end":   col_end,
                "score": col_score,
            },
            columns=['chromsome', 'start', 'end', 'score'])

        return intval_table