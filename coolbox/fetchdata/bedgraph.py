import pandas as pd

from coolbox.fetchdata.base import FetchTrackData
from coolbox.utilities import split_genome_range, change_chrom_names


class FetchBedGraph(FetchTrackData):

    def fetch_data(self, genome_range):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        Return
        ------
        intervals : pandas.core.frame.DataFrame
            Bed graph interval table.
        """
        return self.fetch_intervals(genome_range)

    def fetch_intervals(self, genome_range):
        """
        Fetch intervals within input chromosome range.
        """
        chrom, start, end = split_genome_range(genome_range)
        if chrom not in self.interval_tree:
            chrom_ = change_chrom_names(chrom)
        else:
            chrom_ = chrom

        intervals = sorted(self.interval_tree[chrom_][start:end])

        col_chrom = [chrom] * len(intervals)
        col_start = []
        col_end   = []
        col_score = []
        for itv in intervals:
            col_start.append(itv.begin)
            col_end.append(itv.end)
            col_score.append(itv.data[0])

        intval_table = pd.DataFrame(
            {
                'chromsome': col_chrom,
                'start': col_start,
                'end':   col_end,
                'score': col_score,
            },
            columns=['chromsome', 'start', 'end', 'score']
        )

        return intval_table