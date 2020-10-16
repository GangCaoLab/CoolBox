import pandas as pd

from coolbox.fetchdata.base import FetchTrackData
from coolbox.utilities import split_genome_range, change_chrom_names


class FetchBed(FetchTrackData):

    def fetch_data(self, genome_range):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        Return
        ------
        intervals : pandas.core.frame.DataFrame
            Bed interval table.
        """
        return self.fetch_intervals(genome_range)

    def fetch_intervals(self, genome_range):
        """
        Fetch intervals within input chromosome range.
        """
        chrom, start, end = split_genome_range(genome_range)
        if chrom not in self.interval_tree:
            chrom = change_chrom_names(chrom)
        intervals = sorted(self.interval_tree[chrom][start:end])
        intval_table = self.intervals2dataframe(intervals)
        return intval_table

    def intervals2dataframe(self, intervals):
        """
        Convert intervals list to pandas.DataFrame
        """
        bed_fields = ['chromosome', 'start', 'end',
                      'name', 'score', 'strand',
                      'thick_start', 'thick_end',
                      'rgb', 'block_count',
                      'block_sizes', 'block_starts']

        if self.bed_type == 'bed12':
            fields = bed_fields
        elif self.bed_type == 'bed9':
            fields = bed_fields[:9]
        else:
            fields = bed_fields[:6]

        rows = []
        for intv in intervals:
            bed = intv.data
            rows.append(tuple(bed))
        df = pd.DataFrame(rows, columns=fields)

        return df