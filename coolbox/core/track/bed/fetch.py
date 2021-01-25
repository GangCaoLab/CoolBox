from typing import Tuple, Sequence, Union
import pandas as pd

from coolbox.utilities import (
    get_logger
)
from coolbox.utilities.bed import (
    ReadBed, query_bed
)

from coolbox.utilities.genome import GenomeRange

log = get_logger(__name__)


class FetchBed(object):
    def fetch_intervals(self, bgz_file, gr: GenomeRange) -> pd.DataFrame:
        """
        Fetch intervals within input chromosome range.
        """
        intervals, bed_type = self.load_range(bgz_file, gr)
        if len(intervals) == 0:
            gr.change_chrom_names()
            intervals, bed_type = self.load_range(bgz_file, gr)
            if len(intervals) == 0:
                log.debug(f"No valid intervals were found in file {bgz_file} within range {gr}")
        intval_table = self.intervals2dataframe(intervals, bed_type)
        # attach bed_type info onto the dataframe
        intval_table.bed_type = bed_type
        return intval_table

    @staticmethod
    def load_range(bgz_file, gr: GenomeRange) -> Tuple[Sequence, Union[str, None]]:
        intervals = []
        try:
            bed_iterator = ReadBed(query_bed(bgz_file, gr.chrom, gr.start, gr.end))
        except StopIteration:
            log.info(f"No records in the range {str(gr)}")
            return [], None

        for bed in bed_iterator:
            intervals.append(bed)

        return intervals, bed_iterator.file_type

    @staticmethod
    def intervals2dataframe(intervals, bed_type: Union[str, None]) -> pd.DataFrame:
        """
        Convert intervals list to pandas.DataFrame
        """
        bed_fields = ['chromosome', 'start', 'end',
                      'name', 'score', 'strand',
                      'thick_start', 'thick_end',
                      'rgb', 'block_count',
                      'block_sizes', 'block_starts']

        if bed_type == 'bed12':
            fields = bed_fields
        elif bed_type == 'bed9':
            fields = bed_fields[:9]
        else:
            fields = bed_fields[:6]

        df = pd.DataFrame(intervals, columns=fields)

        return df
