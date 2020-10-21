import os.path as osp
import pandas as pd

from coolbox.fetchdata.base import FetchTrackData
from coolbox.utilities.bed import (
    bgz_bed, index_bed, query_bed
)
from coolbox.utilities import (
    change_chrom_names,
    Interval, IntervalTree,
    ReadBed,
    get_logger, split_genome_range,
    to_gr,
)


log = get_logger(__name__)


def build_bed_index(file):
    if file.endswith(".bgz"):
        bgz_file = file
    else:
        bgz_file = file + '.bgz'
        log.info(f"Bgzip bed file, save to {bgz_file}")
        bgz_bed(file, bgz_file)
    if not osp.exists(bgz_file + '.tbi'):
        log.info(f"Make tabix of bgz file, save to {bgz_file}.tbi")
        index_bed(bgz_file)
    return bgz_file


class FetchBed(FetchTrackData):

    def __init__(self, *args, **kwargs):
        file = self.properties['file']
        self.bgz_file = build_bed_index(file)

    def load_range(self, genome_range):
        genome_range = to_gr(genome_range)
        valid_intervals, min_score, max_score = self.__load(genome_range)
        if valid_intervals == 0:
            genome_range.change_chrom_names()
            valid_intervals, min_score, max_score = self.__load(genome_range)
            if valid_intervals == 0:
                log.debug("No valid intervals were found in file {} within range {}".format(
                    self.properties['file'],
                    f"{str(genome_range)}"
                ))
        self.score_range = (min_score, max_score)

    def __load(self, genome_range):
        valid_intervals = 0
        interval_tree = self.interval_tree
        max_score = float('-inf')
        min_score = float('inf')

        chrom, start, end = split_genome_range(genome_range)
        try:
            bed_file_h = ReadBed(query_bed(self.bgz_file, chrom, start, end))
        except StopIteration:
            log.info(f"No records in the range {str(genome_range)}")
            return valid_intervals, min_score, max_score
        self.bed_type = bed_file_h.file_type

        for bed in bed_file_h:
            if bed.score < min_score:
                min_score = bed.score
            if bed.score > max_score:
                max_score = bed.score

            if bed.chromosome not in interval_tree:
                interval_tree[bed.chromosome] = IntervalTree()

            itv = Interval(bed.start, bed.end, bed)
            if itv not in interval_tree:
                interval_tree[bed.chromosome].add(itv)
            valid_intervals += 1
        return valid_intervals, min_score, max_score

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
        self.load_range(genome_range)
        chrom, start, end = split_genome_range(genome_range)
        if chrom not in self.interval_tree:
            chrom = change_chrom_names(chrom)
        if chrom not in self.interval_tree:
            intervals = []
        else:
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