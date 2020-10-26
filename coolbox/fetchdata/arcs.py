import pandas as pd
from intervaltree import IntervalTree, Interval

from coolbox.fetchdata.base import FetchTrackData
from coolbox.utilities import split_genome_range, change_chrom_names

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class FetchArcs(FetchTrackData):

    def __init__(self, *args, **kwargs):
        valid_intervals = 0
        interval_tree = {}
        line_number = 0
        with open(self.properties['file'], 'r') as file_h:
            for line in file_h.readlines():
                line_number += 1
                if line.startswith('browser') or line.startswith('track') or line.startswith('#') or self.__is_header(line):
                    continue
                try:
                    chrom1, start1, end1, chrom2, start2, end2, score, *other = line.strip().split('\t')
                except Exception as detail:
                    msg = 'File not valid. The format is chrom1 start1, end1, ' \
                          'chrom2, start2, end2, score\nError: {}\n in line\n {}'.format(detail, line)
                    raise IOError(msg)

                try:
                    start1 = int(start1)
                    end1 = int(end1)
                    start2 = int(start2)
                    end2 = int(end2)
                except ValueError as detail:
                    msg = "Error reading line: {}. One of the fields is not " \
                          "an integer.\nError message: {}".format(line_number, detail)
                    raise IOError(msg)

                assert start1 <= end1, "Error in line #{}, end1 larger than start1 in {}".format(line_number, line)
                assert start2 <= end2, "Error in line #{}, end2 larger than start2 in {}".format(line_number, line)
                try:
                    score = float(score)
                except ValueError as detail:
                    msg = "Error reading line: {}. The score is not valid {}. " \
                          "\nError message: {}".format(line_number, score, detail)
                    raise IOError(msg)

                if chrom1 != chrom2:
                    log.warning("Only links in same chromosome are used. Skipping line\n{}\n".format(line))
                    continue

                if chrom1 not in interval_tree:
                    interval_tree[chrom1] = IntervalTree()

                if start2 < start1:
                    start1, start2 = start2, start1
                    end1, end2 = end2, end1

                # each interval spans from the smallest start to the largest end
                interval_tree[chrom1].add(Interval(start1, end2, score))
                valid_intervals += 1

        if valid_intervals == 0:
            log.warning("No valid intervals were found in file {}".format(self.properties['file']))

        file_h.close()
        self.interval_tree = interval_tree

    def fetch_data(self, genome_range):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        Return
        ------
        intervals : pandas.core.frame.DataFrame
            Arcs interval table.
        """
        return self.fetch_intervals(genome_range)

    def fetch_intervals(self, genome_range):
        """
        Fetch intervals within input chromosome range.
        """
        chrom, start, end = split_genome_range(genome_range)
        if chrom not in self.interval_tree:
            chrom = change_chrom_names(chrom)
            if chrom not in self.interval_tree:
                log.warning(f"chromosome of {str(genome_range)} not in Arc dataset.")
                return pd.DataFrame([], columns=['chromsome', 'start', 'end', 'score'])
        intervals = self.interval_tree[chrom][start:end]
        intervals = [(chrom, itv.begin, itv.end, float(itv.data)) for itv in intervals]
        intval_table = pd.DataFrame(intervals, columns=['chromsome', 'start', 'end', 'score'])
        return intval_table

    def __is_header(self, line):
        fields = line.split()
        for idx in [1, 2, 4, 5]:
            if not fields[idx].isdigit():
                return True
        return False