import os.path as osp
import subprocess as subp

import pandas as pd

from coolbox.fetchdata.base import FetchTrackData
from coolbox.utilities import (
    split_genome_range, change_chrom_names,
    GenomeRange, get_logger
)
from coolbox.utilities.bed import bgz_bed, tabix_query


log = get_logger(__name__)


def index_bedgraph(bgz_file):
    cmd = ['tabix', '-b', '2', '-e', '3', bgz_file]
    subp.check_call(cmd)


def build_bedgraph_bgz(file):
    if file.endswith(".bgz"):
        bgz_file = file
    else:
        bgz_file = file + '.bgz'
        log.info(f"Bgzip bedgraph file, save to {bgz_file}")
        bgz_bed(file, bgz_file)
    if not osp.exists(bgz_file + '.tbi'):
        log.info(f"Make tabix of bgz file, save to {bgz_file}.tbi")
        index_bedgraph(bgz_file)
    return bgz_file


class FetchBedGraph(FetchTrackData):

    def __init__(self):
        file = self.properties['file']
        self.bgz_file = build_bedgraph_bgz(file)

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
        gr = GenomeRange(chrom, start, end)

        rows = self.__load(gr)
        if len(rows) == 0:
            chrom = change_chrom_names(chrom)
            rows = self.__load(GenomeRange(chrom, start, end))

        intval_table = pd.DataFrame(
            rows,
            columns=['chromsome', 'start', 'end', 'score']
        )

        return intval_table

    def __load(self, genome_range):
        rows = []
        gr = genome_range
        for it in tabix_query(self.bgz_file, gr.chrom, gr.start, gr.end, split=True):
            rows.append([
                it[0], int(it[1]), int(it[2]), float(it[3])
            ])
        return rows
