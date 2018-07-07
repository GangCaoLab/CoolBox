from collections import OrderedDict

import pandas as pd

from .utilities import GenomeRange, change_chrom_names

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = [
    "FetchBigWig", "FetchBedGraph", "FetchBed",
    "FetchArcs", "FetchCool", "FetchFrame"
]



def split_genome_range(genome_range):
    """
    Little helper func.
    enforce genome_range is a GenomeRange object, and
    split genome_range to chrom, start, end
    """
    if isinstance(genome_range, str):
        genome_range = GenomeRange(genome_range)
    else:
        assert isinstance(genome_range, GenomeRange), \
            "genome_range is a `GenomeRange` object or a genome range str"
    chrom, start, end = genome_range.chrom, genome_range.start, genome_range.end
    return chrom, start, end


class FetchFrame(object):
    """
    FetchFrame used for fetch data from all tracks in the frame.
    """

    def fetch_data(self, genome_range=None):
        if genome_range is None:
            genome_range = self.current_range

        tracks_data = OrderedDict()
        for name, track in self.tracks.items():
            if hasattr(track, 'fetch_data'):
                data = track.fetch_data(genome_range)
            else:
                data = []
            tracks_data.update([(name, data)])

        return tracks_data


class FetchTrackData(object):
    """
    FetchTrackData base class.
    FetchTrackData object used for fetch data from specify data format.

    All `FetchTrackData` must have a `fetch_data` method.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


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


class FetchArcs(FetchTrackData):

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
        intervals = self.interval_tree[chrom][start:end]
        intervals = [(chrom, itv.begin, itv.end, float(itv.data)) for itv in intervals]
        intval_table = pd.DataFrame(intervals, columns=['chromsome', 'start', 'end', 'score'])
        return intval_table


class FetchCool(FetchTrackData):

    def fetch_data(self, genome_range1, genome_range2=None):
        """
        Parameters
        ----------
        genome_range1 : {str, GenomeRange}

        genome_range2 : {str, GenomeRange}, optional.

        balance : {bool, optional}
            Balance matrix or not, default False.

        Return
        ------
        pixels : pandas.core.frame.DataFrame
            Hi-C pixels table.
            The pixel table contains the non-zero upper triangle entries of the contact map.
        """
        return self.fetch_pixels(genome_range1, genome_range2)

    def fetch_array(self, genome_range, balance=False):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        balance : bool, optional
            balance matrix or not, default False.

        Return
        ------
        arr : numpy.ndarray
            Intervals within input chromosome range.
        """
        chrom, start, end = split_genome_range(genome_range)
        arr = self.cool.matrix(balance=balance).fetch(str(GenomeRange(chrom, start, end)))
        return arr

    def fetch_pixels(self, genome_range1,
                           genome_range2=None,
                           join=True):
        mat = self.cool.matrix(as_pixels=True, join=join)
        chrom1, start1, end1 = split_genome_range(genome_range1)

        if chrom1 not in self.cool.chromnames:
            chrom1 = change_chrom_names(chrom1)

        if genome_range2 is not None:
            chrom2, start2, end2 = split_genome_range(genome_range2)
            if chrom2 not in self.cool.chromnames:
                chrom2 = change_chrom_names(chrom2)
            pixels = mat.fetch(str(GenomeRange(chrom1, start1, end1)),
                               str(GenomeRange(chrom2, start2, end2)))
        else:
            pixels = mat.fetch(str(GenomeRange(str(GenomeRange(chrom1, start1, end1)))))

        return pixels
