import os.path as osp
import subprocess as subp
import pandas as pd
import numpy as np

from coolbox.fetchdata.base import FetchTrackData
from coolbox.utilities import split_genome_range

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def is_bam_sorted(bam_path):
    p = subp.Popen(['samtools', 'view', '-H', bam_path], stdout=subp.PIPE)
    for line in p.stdout:
        line = line.decode("utf-8")
        if "SO:unsorted" in line:
            return False
    return True


def process_bam(bam_path):
    if bam_path.endswith(".bam"):
        bai_path = bam_path + '.bai'
        if osp.exists(bai_path):
            return bam_path
        if not is_bam_sorted(bam_path):
            sorted_bam_path = bam_path[:-4] + '.sorted.bam'
            subp.check_call(['samtools', 'sort', bam_path, '-o', sorted_bam_path])
        else:
            sorted_bam_path = bam_path
        subp.check_call(['samtools', 'index', sorted_bam_path])
    elif bam_path.endswith(".sam"):
        sorted_bam_path = bam_path[:-4] + '.sorted.bam'
        subp.check_call(['samtools', 'sort', bam_path, '-o', sorted_bam_path])
        subp.check_call(['samtools', 'index', sorted_bam_path])
    else:
        raise IOError("BAM input file should be in .bam or .sam format")
    return sorted_bam_path


def query_bam(filename, chrom, start, end, split=True):
    """Call tabix and generate an array of strings for each line it returns."""
    query = '{}:{}-{}'.format(chrom, start, end)
    p = subp.Popen(['samtools', 'view', filename, query], stdout=subp.PIPE)
    for line in p.stdout:
        line = line.decode('utf-8')
        if not split:
            yield line
        else:
            items = line.strip().split('\t')
            items_ = items[:11] + ["\t".join(items[12:])]
            yield items_


def coverage_by_samtools(bam_path, region, bins):
    cmd = ["samtools", "coverage", bam_path, "-r", region, "-w", str(bins)]
    p = subp.Popen(cmd, stdout=subp.PIPE)
    lines = []
    for line in p.stdout:
        line = line.decode('utf-8')
        lines.append(line)
    covs = parse_samtools_cov(lines)
    return np.array(covs)


def parse_samtools_cov(lines):
    covs = {}
    for line in lines[1:-1]:
        left, mid, _ = line.split("│")
        percent = float(left.strip("> %"))
        for i, c in enumerate(mid):
            covs.setdefault(i, 0)
            if c != ' ' and covs[i] == 0:
                covs[i] = percent
    covs = [covs[i] for i in sorted(covs.keys())]
    return covs


class FetchBAM(FetchTrackData):

    def __init__(self, *args, **kwargs):
        bam_path = self.properties['file']
        self.indexed_bam = process_bam(bam_path)

    def fetch_data(self, genome_range):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        Return
        ------
        intervals : pandas.core.frame.DataFrame
            Sam interval table.
        """
        return self.fetch_intervals(genome_range)

    def fetch_intervals(self, genome_range):
        chrom, start, end = split_genome_range(genome_range)
        rows = []
        for row_items in query_bam(self.indexed_bam, chrom, start, end, split=True):
            rows.append(row_items)
        # https://samtools.github.io/hts-specs/SAMv1.pdf
        fields = ["qname", "flag", "rname", "pos", "mapq", "cigar",
                  "rnext", "pnext", "tlen", "seq", "qual", "options"]
        df = pd.DataFrame(rows, columns=fields)
        if df.shape[0] > 0:
            df['flag'] = df['flag'].astype(int)
            df['pos'] = df['pos'].astype(int)
            df['mapq'] = df['mapq'].astype(int)
        return df

    def fetch_coverage(self, genome_range, bins=100):
        scores_per_bin = coverage_by_samtools(
            self.indexed_bam,
            str(genome_range),
            bins
        )
        return scores_per_bin
