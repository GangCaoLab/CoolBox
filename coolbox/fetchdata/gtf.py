import os.path as osp
import subprocess as subp

import pandas as pd

from coolbox.fetchdata.base import FetchTrackData
from coolbox.utilities import split_genome_range, tabix_query

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def process_gtf(gtf_path, out_path):
    cmd = f'(grep ^"#" {gtf_path}; grep -v ^"#" {gtf_path} | sort -k1,1 -k4,4n) | bgzip > {out_path}'
    subp.check_call(cmd, shell=True)


def gtf_gz_to_bgz(gz, bgz):
    cmd = f'gunzip -c {gz} | grep -v ^"#" | bgzip > {bgz}'
    subp.check_call(cmd, shell=True)


def tabix_index(filename, preset="gff"):
    """Call tabix to create an index for a bgzip-compressed file."""
    subp.check_call([
        'tabix', '-p', preset, filename
    ])


def build_gtf_index(file):
    if file.endswith(".gtf"):
        bgz_file = file + ".bgz"
        if not osp.exists(bgz_file):
            log.info(f"Process the gtf and do bgzip, save to {bgz_file}.")
            process_gtf(file, bgz_file)
    elif file.endswith(".gtf.gz"):
        bgz_file = file[-7:] + ".bgz"
        log.info(f"Convert .gtf.gz to .gtf.bgz, save to {bgz_file}.")
        gtf_gz_to_bgz(file, bgz_file)
    elif file.endswith(".gtf.bgz"):
        bgz_file = file
    else:
        raise IOError(f"GTF track only support GTF file(.gtf or .gtf.gz), got {file}.")

    idx_file = bgz_file + ".tbi"
    if not osp.exists(idx_file):
        log.info(f"Tabix index not found, build it in {idx_file}")
        tabix_index(bgz_file)
    return bgz_file


class FetchGTF(FetchTrackData):

    def __init__(self, *args, **kwargs):
        file = self.properties['file']
        self.bgz_file = build_gtf_index(file)

    def fetch_data(self, genome_range):
        return self.fetch_intervals(genome_range)

    def fetch_intervals(self, genome_range):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        Return
        ------
        intervals : pandas.core.frame.DataFrame
            Annotation interval table.
        """
        chrom, start, end = split_genome_range(genome_range)
        rows = []
        for row in tabix_query(self.bgz_file, chrom, start, end):
            rows.append(row)
        columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        df = pd.DataFrame(rows, columns=columns)
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)
        df['gene_name'] = df['attribute'].str.extract(".*gene_name (.*?) ").iloc[:, 0].str.strip('\";')
        df['gene_name'][df['gene_name'].isna()] = ""
        return df


