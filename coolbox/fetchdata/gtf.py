import os
import os.path as osp
import subprocess as subp

from coolbox.fetchdata.base import FetchTrackData
from coolbox.utilities import split_genome_range, change_chrom_names

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def process_gtf(gtf_path, out_path):
    cmd = f'(grep ^"#" {gtf_path}; grep -v ^"#" {gtf_path} | sort -k1,1 -k4,4n) | bgzip > {out_path}'
    subp.check_call(cmd, shell=True)


def tabix_index(filename, preset="gff"):
    """Call tabix to create an index for a bgzip-compressed file."""
    subp.check_call([
        'tabix', '-p', preset, filename
    ])


def build_gtf_index(file):
    if file.endswith(".gtf"):
        bgz_file = file + ".gz"
        log.info(f"Process the gtf and do bgzip, save to {bgz_file}.")
        process_gtf(file, bgz_file)
    elif file.endswith(".gtf.gz"):
        bgz_file = file
    else:
        raise IOError(f"GTF track only support GTF file(.gtf or .gtf.gz), got {file}.")

    idx_file = bgz_file + ".tbi"
    if osp.exists(idx_file):
        log.info(f"Tabix index not found, build it in {idx_file}")
        tabix_index(bgz_file)


class FetchGTF(FetchTrackData):

    def __init__(self, *args, **kwargs):
        file = self.properties['file']
        build_gtf_index(file)
