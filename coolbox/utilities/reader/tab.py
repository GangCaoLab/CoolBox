import abc
import platform
import typing as T
from pathlib import Path
import subprocess as subp
from functools import partial
from os import path as osp

import numpy as np
import oxbow as ox
import pandas as pd

from ..logtools import get_logger
from ..genome import GenomeRange
from ..cmd import check_tool, ensure_unix, ensure_tool_installed

log = get_logger(__name__)


def bgz_bed(bed_path, bgz_path):
    ensure_unix()
    cmd = ""
    cmd += "zcat" if bed_path.endswith(".gz") else "cat"
    subp.check_call(cmd + f" {bed_path} | sort -k1,1 -k2,2n | bgzip > {bgz_path}",
                    shell=True)
    return bgz_path


def index_bed(bgz_path):
    ensure_unix()
    cmd = ['tabix', '-p', 'bed', bgz_path]
    subp.check_call(cmd)


def tabix_query(filename, chrom, start, end, split=True):
    """Call tabix and generate an array of strings for each line it returns."""
    ensure_unix()
    query = '{}:{}-{}'.format(chrom, start, end)
    p = subp.Popen(['tabix', '-f', filename, query], stdout=subp.PIPE)
    for line in p.stdout:
        line = line.decode('utf-8')
        if not split:
            yield line
        else:
            yield line.strip().split('\t')


query_bed = partial(tabix_query, split=False)


def index_bedgraph(bgz_file):
    cmd = ['tabix', '-b', '2', '-e', '3', bgz_file]
    subp.check_call(cmd)


def build_bedgraph_bgz(file):
    file = osp.expanduser(file)
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


def build_bed_index(file):
    file = osp.expanduser(file)
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


def build_snp_index(file, col_chrom, col_pos):
    file = osp.expanduser(file)
    c = col_chrom + 1
    p = col_pos + 1
    if file.endswith(".bgz"):
        bgz_file = file
    elif osp.exists(file + '.bgz'):
        bgz_file = file + '.bgz'
    else:
        bgz_file = file + '.bgz'
        cmd = "zcat" if file.endswith('.gz') else "cat"
        cmd += f" {file} | sort -k{c},{c} -k{p},{p}n | bgzip > {bgz_file}"
        subp.check_call(cmd, shell=True)
    index_file = bgz_file + '.tbi'
    if not osp.exists(index_file):
        cmd = ['tabix', '-s', str(c), '-b', str(p), '-e', str(p), bgz_file]
        subp.check_call(cmd)
    return bgz_file


def bgz_bedpe(bedpe_path, bgz_path):
    if not osp.exists(bgz_path):
        cmd = f"sort -k1,1 -k4,4 -k2,2n -k5,5n {bedpe_path} | bgzip > {bgz_path}"
        subp.check_call(cmd, shell=True)


def index_bedpe(bgz_path):
    cmd = f"pairix -f -s 1 -d 4 -b 2 -e 3 -u 5 -v 6 {bgz_path}".split(" ")
    subp.check_call(cmd)


def pairix_query(bgz_file, query: GenomeRange, second: T.Optional[GenomeRange] = None,
                 open_region: bool = False, split: bool = True):
    if second:
        query = f"{query}|{second}"
    else:
        if open_region:
            query = f"{query}|{query.chrom}"
    cmd = ['pairix', str(bgz_file), str(query)]
    p = subp.Popen(cmd, stdout=subp.PIPE)
    for line in p.stdout:
        line = line.decode('utf-8')
        if not split:
            yield line
        else:
            yield line.strip().split('\t')


def process_bedpe(path):
    if path.endswith(".bgz"):
        bgz_file = path
    else:
        bgz_file = path + ".bgz"
        bgz_bedpe(path, bgz_file)
    if not osp.exists(f"{bgz_file}.px2"):
        index_bedpe(bgz_file)
    return bgz_file


def bgz_pairs(pairs_path, bgz_path):
    if not osp.exists(bgz_path):
        cmd = f"grep -v '#' {pairs_path} | sort -k2,2 -k4,4 -k3,3n -k5,5n | bgzip > {bgz_path}"
        subp.check_call(cmd, shell=True)


def index_pairs(bgz_path):
    cmd = f"pairix -f -p pairs {bgz_path}".split(" ")
    subp.check_call(cmd)


def process_pairs(path):
    if path.endswith(".bgz"):
        bgz_file = path
    else:
        bgz_file = path + ".bgz"
        bgz_pairs(path, bgz_file)
    if not osp.exists(f"{bgz_file}.px2"):
        index_pairs(bgz_file)
    return bgz_file


def _process_gtf(gtf_path, out_path):
    cmd = f'(grep ^"#" {gtf_path}; grep -v ^"#" {gtf_path} | sort -k1,1 -k4,4n) | bgzip > {out_path}'
    subp.check_call(cmd, shell=True)


def _gtf_gz_to_bgz(gz, bgz):
    cmd = f'gunzip -c {gz} | grep -v ^"#" | bgzip > {bgz}'
    subp.check_call(cmd, shell=True)


def tabix_index(filename, preset="gff"):
    """Call tabix to create an index for a bgzip-compressed file."""
    subp.check_call([
        'tabix', '-p', preset, filename
    ])


def build_gtf_index(file):
    file = osp.expanduser(file)
    if file.endswith(".gtf"):
        bgz_file = file + ".bgz"
        if not osp.exists(bgz_file):
            log.info(f"Process the gtf and do bgzip, save to {bgz_file}.")
            _process_gtf(file, bgz_file)
    elif file.endswith(".gtf.gz"):
        bgz_file = file.rstrip(".gz") + ".bgz"
        log.info(f"Convert .gtf.gz to .gtf.bgz, save to {bgz_file}.")
        if not osp.exists(bgz_file):
            _gtf_gz_to_bgz(file, bgz_file)
    elif file.endswith(".gtf.bgz"):
        bgz_file = file
    else:
        raise IOError(f"GTF track only support GTF file(.gtf or .gtf.gz), got {file}.")

    idx_file = bgz_file + ".tbi"
    if not osp.exists(idx_file):
        log.info(f"Tabix index not found, build it in {idx_file}")
        tabix_index(bgz_file)
    return bgz_file


def _is_bam_sorted(bam_path):
    p = subp.Popen(['samtools', 'view', '-H', bam_path], stdout=subp.PIPE)
    for line in p.stdout:
        line = line.decode("utf-8")
        if "SO:unsorted" in line:
            return False
    return True


def process_bam(bam_path):
    """Sort and index a BAM file.
    If input is a SAM file, it will be converted to a BAM file first."""
    if bam_path.endswith(".bam"):
        bai_path = bam_path + '.bai'
        if osp.exists(bai_path):
            return bam_path
        if osp.exists(bam_path[:-4] + '.sorted.bam.bai'):
            return bam_path[:-4] + '.sorted.bam.bai'

        ensure_unix()
        ensure_tool_installed("samtools")
        if not _is_bam_sorted(bam_path):
            sorted_bam_path = bam_path[:-4] + '.sorted.bam'
            subp.check_call(['samtools', 'sort', bam_path, '-o', sorted_bam_path])
        else:
            sorted_bam_path = bam_path
        subp.check_call(['samtools', 'index', sorted_bam_path])
    elif bam_path.endswith(".sam"):
        ensure_unix()
        ensure_tool_installed("samtools")
        sorted_bam_path = bam_path[:-4] + '.sorted.bam'
        subp.check_call(['samtools', 'sort', bam_path, '-o', sorted_bam_path])
        subp.check_call(['samtools', 'index', sorted_bam_path])
    else:
        raise IOError("BAM input file should be in .bam or .sam format")
    return sorted_bam_path


def query_bam(filename, chrom, start, end, split=True):
    """Call tabix and generate an array of strings for each line it returns."""
    ensure_unix()
    ensure_tool_installed("samtools")
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


def _parse_samtools_cov(lines):
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


def coverage_by_samtools(bam_path, region, bins):
    ensure_unix()
    ensure_tool_installed("samtools")
    cmd = ["samtools", "coverage", bam_path, "-r", region, "-w", str(bins)]
    p = subp.Popen(cmd, stdout=subp.PIPE)
    lines = []
    for line in p.stdout:
        line = line.decode('utf-8')
        lines.append(line)
    covs = _parse_samtools_cov(lines)
    return np.array(covs)


class TabFileReader(abc.ABC):
    """Generic tab-separated file reader.

    Including:
    - BED
    - BedGraph
    - bigWig
    - bigBED
    - GTF
    - BAM
    """
    def __init__(self, path: Path, **params):
        self.path = path
        self.params = params

    @abc.abstractmethod
    def query(self, chrom: str, start: int, end: int) -> pd.DataFrame:
        """Query the file"""
        pass


class TabFileReaderWithOxbow(TabFileReader):
    def __init__(self, path: Path):
        super().__init__(path)
        self.bgz_path = Path(str(path) + ".bgz")
        suffix = path.suffix
        if suffix == ".gtf":
            ds = ox.from_gtf(self.bgz_path)
        elif suffix in [".bed", ".bedgraph"]:
            ds = ox.from_bed(self.bgz_path)
        elif suffix in ['.bw', '.bigwig']:
            ds = ox.from_bigwig(self.path)
        elif suffix == 'bam':
            ds = ox.from_bam(self.path)
        else:
            raise ValueError(f"Unsupported file type: {suffix}")
        self.ds = ds

    def query(self, chrom: str, start: int, end: int) -> pd.DataFrame:
        sub = self.ds.regions(f"{chrom}:{start}-{end}")
        return sub.pd()


def index_tab_file(path: Path):
    path_str = str(path)


def get_indexed_tab_reader(path: str) -> TabFileReader:
    if path.endswith(".bam"):
        sorted_bam_path = process_bam(path)
        return TabFileReaderWithOxbow(sorted_bam_path)
    elif path.endswith(".bw") or path.lower().endswith(".bigwig"):
        return TabFileReaderWithOxbow(path)


