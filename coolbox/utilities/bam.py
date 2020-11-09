import subprocess as subp
import os.path as osp

import numpy as np


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
