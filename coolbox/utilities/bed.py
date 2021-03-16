import typing as T
import subprocess as subp
from functools import partial
from os import path as osp

from intervaltree import Interval, IntervalTree
import collections
import types

from .filetool import opener, to_string
from .logtools import get_logger
from .genome import GenomeRange

log = get_logger(__name__)


def file_to_intervaltree(file_name):
    """
    converts a BED like file into a bx python interval tree

    Parameters
    ----------
    file_name : str
        Path to file.

    Return
    ------
    interval_tree : dict
        Interval tree dictionary. They key is the chromosome/contig name and the
        value is an IntervalTree. Each of the intervals have as 'value' the fields[3:] if any.
    """
    import sys
    import numpy as np

    # iterate over a BED like file
    # saving the data into an interval tree
    # for quick retrieval
    file_h = opener(file_name)
    line_number = 0
    valid_intervals = 0
    prev_chrom = None
    prev_start = -1
    prev_line = None
    interval_tree = {}
    min_value = np.inf
    max_value = -np.inf

    for line in file_h.readlines():
        line_number += 1
        line = to_string(line)
        if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        try:
            chrom, start, end = fields[0:3]
        except Exception as detail:
            msg = "Error reading line: {}\nError message: {}".format(line_number, detail)
            sys.exit(msg)

        try:
            start = int(start)
        except ValueError as detail:
            msg = "Error reading line: {}. The start field is not " \
                  "an integer.\nError message: {}".format(line_number, detail)
            sys.exit(msg)

        try:
            end = int(end)
        except ValueError as detail:
            msg = "Error reading line: {}. The end field is not " \
                  "an integer.\nError message: {}".format(line_number, detail)
            sys.exit(msg)

        if prev_chrom == chrom:
            assert prev_start <= start, \
                "BED file not sorted. Please use a sorted bed file.\n{}{} ".format(prev_line, line)

        if chrom not in interval_tree:
            interval_tree[chrom] = IntervalTree()

        value = None

        if len(fields) > 3:
            value = fields[3:]
            try:
                line_min = min(map(float, value))
                if line_min < min_value:
                    min_value = line_min

                line_max = max(map(float, value))
                if line_max > max_value:
                    max_value = line_max
            except ValueError:
                pass

        assert end > start, "Start position larger or equal than end for line\n{} ".format(line)

        interval_tree[chrom].add(Interval(start, end, value))
        valid_intervals += 1

    if valid_intervals == 0:
        log.warning("No valid intervals were found in file {}".format(file_name))
    file_h.close()

    return interval_tree, min_value, max_value


class ReadBed(object):
    """
    Reads a bed file. Based on the number of fields
    it tries to guess the type of bed file used. Current options
    are bed3, bed6 and bed12

    Example:
    bed = readBed(open("file.bed", 'r'))
    for interval in bed:
        print(interval['start'])

    """

    def __init__(self, file_handle):

        self.file_type = None
        self.file_handle = file_handle
        self.line_number = 0
        # guess file type
        fields = self.get_no_comment_line()
        fields = to_string(fields)
        fields = fields.split()

        self.guess_file_type(fields)
        if type(file_handle) is types.GeneratorType:
            self._file_name = "<from generator>"
        else:
            self._file_name = file_handle.name
            self.file_handle.seek(0)
        self.prev_chrom = None
        self.prev_start = -1
        self.prev_line = None

        # list of bed fields
        self.fields = ['chromosome', 'start', 'end',
                       'name', 'score', 'strand',
                       'thick_start', 'thick_end',
                       'rgb', 'block_count',
                       'block_sizes', 'block_starts']

        if self.file_type == 'bed12':
            self.BedInterval = collections.namedtuple('BedInterval', self.fields)
        elif self.file_type == 'bed9':
            self.BedInterval = collections.namedtuple('BedInterval', self.fields[:9])
        else:
            self.BedInterval = collections.namedtuple('BedInterval', self.fields[:6])

    def __iter__(self):
        return self

    def get_no_comment_line(self):
        """
        Skips comment lines starting with '#'
        "track" or "browser" in the bed files
        """
        line = next(self.file_handle)
        line = to_string(line)
        if line.startswith("#") or line.startswith("track") or \
                line.startswith("browser") or line.strip() == '':
            line = self.get_no_comment_line()

        self.line_number += 1
        return line

    def guess_file_type(self, line_values):
        """try to guess type of bed file by counting the fields
        """
        if len(line_values) == 3:
            self.file_type = 'bed3'
        elif len(line_values) == 4:
            self.file_type = 'bedgraph'
        elif len(line_values) == 6:
            self.file_type = 'bed6'
        elif len(line_values) == 12:
            self.file_type = 'bed12'
        elif len(line_values) == 9:
            # this is a case where a specific color is encoded in the 10 field of the bed file
            self.file_type = 'bed9'
        elif len(line_values) > 6:
            # assume bed6
            self.file_type = 'bed6'
            log.warning("Number of fields in BED file is not standard. Assuming bed6\n")
        else:
            # assume bed3
            self.file_type = 'bed3'
            log.warning("Number of fields in BED file is not standard. Assuming bed3\n")
        return self.file_type

    def next(self):
        """
        Return
        ------
        bedInterval object
        """
        line = self.get_no_comment_line()

        bed = self.get_bed_interval(line)
        if self.prev_chrom == bed.chromosome:
            assert self.prev_start <= bed.start, \
                "BED file not sorted. Please use a sorted bed file.\n" \
                "File: {}\n" \
                "Previous line: {}\n Current line{} ".format(self._file_name, self.prev_line, line)

        self.prev_chrom = bed.chromosome
        self.prev_start = bed.start
        self.prev_line = line

        return bed

    def __next__(self):
        """
        Return
        ------
        bedInterval object
        """
        line = self.get_no_comment_line()

        bed = self.get_bed_interval(line)
        if self.prev_chrom == bed.chromosome:
            assert self.prev_start <= bed.start, \
                "BED file not sorted. Please use a sorted bed file.\n" \
                "File: {}\n" \
                "Previous line: {}\n Current line{} ".format(self._file_name, self.prev_line, line)

        self.prev_chrom = bed.chromosome
        self.prev_start = bed.start
        self.prev_line = line

        return bed

    def get_bed_interval(self, bed_line):
        r"""
        Processes each bed line from a bed file, casts the values and returns
        a namedtuple object

        >>> bed_line="chr1\t0\t1000\tgene_1\t0.5\t-\t0\t1000\t0\t3\t10,20,100\t20,200,700"
        >>> with open('/tmp/test.bed', 'w') as fh:
        ...     foo = fh.write(bed_line)
        >>> bed_f = ReadBed(open('/tmp/test.bed','r'))
        >>> bed = bed_f.get_bed_interval(bed_line)
        >>> bed.chromosome
        'chr1'
        >>> bed.block_starts
        [20, 200, 700]

        >>> bed_line="chr2\t0\t1000\tgene_1\t0.5\t-\n"
        >>> with open('/tmp/test.bed', 'w') as fh:
        ...     foo = fh.write(bed_line)
        >>> bed_f = ReadBed(open('/tmp/test.bed','r'))
        >>> bed_f.get_bed_interval(bed_line)
        BedInterval(chromosome='chr2', start=0, end=1000, name='gene_1', score=0.5, strand='-')
        """

        line_data = bed_line.strip()
        line_data = to_string(line_data)
        line_data = line_data.split("\t")

        if self.file_type == 'bed12':
            assert len(line_data) == 12, "File type detected is bed12 but line {}: {} does " \
                                         "not have 12 fields.".format(self.line_number, bed_line)

        elif self.file_type == 'bed3':
            assert len(line_data) == 3, "File type detected is bed3 but line {}: {} does " \
                                        "not have 3 fields.".format(self.line_number, bed_line)

        elif self.file_type == 'bed6':
            assert len(line_data) == 6, "File type detected is bed6 but line {}: {} does " \
                                        "not have 6 fields.".format(self.line_number, bed_line)
        line_values = []
        for idx, r in enumerate(line_data):
            # first field is always chromosome/contig name
            # and should be cast as a string
            # same for field 3 (name)
            if idx in [0, 3]:
                line_values.append(r)
            # check field strand
            elif idx == 5:
                if r not in ['+', '-', '.']:
                    if r == '1':
                        r = '+'
                    elif r == '-1':
                        r = '-'
                    else:
                        log.warning("*Warning, invalid strand value found {} for line #{}:\n{}\n "
                                    "Setting strand to '.'\n".format(r, bed_line, self.line_number))
                        r = '.'
                line_values.append(r)

            elif idx in [1, 2, 6, 7, 9]:
                # start and end fields must be integers, same for thichStart(6),
                # and thickEnd(7) and blockCount(9) fields
                try:
                    line_values.append(int(r))
                except ValueError:
                    log.warning("Value: {} in field {} at line {} is not an integer\n".format(r, idx + 1,
                                                                                              self.line_number))
                    return dict()
            # check item rgb
            elif idx == 8:
                r = to_string(r)
                rgb = r.split(",")
                if len(rgb) == 3:
                    try:
                        r = list(map(int, rgb))
                    except ValueError as detail:
                        log.warning("Error reading line: #{}. The rgb field {} is not "
                                    "valid.\nError message: {}\n".format(self.line_number, r, detail))
                line_values.append(r)

            elif idx in [10, 11]:
                # this are the block sizes and block start positions
                r = to_string(r)
                r_parts = r.split(',')
                try:
                    r = [int(x) for x in r_parts if x != '']
                except ValueError as detail:
                    log.warning("Error reading line #{}. The block field {} is not "
                                "valid.\nError message: {}\n".format(self.line_number, r, detail))
                line_values.append(r)

            else:
                try:
                    tmp = float(r)
                except ValueError:
                    tmp = r
                except TypeError:
                    tmp = r
                line_values.append(tmp)

        assert line_values[2] > line_values[1], \
            "Start position larger or equal than end for line #{}:\n{}\n".format(self.line_number,
                                                                                 bed_line)

        if self.file_type == 'bed3':
            line_values = line_values[0:3]
            # in case of a bed3, the id, score and strand
            # values are added as ".", 0, "." respectively
            line_values.extend([".", 0, "."])
        elif self.file_type == 'bed6':
            line_values = line_values[0:6]

        return self.BedInterval._make(line_values)


def bgz_bed(bed_path, bgz_path):
    cmd = ""
    if bed_path.endswith(".gz"):
        cmd += "zcat"
    else:
        cmd += "cat"
    subp.check_call(cmd + f" {bed_path} | sort -k1,1 -k2,2n | bgzip > {bgz_path}",
                    shell=True)
    return bgz_path


def index_bed(bgz_path):
    cmd = ['tabix', '-p', 'bed', bgz_path]
    subp.check_call(cmd)


def tabix_query(filename, chrom, start, end, split=True):
    """Call tabix and generate an array of strings for each line it returns."""
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
        if file.endswith('.gz'):
            cmd = "zcat"
        else:
            cmd = "cat"
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
    file = osp.expanduser(file)
    if file.endswith(".gtf"):
        bgz_file = file + ".bgz"
        if not osp.exists(bgz_file):
            log.info(f"Process the gtf and do bgzip, save to {bgz_file}.")
            process_gtf(file, bgz_file)
    elif file.endswith(".gtf.gz"):
        bgz_file = file.rstrip(".gz") + ".bgz"
        log.info(f"Convert .gtf.gz to .gtf.bgz, save to {bgz_file}.")
        if not osp.exists(bgz_file):
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

