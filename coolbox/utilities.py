import os
import sys
import collections
from os.path import abspath, dirname, join
import io

import numpy as np
from intervaltree import IntervalTree, Interval

import logging


class GenomeRange(object):
    """
    Express a range on the genome.

    Attributes
    ----------
    chrom : str
        chromosome

    start : int
        start position

    end : int
        end position

    """

    def __init__(self, *args):
        """
        >>> range1 = GenomeRange("chr1", 1000, 2000)
        >>> str(range1)
        'chr1:1000-2000'
        >>> range2 = GenomeRange("chr2:2000-4000")
        >>> (range2.chrom, range2.start, range2.end)
        ('chr2', 2000, 4000)
        >>> range3 = GenomeRange("chr1", 2000, 1000)
        Traceback (most recent call last):
        ...
        ValueError: Please check that the region end is larger than the region start. Values given: start: 2000, end: 1000
        """
        if len(args) == 1:
            chrom, start, end = GenomeRange.parse_region_string(args[0])
        elif len(args) == 3:
            chrom, start, end = args
        else:
            raise ValueError("inappropriate init arguments. "
                "correct example: `range1 = GenomeRange(\"chr1:1000-2000\")` or "
                "`range1 = GenomeRange(\"chr1\", 1000, 2000)`")

        if end <= start:
            raise ValueError("Please check that the region end is larger than the region start. "
                "Values given: start: {}, end: {}".format(start, end))

        self.chrom = chrom
        self.start = start
        self.end = end

    @staticmethod
    def parse_region_string(region_string):
        """
        splits a region string into
        a (chrom, start, end) tuple

        Parameters
        ----------
        region_string : str
            Region string to be parsed, like: "chr:start-end"

        Return
        ------
        result : tuple of str
            Result tuple (chrom, start, end)

        >>> GenomeRange.parse_region_string("chr1:10-20")
        ('chr1', 10, 20)
        >>> GenomeRange.parse_region_string("chr1:0")
        Traceback (innermost last):
         ...
        ValueError: Failure to parse region string, please check that region format should be like "chr:start-end".
        """
        if region_string:
            # separate the chromosome name and the location using the ':' character
            chrom, position = region_string.strip().split(":")

            # clean up the position
            for char in ",.;|!{}()":
                position = position.replace(char, '')

            position_list = position.split("-")
            try:
                region_start = int(position_list[0])
                region_end = int(position_list[1])
            except:
                raise ValueError("Failure to parse region string, please check that region format "
                    "should be like \"chr:start-end\".")

            return chrom, region_start, region_end
    
    def change_chrom_names(self):
        """
        >>> range1 = GenomeRange("chr1", 1000, 2000)
        >>> range1.chrom
        'chr1'
        >>> range1.change_chrom_names()
        >>> range1.chrom
        '1'
        >>> range1.change_chrom_names()
        >>> range1.chrom
        'chr1'
        """
        self.chrom = change_chrom_names(self.chrom)
    
    def __str__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end)

    @property
    def length(self):
        """
        >>> range1 = GenomeRange("chr1", 0, 1000)
        >>> range1.length
        1000
        """
        return self.end - self.start

    def __eq__(self, other):
        """
        >>> GenomeRange('chr1', 1000, 2000) == GenomeRange("chr1:1000-2000")
        True
        >>> GenomeRange("chr1:1000-2000") == GenomeRange("1:1000-2000")
        False
        """
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))


def op_err_msg(a, b, op='+'):
    """
    Generate the error message of error operand type.
    """
    msg = "unsupported operand type(s) for {}: '{}' and '{}'".format(op, str(type(a)), str(type(b)))
    return msg


def cm2inch(*tupl):
    """
    convert length unit from cm to inch. 

    >>> cm2inch(10)
    3.937007874015748
    >>> cm2inch(2, 2)
    (0.7874015748031495, 0.7874015748031495)
    >>> cm2inch((1, 5))
    (0.39370078740157477, 1.968503937007874)
    """
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        if len(tupl) != 1:
            return tuple(i / inch for i in tupl)
        else:
            return tupl[0] / inch


def opener(filename):
    """
    Determines if a file is compressed or not

    >>> import gzip
    >>> msg = "hello blablabla"
    >>> tmp_f_raw  = open("/tmp/test_opener.txt", 'w')
    >>> tmp_f_raw.write(msg)
    15
    >>> tmp_f_raw.close()
    >>> tmp_f_gzip = gzip.open('/tmp/test_opener.txt.gz', 'wb')
    >>> tmp_f_gzip.write(to_bytes(msg))
    15
    >>> tmp_f_gzip.close()

    >>> test_raw = opener(tmp_f_raw.name)
    >>> type(test_raw)
    <class '_io.BufferedReader'>
    >>> test_gzip = opener(tmp_f_gzip.name)
    >>> type(test_gzip)
    <class 'gzip.GzipFile'>

    >>> test_raw.close()
    >>> test_gzip.close()

    >>> import os
    >>> os.remove(test_raw.name)
    >>> os.remove(test_gzip.name)

    """
    import gzip
    f = open(filename, 'rb')
    if f.read(2) == b'\x1f\x8b':
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.seek(0)
        return f


def to_string(s):
    """
    Convert bytes, bytes list to string, string list.

    >>> to_string("hello")
    'hello'
    >>> to_string(b"hello")
    'hello'
    >>> to_string([b'hello', b'world'])
    ['hello', 'world']
    """
    if isinstance(s, str):
        return s
    if isinstance(s, bytes):
        return s.decode('ascii')
    if isinstance(s, list):
        return [to_string(x) for x in s]
    return s


def to_bytes(s):
    """
    Like toString, to bytes.

    >>> to_bytes('hello')
    b'hello'
    >>> to_bytes(['hello', 'world'])
    [b'hello', b'world']
    """
    if isinstance(s, bytes):
        return s
    if isinstance(s, str):
        return bytes(s, 'ascii')
    if isinstance(s, list):
        return [to_bytes(x) for x in s]
    return s


def rgb2hex(r, g, b):
    """
    Convert rgb color to hex format.

    >>> rgb2hex(129, 154, 70)
    '#819a46'
    >>> rgb2hex(-10, 256, -1)
    Traceback (most recent call last):
    ...
    AssertionError: (r, g, b) value must within range 0 ~ 255.
    """
    assert (0, 0, 0) <= (r, g, b) <= (255, 255, 255), \
        "(r, g, b) value must within range 0 ~ 255."
    hex = "#{:02x}{:02x}{:02x}".format(r,g,b)
    return hex


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
                "Bed file not sorted. Please use a sorted bed file.\n{}{} ".format(prev_line, line)

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


def change_chrom_names(chrom):
    """
    Changes UCSC chromosome names to ensembl chromosome names
    and vice versa.

    >>> change_chrom_names("chr1")
    '1'
    >>> change_chrom_names("1")
    'chr1'
    """
    # TODO: mapping from chromosome names like mithocondria is missing
    if chrom.startswith('chr'):
        # remove the chr part from chromosome name
        chrom = chrom[3:]
    else:
        # prefix with 'chr' the chromosome name
        chrom = 'chr' + chrom

    return chrom


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
                "Bed file not sorted. Please use a sorted bed file.\n" \
                "File: {}\n" \
                "Previous line: {}\n Current line{} ".format(self.file_handle.name, self.prev_line, line)

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
                "Bed file not sorted. Please use a sorted bed file.\n" \
                "File: {}\n" \
                "Previous line: {}\n Current line{} ".format(self.file_handle.name, self.prev_line, line)

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


class GenomeLength(dict):

    def __init__(self, length_file, genome_name=""):
        self.name = genome_name
        self.length_file = length_file
        self.parse_file(length_file)

    def parse_file(self, length_file):
        with opener(length_file) as f:
            for idx, line in enumerate(f):
                line = to_string(line)
                chrom, length, *_ = line.strip().split()
                try:
                    length = int(length)
                except ValueError as detail:
                    log.warning("Error reading line #{}. The field {} is not a integer.\n"
                                "Error message: {}\n".format(idx+1, length, detail))
                self.__setitem__(chrom, length)

    def check_range(self, genome_range):
        """
        Check a genome range is valid or not.

        (both start and end position is one based.)
        """
        if genome_range.chrom not in self:
            return False 
        if genome_range.start < 1:
            return False
        if genome_range.end > self[genome_range.chrom]:
            return False
        return True

    def bound_range(self, genome_range):
        """
        Bound a genome range within reference.

        (both start and end position is one based.) 
        """
        if self.check_range(genome_range):
            return genome_range
        else:
            chrom, start, end = "", 0, 0

            if genome_range.chrom not in self:
                raise ValueError("{} not in chromosome file: {}".format(
                    genome_range.chrom, self.length_file))
            else:
                chrom = genome_range.chrom

            if genome_range.start < 1:
                start = 1
            else:
                start = genome_range.start

            if genome_range.end > self[genome_range.chrom]:
                end = self[genome_range.chrom]
            else:
                end = genome_range.end

            bounded_range = GenomeRange(chrom, start, end)
            return bounded_range


def get_size(obj, seen=None):
    """
    Recursively finds size of objects

    From:
        https://stackoverflow.com/a/40880923/8500469
    """
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


def fig2bytes(fig, encode='png', dpi=None):
    """
    Convert matplotlib.figure.Figure object to image bytes.
    """
    buf = io.BytesIO()
    fig.savefig(buf, format=encode, dpi=None)
    buf.seek(0)
    img_bytes = buf.read()
    buf.close()
    return img_bytes


#
# read genome length file #

here = dirname(abspath(__file__))

HG19 = GenomeLength(join(here, "genome/hg19.txt"), genome_name="hg19")
HG38 = GenomeLength(join(here, "genome/hg38.txt"), genome_name="hg38")
MM9  = GenomeLength(join(here, "genome/mm9.txt"),  genome_name="mm9")
MM10 = GenomeLength(join(here, "genome/mm10.txt"), genome_name="mm10")

BUILT_IN_GENOMES = {
    'hg19': HG19,
    'hg38': HG38,
    'mm9':  MM9,
    'mm10': MM10,
}


#
# format convert #

_fields_rg = ("bin", "name", "chrom", "strand", "txStart", "txEnd",
              "cdsStart", "cdsEnd", "exonCount", "exonStart", "exonEnds",
              "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")

_refGeneRec = collections.namedtuple("refGeneRec", _fields_rg)


class refGeneRec(_refGeneRec):

    def to_bed12_line(self):
        chrom = self.chrom
        start = self.txStart
        end = self.txEnd
        name = self.name2
        score = self.score
        strand = self.strand

        thick_start = start
        thick_end = end

        item_rgb = "0"

        block_count = self.exonCount
        block_starts = self.offset_zero(self.exonStart, start) + ","
        block_ends = self.exonEnds
        block_sizes = self.get_exons_size(self.exonStart, block_ends) + ","

        bed12_item = (chrom, start, end, name, score, strand, thick_start, thick_end, item_rgb,
                      block_count, block_sizes, block_starts)

        bed12_line = "\t".join(bed12_item)
        return bed12_line

    def to_line(self):
        return "\t".join(self)

    @staticmethod
    def offset_zero(block, offset):
        offset = int(offset)
        block = [int(i) for i in block.split(",") if i]
        block = [str(i - offset) for i in block]
        return ",".join(block)

    @staticmethod
    def get_exons_size(exons_start, exons_end):
        starts = exons_start.split(",")
        starts = [int(i) for i in starts if i]
        ends = exons_end.split(",")
        ends = [int(i) for i in ends if i]
        sizes = [ends[i] - starts[i] for i, _ in enumerate(starts)]
        exons_size = ",".join([str(i) for i in sizes])
        return exons_size


def refgene_txt_to_bed12(txt_file, bed_file):
    with opener(txt_file) as f_in, open(bed_file, 'w') as f_out:
        for line in f_in:
            line = to_string(line)
            items = line.strip().split("\t")
            refg_rec = refGeneRec._make(items)
            out_line = refg_rec.to_bed12_line() + "\n"
            f_out.write(out_line)


LOG_LEVEL = logging.WARNING


def get_logger(name, file_=sys.stderr, level=LOG_LEVEL):
    FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s"
    formatter = logging.Formatter(fmt=FORMAT)
    if isinstance(file_, str):
        handler = logging.FileHandler(file_)
    else:
        handler = logging.StreamHandler(file_)
    handler.setFormatter(formatter)
    log = logging.getLogger(name)
    log.addHandler(handler)
    log.setLevel(level)
    return log

log = get_logger(__name__)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
