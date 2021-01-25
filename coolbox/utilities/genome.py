from typing import Iterator, Union
from .logtools import get_logger
from .filetool import opener, to_string

log = get_logger(__name__)


def to_gr(obj):
    """
    Convert object to GenomeRange.

    >>> gr = to_gr("chr1:0-1000")
    >>> (gr.chrom, gr.start, gr.end)
    ("chr1", 0, 1000)

    """
    if type(obj) is GenomeRange:
        return obj
    else:
        return GenomeRange(obj)


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
        >>> range3 = GenomeRange(range2)
        >>> tuple(range3)
        >>> range4 = GenomeRange("chr1", 2000, 1000)
        Traceback (most recent call last):
        ...
        ValueError: Please check that the region end is larger than the region start. Values given: start: 2000, end: 1000
        """
        if len(args) == 1:
            if isinstance(args[0], GenomeRange):
                chrom, start, end = tuple(args[0])
            # str format
            else:
                chrom, start, end = GenomeRange.parse_region_string(args[0])

        elif len(args) == 3:
            chrom, start, end = args
        else:
            raise ValueError("inappropriate init arguments. "
                             "correct example: `range1 = GenomeRange(\"chr1:1000-2000\")` or "
                             "`range1 = GenomeRange(\"chr1\", 1000, 2000)`")

        if end < start:
            raise ValueError("Please check that the region end is larger than the region start. "
                             "Values given: start: {}, end: {}".format(start, end))

        self.chrom = chrom
        self.start = start
        self.end = end

    def __iter__(self) -> Iterator[Union[str, int]]:
        yield self.chrom
        yield self.start
        yield self.end

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
        try:
            # separate the chromosome name and the location using the ':' character
            chrom, position = region_string.strip().split(":")

            # clean up the position
            for char in ",.;|!{}()":
                position = position.replace(char, '')

            position_list = position.split("-")
            region_start = int(position_list[0])
            region_end = int(position_list[1])
            return chrom, region_start, region_end

        except Exception:
            raise ValueError(f"Failure to parse region string {region_string}, please check that region format "
                             f"should be like \"chr:start-end\".")

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

    def __contains__(self, another):
        if another.chrom != self.chrom:
            return False
        if another.start < self.start:
            return False
        if another.end > self.end:
            return False
        return True

    def __str__(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"

    def __repr__(self) -> str:
        return f"GenomeRange('{self}')"

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
                                "Error message: {}\n".format(idx + 1, length, detail))
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
                if start >= end:
                    start = end - genome_range.length
                    if start < 1:
                        start = 1
            else:
                end = genome_range.end

            bounded_range = GenomeRange(chrom, start, end)
            return bounded_range


def split_genome_range(genome_range):
    """
    Little helper func.
    enforce gr is a GenomeRange object, and
    split gr to chrom, start, end
    """
    if isinstance(genome_range, str):
        genome_range = GenomeRange(genome_range)
    else:
        assert isinstance(genome_range, GenomeRange), \
            "gr is a `GenomeRange` object or a genome range str"
    chrom, start, end = genome_range.chrom, genome_range.start, genome_range.end
    return chrom, start, end
