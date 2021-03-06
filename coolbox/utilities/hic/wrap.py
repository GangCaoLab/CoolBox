from coolbox.utilities.logtools import get_logger
from coolbox.utilities.genome import to_gr

log = get_logger(__name__)


class StrawWrap(object):
    """
    A wrap for straw Python API, for read .hic file.

    Parameters
    ----------
    path : str
        Path to the '.hic' file

    normalization : {False, 'VC', 'VC_SQRT', 'KR'}
        Method for the matrix normalization.
        default 'KR'

    binsize : {'auto', int}
        resolution of the data. for example 5000.
        'auto' for calculate resolution automatically.
        default 'auto'

    """

    def __init__(self, path, normalization='KR', binsize='auto'):
        self.hic_file = path
        if normalization == 'no' or normalization is False:
            normalization = 'NONE'
        elif normalization is True:
            normalization = 'KR'
        self.normalization = normalization
        self.binsize = binsize
        self.chromosomes, self.resolutions, self.masterindex, self.genome, self.metadata = self.__info()
        self.fetched_binsize = None

    def fetch(self, genome_range1, genome_range2=None):
        """
        Return
        ------
        matrix : numpy.ndarray
        """
        from coolbox.utilities.genome import GenomeRange

        if genome_range2 is None:
            genome_range2 = genome_range1
        genome_range1 = to_gr(genome_range1)
        genome_range2 = to_gr(genome_range2)

        if isinstance(genome_range1, str):
            genome_range1 = GenomeRange(genome_range1)
        if isinstance(genome_range2, str):
            genome_range2 = GenomeRange(genome_range2)

        if genome_range1.chrom.startswith("chr"):
            genome_range1.change_chrom_names()
        if genome_range2.chrom.startswith("chr"):
            genome_range2.change_chrom_names()

        binsize = self.infer_binsize(genome_range1)
        self.fetched_binsize = binsize  # expose fetched binsize

        straw_list = self.__fetch_straw_list(genome_range1, genome_range2, binsize)
        matrix = self.__list_to_matrix(straw_list, genome_range1, genome_range2, binsize)
        return matrix

    def infer_binsize(self, genome_range):
        from .tools import infer_resolution
        if self.binsize == 'auto':
            binsize = infer_resolution(genome_range, self.resolutions)
        else:
            binsize = self.binsize
        return binsize

    def __fetch_straw_list(self, genome_range1, genome_range2, binsize):
        from coolbox.utilities.hic.straw import straw
        chr1loc = str(genome_range1).replace('-', ':')
        chr2loc = str(genome_range2).replace('-', ':')
        slist = []
        try:
            slist = straw(self.normalization, self.hic_file, chr1loc, chr2loc, 'BP', binsize)
        except Exception as e:
            log.warning("Error occurred when reading the dothic file with straw:")
            log.warning(str(e))
            if self.normalization != "NONE":
                log.warning("Try to read unbalanced matrix.")
                self.normalization = "NONE"
                try:
                    slist = straw(self.normalization, self.hic_file, chr1loc, chr2loc, 'BP', binsize)
                except Exception as e:
                    log.warning("Failed.")
                    log.warning(str(e))
                log.warning("Unbalanced matrix is readed.")
        return slist

    def fetch_pixels(self, genome_range1, genome_range2=None):
        from pandas import DataFrame
        if genome_range2 is None:
            genome_range2 = genome_range1
        if genome_range1.chrom.startswith("chr"):
            genome_range1.change_chrom_names()
        if genome_range2.chrom.startswith("chr"):
            genome_range2.change_chrom_names()
        binsize = self.infer_binsize(genome_range1)
        slist = self.__fetch_straw_list(genome_range1, genome_range2, binsize)
        pixels = DataFrame(slist)
        pixels = pixels.T
        pixels.columns = ['start1', 'start2', 'value']
        pixels['start1'] = pixels['start1'].astype('int32')
        pixels['start2'] = pixels['start2'].astype('int32')
        pixels['end1'] = pixels['start1'] + binsize
        pixels['end2'] = pixels['start2'] + binsize
        pixels['chrom1'] = genome_range1.chrom
        pixels['chrom2'] = genome_range2.chrom
        pixels = pixels[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'value']]
        return pixels

    def __list_to_matrix(self, straw_list, genome_range1, genome_range2, binsize):
        import numpy as np
        binlen1 = (genome_range1.length // binsize) + 1
        binlen2 = (genome_range2.length // binsize) + 1
        mat = np.zeros((binlen1, binlen2), dtype=np.float64)
        is_cis = (genome_range1 == genome_range2)
        for rec in zip(*straw_list):
            loc1, loc2, c = rec[0], rec[1], rec[2]
            bin1id = (loc1 - genome_range1.start) // binsize
            bin2id = (loc2 - genome_range2.start) // binsize
            if is_cis:
                mat[bin1id, bin2id] = c
                mat[bin2id, bin1id] = c
            else:
                if (0 <= bin1id < mat.shape[0]) and (0 <= bin2id < mat.shape[1]):
                    mat[bin1id, bin2id] = c
        return mat

    def __info(self):
        """
        from hic2cool code:
            https://github.com/4dn-dcic/hic2cool/blob/master/hic2cool/hic2cool_utils.py#L73

        Takes in a .hic file and returns a dictionary containing information about
        the chromosome. Keys are chromosome index numbers (0 through # of chroms
        contained in file) and values are [chr idx (int), chr name (str), chrom
        length (str)]. Returns the masterindex used by the file as well as the open
        file object.
        """
        import sys
        import struct
        from .straw import readcstr

        req = open(self.hic_file, 'rb')

        chrs = {}
        resolutions = []
        magic_string = struct.unpack(b'<3s', req.read(3))[0]
        req.read(1)
        if (magic_string != b"HIC"):
            print('This does not appear to be a HiC file; '
                  'magic string is incorrect')
            sys.exit()
        global version
        version = struct.unpack(b'<i', req.read(4))[0]
        masterindex = struct.unpack(b'<q', req.read(8))[0]
        genome = b""
        c = req.read(1)
        while (c != b'\0'):
            genome += c
            c = req.read(1)
        genome = genome.decode('ascii')
        # metadata extraction
        metadata = {}
        nattributes = struct.unpack(b'<i', req.read(4))[0]
        for x in range(nattributes):
            key = readcstr(req)
            value = readcstr(req)
            metadata[key] = value
        nChrs = struct.unpack(b'<i', req.read(4))[0]
        for i in range(0, nChrs):
            name = readcstr(req)
            length = struct.unpack(b'<i', req.read(4))[0]
            if name and length:
                chrs[i] = [i, name, length]
        nBpRes = struct.unpack(b'<i', req.read(4))[0]
        # find bp delimited resolutions supported by the hic file
        for x in range(0, nBpRes):
            res = struct.unpack(b'<i', req.read(4))[0]
            resolutions.append(res)
        return chrs, resolutions, masterindex, genome, metadata


class CoolerWrap(object):
    """
    wrap for cooler file,
    deal with multi resolution.

    Parameters
    ----------
    path : str
        Path to the cooler file.

    binsize : {'auto', int}
        resolution of the data. for example 5000.
        'auto' for calculate resolution automatically.
        default 'auto'

    balance : bool
        Balance the matrix or not.
        default True
    """

    def __init__(self, path, binsize='auto', balance=True):
        import cooler
        from .tools import is_multi_cool, get_cooler_resolutions
        self.path = path

        self.is_multi = is_multi_cool(path)
        self.resolutions = get_cooler_resolutions(path, self.is_multi)
        if self.is_multi:
            self.coolers = self.__load_multi_coolers(path)
        else:
            self.cool = cooler.Cooler(path)

        self.binsize = binsize
        self.balance = balance

    def __load_multi_coolers(self, path):
        from collections import OrderedDict
        import cooler

        def cooler_reso(resolution):
            from h5py import File
            with File(path, 'r') as f:
                if "resolutions" in f:
                    c = cooler.Cooler(path + "::/resolutions/{}".format(resolution))
                else:
                    for grp_name in f:
                        grp = f[grp_name]
                        if str(grp.attrs['bin-size']) == str(resolution):
                            c = cooler.Cooler(path + "::{}".format(grp_name))
                            break
            return c

        coolers = OrderedDict([(reso, cooler_reso(reso)) for reso in self.resolutions])
        return coolers

    def get_cool(self, genome_range):
        binsize = self.infer_binsize(genome_range)
        if self.is_multi:
            cool = self.coolers[binsize]
        else:
            cool = self.cool
        self.fetched_binsize = binsize  # expose fetched binsize
        return cool

    def infer_binsize(self, genome_range):
        from .tools import infer_resolution
        genome_range = to_gr(genome_range)
        if self.is_multi:
            resolutions = [k for k in self.coolers.keys()]
            if self.binsize == 'auto':
                binsize = infer_resolution(genome_range, resolutions)
            else:
                assert self.binsize in resolutions, \
                    "Multi-Cooler file not contain the resolution {}.".format(self.binsize)
                binsize = int(self.binsize)
        else:
            binsize = self.cool.binsize
        return binsize

    def fetch(self, genome_range1, genome_range2=None):
        # TODO what if genome_ranges are invalid
        if genome_range2 is None:
            genome_range2 = genome_range1

        genome_range1 = to_gr(genome_range1)
        genome_range2 = to_gr(genome_range2)

        cool = self.get_cool(genome_range1)

        if genome_range1.chrom not in cool.chromnames:
            genome_range1.change_chrom_names()
        if genome_range2.chrom not in cool.chromnames:
            genome_range2.change_chrom_names()

        try:
            mat = cool.matrix(balance=self.balance).fetch(str(genome_range1), str(genome_range2))
        except ValueError as e:
            log.warning(str(e))
            log.warning("Data is not balanced, force to use unbalanced matrix.")
            mat = cool.matrix(balance=False).fetch(str(genome_range1), str(genome_range2))

        return mat

    def fetch_pixels(self, genome_range1, genome_range2=None, join=True):
        cool = self.get_cool(genome_range1)

        if genome_range2 is None:
            genome_range2 = genome_range1

        if genome_range1.chrom not in cool.chromnames:
            genome_range1.change_chrom_names()
        if genome_range2.chrom not in cool.chromnames:
            genome_range2.change_chrom_names()

        mat = cool.matrix(as_pixels=True, balance=self.balance, join=join)
        pixels = mat.fetch(str(genome_range1), str(genome_range2))

        return pixels
