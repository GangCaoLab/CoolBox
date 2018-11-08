

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
        self.normalization = normalization
        self.binsize = binsize
        self.chromosomes, self.resolutions, self.masterindex, self.genome, self.metadata = self.__info()

    def fetch(self, genome_range1, genome_range2=None):
        """
        Return
        ------
        matrix : numpy.ndarray
        """
        from coolbox.utilities.genome import GenomeRange
        from .tools import infer_resolution
        from .straw import straw

        if genome_range2 is None:
            genome_range2 = genome_range1

        if isinstance(genome_range1, str):
            genome_range1 = GenomeRange(genome_range1)
        if isinstance(genome_range2, str):
            genome_range2 = GenomeRange(genome_range2)

        if genome_range1.chrom.startswith("chr"):
            genome_range1.change_chrom_names()
        if genome_range2.chrom.startswith("chr"):
            genome_range2.change_chrom_names()

        if self.binsize == 'auto':
            binsize = infer_resolution(genome_range1, self.resolutions)
        else:
            binsize = self.binsize

        chr1loc = str(genome_range1).replace('-', ':')
        chr2loc = str(genome_range2).replace('-', ':')
        straw_list = straw(self.normalization, self.hic_file, chr1loc, chr2loc, 'BP', binsize)
        matrix = self.__list_to_matrix(straw_list, genome_range1, genome_range2, binsize)
        return matrix

    def __list_to_matrix(self, straw_list, genome_range1, genome_range2, binsize):
        import numpy as np
        binlen1 = (genome_range1.length // binsize) + 1
        binlen2 = (genome_range2.length // binsize) + 1
        mat = np.zeros((binlen1, binlen2), dtype=np.float64)
        for loc1, loc2, c in zip(*straw_list):
            bin1id = (loc1 - genome_range1.start) // binsize
            bin2id = (loc1 - genome_range2.start) // binsize
            mat[bin1id, bin2id] = c
            mat[bin2id, bin1id] = c
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

