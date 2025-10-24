"""
Straw module

Straw enables programmatic access to .hic files.
.hic files store the contact matrices from Hi-C experiments and the
normalization and expected vectors, along with meta-data in the header.

Usage: strawObj = straw <hicFile(s)>
       matrixObj = strawObj.getNormalizedMatrix <chr1> <chr2> <NONE/VC/VC_SQRT/KR> <BP/FRAG> <binsize>
       data = matrixObj.getDataFromBinRegion <x1,x2,y1,y2>

Example:
   import straw
   strawObj = straw(filename)
   matrixObj = strawObj.getNormalizedMatrix('5', '5', 'KR', 'BP', 5000)
   result = matrixObj.getDataFromBinRegion(0,500,0,500)
   for i in range(len(result[0])):
...   print("{0}\t{1}\t{2}".format(result[0][i], result[1][i], result[2][i]))

See https://github.com/theaidenlab/straw/wiki/Python for more documentation
"""
from __future__ import absolute_import, division, print_function, unicode_literals

__author__ = "Yue Wu, Neva Durand, Yossi Eliaz, Muhammad Shamim, Erez Aiden"
__license__ = "MIT"

import struct
import zlib
import requests
import io
import concurrent.futures
import math
import sys


def __readcstr(f):
    """ Helper function for reading in C-style string from file
    """
    buf = b""
    while True:
        b = f.read(1)
        if b is None or b == b"\0":
            return buf.decode("utf-8")
        elif b == "":
            raise EOFError("Buffer unexpectedly empty while trying to read null-terminated string")
        else:
            buf += b


readcstr = __readcstr


"""
functions for chrom.sizes
internal representation is a dictionary with 
chromosome name as the key
value maps to a tuple containing the index and chromosome length
"""


class ChromDotSizes:
    def __init__(self, data):
        self.data = data

    def getLength(self, chrom):
        try:
            return int(self.data[chrom][1])
        except:
            print(str(chrom) + " not in chrom.sizes. Check that the chromosome name matches the genome.\n")
            return None

    def getIndex(self, chrom):
        try:
            return int(self.data[chrom][0])
        except:
            print(str(chrom) + " not in chrom.sizes. Check that the chromosome name matches the genome.\n")
            return None

    def figureOutEndpoints(self, chrAndPositions):
        chrAndPositionsArray = chrAndPositions.split(":")
        chrom = chrAndPositionsArray[0]

        indx1 = 0
        indx2 = self.getLength(chrom)

        if len(chrAndPositionsArray) == 3:
            indx1 = int(chrAndPositionsArray[1])
            indx2 = int(chrAndPositionsArray[2])

        return chrom, indx1, indx2


def read_metadata(infile, verbose=False):
    """
    Reads the metadata of HiC file from header.

    Args
    infile: str, path to the HiC file
    verbose: bool

    Returns
    metadata: dict, containing the metadata.
                Keys of the metadata:
                HiC version,
                Master index,
                Genome ID (str),
                Attribute dictionary (dict),
                Chromosomes (dict),
                Base pair-delimited resolutions (list),
                Fragment-delimited resolutions (list).
    """
    metadata = {}
    import io
    import struct
    if (infile.startswith("http")):
        # try URL first. 100K should be sufficient for header
        headers = {'range': 'bytes=0-100000', 'x-amz-meta-requester': 'straw'}
        s = requests.Session()
        r = s.get(infile, headers=headers)
        if (r.status_code >= 400):
            print("Error accessing " + infile)
            print("HTTP status code " + str(r.status_code))
            sys.exit(1)
        req = io.BytesIO(r.content)
        myrange = r.headers['content-range'].split('/')
        totalbytes = myrange[1]
    else:
        req = open(infile, 'rb')
    magic_string = struct.unpack('<3s', req.read(3))[0]
    req.read(1)
    if (magic_string != b"HIC"):
        sys.exit('This does not appear to be a HiC file magic string is incorrect')
    version = struct.unpack('<i', req.read(4))[0]
    metadata['HiC version'] = version
    masterindex = struct.unpack('<q', req.read(8))[0]
    metadata['Master index'] = masterindex
    genome = ""
    c = req.read(1).decode("utf-8")
    while (c != '\0'):
        genome += c
        c = req.read(1).decode("utf-8")
    metadata['Genome ID'] = genome
    if (version > 8):
        nvi = struct.unpack('<q', req.read(8))[0]
        nvisize = struct.unpack('<q', req.read(8))[0]
        metadata['NVI'] = nvi
        metadata['NVI size'] = nvisize
    ## read and throw away attribute dictionary (stats+graphs)
    nattributes = struct.unpack('<i', req.read(4))[0]
    d = {}
    for x in range(0, nattributes):
        key = __readcstr(req)
        value = __readcstr(req)
        d[key] = value
    metadata['Attribute dictionary'] = d
    nChrs = struct.unpack('<i', req.read(4))[0]
    d = {}
    for x in range(0, nChrs):
        key = __readcstr(req)
        if (version > 8):
            value = struct.unpack('q', req.read(8))[0]
        else:
            value = struct.unpack('<i', req.read(4))[0]
        d[key] = value
    metadata["Chromosomes"] = d
    nBpRes = struct.unpack('<i', req.read(4))[0]
    l = []
    for x in range(0, nBpRes):
        res = struct.unpack('<i', req.read(4))[0]
        l.append(res)
    metadata["Base pair-delimited resolutions"] = l
    nFrag = struct.unpack('<i', req.read(4))[0]
    l = []
    for x in range(0, nFrag):
        res = struct.unpack('<i', req.read(4))[0]
        l.append(res)
    metadata["Fragment-delimited resolutions"] = l
    for k in metadata:
        if k != 'Attribute dictionary':
            print(k, ':', metadata[k])
    if verbose:
        print('Attribute dictionary', ':', metadata['Attribute dictionary'])
    return metadata


def readHeader(infile, is_synapse):
    """ Reads the header

    Args:
       input file, is_synapse

    Returns:
       list: master index, version number, size of totalbytes, chromDotSizes
    """

    if infile.startswith("http"):
        # try URL first. 100K should be sufficient for header
        headers = getHttpHeader('bytes=0-100000', is_synapse)
        s = requests.Session()
        r = s.get(infile, headers=headers)
        if r.status_code >= 400:
            print("Error accessing " + infile)
            print("HTTP status code " + str(r.status_code))
            return -1
        req = io.BytesIO(r.content)
        myrange = r.headers['content-range'].split('/')
        totalbytes = myrange[1]
    else:
        req = open(infile, 'rb')
        totalbytes = None

    magic_string = struct.unpack('<3s', req.read(3))[0]
    req.read(1)
    if magic_string != b"HIC":
        print('This does not appear to be a HiC file magic string is incorrect')
        return -1
    version = struct.unpack('<i', req.read(4))[0]
    if version < 6:
        print("Version {0} no longer supported".format(str(version)))
        return -1
    #print('HiC version:' + '  {0}'.format(str(version)))
    master = struct.unpack('<q', req.read(8))[0]
    genome = b""
    c = req.read(1)
    while c != b'\0':
        genome += c
        c = req.read(1)

    # read and throw away attribute dictionary (stats+graphs)
    nattributes = struct.unpack('<i', req.read(4))[0]
    for x in range(nattributes):
        key = __readcstr(req)
        value = __readcstr(req)
    nChrs = struct.unpack('<i', req.read(4))[0]
    chromDotSizes = {}
    for i in range(0, nChrs):
        name = __readcstr(req)
        length = struct.unpack('<i', req.read(4))[0]
        chromDotSizes[name] = (i, length)
    return master, version, totalbytes, ChromDotSizes(chromDotSizes)


def readFooter(infile, is_synapse, master, totalbytes):
    """Reads the footer, which contains all the expected and normalization
    vectors. Presumes file pointer is in correct position
    Args:
       req (file): File to read from; presumes file pointer is in correct
       position
       chr1 (str): Chromosome 1
       chr2 (str): Chromosome 2
       norm (str): Normalization type, one of NONE, VC, KR, VC_SQRT
       unit (str): One of BP or FRAG
       resolution (int): Bin size

    Returns:
       list: File position of matrix, position+size chr1 normalization vector,
             position+size chr2 normalization vector
    """
    if infile.startswith("http"):
        headers = getHttpHeader('bytes={0}-{1}'.format(master, totalbytes), is_synapse)
        s = requests.Session()
        r = s.get(infile, headers=headers)
        req = io.BytesIO(r.content)
    else:
        req = open(infile, 'rb')
        req.seek(master)

    filePositions = dict()
    nBytes = struct.unpack('<i', req.read(4))[0]
    nEntries = struct.unpack('<i', req.read(4))[0]

    for i in range(nEntries):
        key = __readcstr(req)
        fpos = struct.unpack('<q', req.read(8))[0]
        sizeinbytes = struct.unpack('<i', req.read(4))[0]
        filePositions[key] = (fpos, sizeinbytes)

    # later save these
    nExpectedValues = struct.unpack('<i', req.read(4))[0]
    for i in range(nExpectedValues):
        key = __readcstr(req)
        binSize = struct.unpack('<i', req.read(4))[0]
        nValues = struct.unpack('<i', req.read(4))[0]
        for j in range(nValues):
            # replace with vector.append
            v = struct.unpack('<d', req.read(8))[0]
        nNormalizationFactors = struct.unpack('<i', req.read(4))[0]
        for j in range(nNormalizationFactors):
            # replace with vector.append
            chrIdx = struct.unpack('<i', req.read(4))[0]
            v = struct.unpack('<d', req.read(8))[0]
    nExpectedValues = struct.unpack('<i', req.read(4))[0]
    for i in range(nExpectedValues):
        str_ = __readcstr(req)
        str_ = __readcstr(req)
        binSize = struct.unpack('<i', req.read(4))[0]
        nValues = struct.unpack('<i', req.read(4))[0]
        for j in range(nValues):
            v = struct.unpack('<d', req.read(8))[0]
        nNormalizationFactors = struct.unpack('<i', req.read(4))[0]
        for j in range(nNormalizationFactors):
            chrIdx = struct.unpack('<i', req.read(4))[0]
            v = struct.unpack('<d', req.read(8))[0]

    normMap = dict()
    nEntries = struct.unpack('<i', req.read(4))[0]
    for i in range(nEntries):
        normtype = __readcstr(req)
        if normtype not in normMap:
            normMap[normtype] = {}
        chrIdx = struct.unpack('<i', req.read(4))[0]
        if chrIdx not in normMap[normtype]:
            normMap[normtype][chrIdx] = {}
        unit = __readcstr(req)
        if unit not in normMap[normtype][chrIdx]:
            normMap[normtype][chrIdx][unit] = {}
        resolution = struct.unpack('<i', req.read(4))[0]
        if resolution not in normMap[normtype][chrIdx][unit]:
            normMap[normtype][chrIdx][unit][resolution] = {}
        filePosition = struct.unpack('<q', req.read(8))[0]
        sizeInBytes = struct.unpack('<i', req.read(4))[0]

        normMap[normtype][chrIdx][unit][resolution]['position'] = filePosition
        normMap[normtype][chrIdx][unit][resolution]['size'] = sizeInBytes

    return filePositions, normMap


def readMatrixZoomData(req, myunit, mybinsize, blockMap):
    """ Reads the Matrix Zoom Data, which gives pointer list for blocks for
    the data. Presumes file pointer is in correct position

    Args:
       req (file): File to read from; presumes file pointer is in correct
       position
       myunit (str): Unit (BP or FRAG) we're searching for
       mybinsize (int): Resolution we're searching for

    Returns:
       list containing boolean indicating if we found appropriate matrix,
       and if so, the counts for the bins and columns
    """
    unit = __readcstr(req)
    temp = struct.unpack('<i', req.read(4))[0]
    temp = struct.unpack('<f', req.read(4))[0]
    temp = struct.unpack('<f', req.read(4))[0]
    temp = struct.unpack('<f', req.read(4))[0]
    temp = struct.unpack('<f', req.read(4))[0]
    binSize = struct.unpack('<i', req.read(4))[0]
    blockBinCount = struct.unpack('<i', req.read(4))[0]
    blockColumnCount = struct.unpack('<i', req.read(4))[0]
    storeBlockData = False
    # for the initial
    myBlockBinCount = -1
    myBlockColumnCount = -1
    if myunit == unit and mybinsize == binSize:
        myBlockBinCount = blockBinCount
        myBlockColumnCount = blockColumnCount
        storeBlockData = True
    nBlocks = struct.unpack('<i', req.read(4))[0]
    for b in range(nBlocks):
        blockNumber = struct.unpack('<i', req.read(4))[0]
        filePosition = struct.unpack('<q', req.read(8))[0]
        blockSizeInBytes = struct.unpack('<i', req.read(4))[0]
        entry = dict()
        entry['size'] = blockSizeInBytes
        entry['position'] = filePosition
        if storeBlockData:
            blockMap[blockNumber] = entry
    return storeBlockData, myBlockBinCount, myBlockColumnCount


def readMatrix(req, unit, binsize, blockMap):
    """ Reads the matrix - that is, finds the appropriate pointers to block
    data and stores them. Needs to read through headers of zoom data to find
    appropriate matrix. Presumes file pointer is in correct position.

    Args:
       req (file): File to read from; presumes file pointer is in correct
       position
       unit (str): Unit to search for (BP or FRAG)
       binsize (int): Resolution to search for

    Returns:
       list containing block bin count and block column count of matrix

    Raises:
       ValueError if the .hic file can't be parsed with the specified resolution (binsize)
    """
    c1 = struct.unpack('<i', req.read(4))[0]
    c2 = struct.unpack('<i', req.read(4))[0]
    nRes = struct.unpack('<i', req.read(4))[0]
    i = 0
    found = False
    blockBinCount = -1
    blockColumnCount = -1
    while i < nRes and (not found):
        found, blockBinCount, blockColumnCount = readMatrixZoomData(req, unit, binsize, blockMap)
        i = i + 1
    if not found:
        raise ValueError(f"Error: could not parse .hic file using specified resolution/bin-size ({binsize})")

    return blockBinCount, blockColumnCount


def getBlockNumbersForRegionFromBinPosition(regionIndices, blockBinCount, blockColumnCount, intra):
    """ Gets the block numbers we will need for a specific region; used when
    the range to extract is sent in as a parameter

    Args:
       regionIndices (array): Array of ints giving range
       blockBinCount (int): The block bin count of the matrix
       blockColumnCount (int): The block column count of the matrix
       intra: Flag indicating if this is an intrachromosomal matrix

    Returns:
       blockSet (set): A set of blocks to print
    """
    col1 = int(regionIndices[0] / blockBinCount)
    col2 = int((regionIndices[1] + 1) / blockBinCount)
    row1 = int(regionIndices[2] / blockBinCount)
    row2 = int((regionIndices[3] + 1) / blockBinCount)
    blocksSet = set()

    for r in range(row1, row2 + 1):
        for c in range(col1, col2 + 1):
            blockNumber = r * blockColumnCount + c
            blocksSet.add(blockNumber)
    # in Java code, this is "if getBelowDiagonal"
    if intra and col2 > row1:
        for r in range(col1, col2 + 1):
            for c in range(row1, row2 + 1):
                blockNumber = r * blockColumnCount + c
                blocksSet.add(blockNumber)
    return blocksSet


def readBlock(req, size, version):
    """ Reads the block - reads the compressed bytes, decompresses, and stores
    results in array. Presumes file pointer is in correct position.

    Args:
       req (file): File to read from. Presumes file pointer is in correct
       position
       size (int): How many bytes to read

    Returns:
       array containing row, column, count data for this block
    """
    compressedBytes = req.read(size)
    uncompressedBytes = zlib.decompress(compressedBytes)
    nRecords = struct.unpack('<i', uncompressedBytes[0:4])[0]
    v = []
    if version < 7:
        for i in range(nRecords):
            binX = struct.unpack('<i', uncompressedBytes[(12 * i + 4):(12 * i + 8)])[0]
            binY = struct.unpack('<i', uncompressedBytes[(12 * i + 8):(12 * i + 12)])[0]
            counts = struct.unpack('<f', uncompressedBytes[(12 * i + 12):(12 * i + 16)])[0]
            record = dict()
            record['binX'] = binX
            record['binY'] = binY
            record['counts'] = counts
            v.append(record)
    else:
        binXOffset = struct.unpack('<i', uncompressedBytes[4:8])[0]
        binYOffset = struct.unpack('<i', uncompressedBytes[8:12])[0]
        useShort = struct.unpack('<b', uncompressedBytes[12:13])[0]
        type_ = struct.unpack('<b', uncompressedBytes[13:14])[0]
        index = 0
        if type_ == 1:
            rowCount = struct.unpack('<h', uncompressedBytes[14:16])[0]
            temp = 16
            for i in range(rowCount):
                y = struct.unpack('<h', uncompressedBytes[temp:(temp + 2)])[0]
                temp = temp + 2
                binY = y + binYOffset
                colCount = struct.unpack('<h', uncompressedBytes[temp:(temp + 2)])[0]
                temp = temp + 2
                for j in range(colCount):
                    x = struct.unpack('<h', uncompressedBytes[temp:(temp + 2)])[0]
                    temp = temp + 2
                    binX = binXOffset + x
                    if useShort == 0:
                        c = struct.unpack('<h', uncompressedBytes[temp:(temp + 2)])[0]
                        temp = temp + 2
                        counts = c
                    else:
                        counts = struct.unpack('<f', uncompressedBytes[temp:(temp + 4)])[0]
                        temp = temp + 4
                    record = dict()
                    record['binX'] = binX
                    record['binY'] = binY
                    record['counts'] = counts
                    v.append(record)
                    index = index + 1
        elif type_ == 2:
            temp = 14
            nPts = struct.unpack('<i', uncompressedBytes[temp:(temp + 4)])[0]
            temp = temp + 4
            w = struct.unpack('<h', uncompressedBytes[temp:(temp + 2)])[0]
            temp = temp + 2
            for i in range(nPts):
                row = int(i / w)
                col = i - row * w
                bin1 = int(binXOffset + col)
                bin2 = int(binYOffset + row)
                if useShort == 0:
                    c = struct.unpack('<h', uncompressedBytes[temp:(temp + 2)])[0]
                    temp = temp + 2
                    if c != -32768:
                        record = dict()
                        record['binX'] = bin1
                        record['binY'] = bin2
                        record['counts'] = c
                        v.append(record)
                        index = index + 1
                else:
                    counts = struct.unpack('<f', uncompressedBytes[temp:(temp + 4)])[0]
                    temp = temp + 4
                    if counts != 0x7fc00000:
                        record = dict()
                        record['binX'] = bin1
                        record['binY'] = bin2
                        record['counts'] = counts
                        v.append(record)
                        index = index + 1
    return v


def readBlockWorker(infile, is_synapse, blockNum, binsize, blockMap, norm, c1Norm, c2Norm, binPositionBox, isIntra,
                    version):
    yActual = []
    xActual = []
    counts = []
    idx = dict()
    if blockNum in blockMap:
        idx = blockMap[blockNum]
    else:
        idx['size'] = 0
        idx['position'] = 0

    if idx['size'] == 0:
        records = []
    else:
        if infile.startswith("http"):
            headers = getHttpHeader('bytes={0}-{1}'.format(idx['position'], idx['position'] + idx['size']), is_synapse)
            s = requests.Session()
            r = s.get(infile, headers=headers);
            req = io.BytesIO(r.content)
        else:
            req = open(infile, 'rb')
            req.seek(idx['position'])
        records = readBlock(req, idx['size'], version)

    # No caching currently; in Java code we keep all records and check positions later
    if norm != "NONE":
        for record in records:
            binX = record['binX']
            binY = record['binY']

            if ((binPositionBox[0] <= binX <= binPositionBox[1] and binPositionBox[2] <= binY <=
                 binPositionBox[3]) or (
                    isIntra and binPositionBox[0] <= binY <= binPositionBox[1] and binPositionBox[2] <= binX <=
                    binPositionBox[3])):
                c = record['counts']
                a = c1Norm[binX] * c2Norm[binY]
                if a != 0.0:
                    c = (c / a)
                else:
                    c = "inf"
                xActual.append(binX)
                yActual.append(binY)
                counts.append(c)
    else:
        for record in records:
            binX = record['binX']
            binY = record['binY']
            if ((binPositionBox[0] <= binX <= binPositionBox[1] and binPositionBox[2] <= binY <=
                 binPositionBox[3]) or (
                    isIntra and binPositionBox[0] <= binY <= binPositionBox[1] and binPositionBox[2] <= binX <=
                    binPositionBox[3])):
                c = record['counts']
                xActual.append(binX)
                yActual.append(binY)
                counts.append(c)
    return xActual, yActual, counts


def readNormalizationVector(req):
    """ Reads the normalization vector from the file; presumes file pointer is
    in correct position

    Args:
       req (file): File to read from; presumes file pointer is in correct
       position

    Returns:
      Array of normalization values

    """
    value = []
    nValues = struct.unpack('<i', req.read(4))[0]
    for i in range(nValues):
        d = struct.unpack('<d', req.read(8))[0]
        value.append(d)
    return value


def getHttpHeader(endrange, is_synapse):
    if is_synapse:
        return {'range': endrange}
    return {'range': endrange, 'x-amz-meta-requester': 'straw'}


def readLocalNorm(infile, position):
    req = open(infile, 'rb')
    req.seek(position)
    return readNormalizationVector(req)


def readHttpNorm(infile, normEntry, is_synapse):
    endrange = 'bytes={0}-{1}'.format(normEntry['position'], normEntry['position'] + normEntry['size'])
    headers = getHttpHeader(endrange, is_synapse)
    s = requests.Session()
    r = s.get(infile, headers=headers);
    req = io.BytesIO(r.content);
    return readNormalizationVector(req)


class straw:
    def __init__(self, infile, is_synapse=False):
        """ This is the main workhorse method of the module. Reads a .hic file and
        extracts the given contact matrix. Stores in an array in sparse upper
        triangular format: row, column, (normalized) count

        Args:
           norm(str): Normalization type, one of VC, KR, VC_SQRT, or NONE
           infile(str): File name or URL of .hic file
           chr1loc(str): Chromosome name and (optionally) range, i.e. "1" or "1:10000:25000"
           chr2loc(str): Chromosome name and (optionally) range, i.e. "1" or "1:10000:25000"
           unit(str): One of BP or FRAG
           binsize(int): Resolution, i.e. 25000 for 25K
        """

        self.isHttpFile = infile.startswith("http")
        self.infile = infile
        self.is_synapse = is_synapse
        self.master, self.version, totalbytes, self.chromDotSizes = readHeader(infile, is_synapse)
        self.myFilePositions, self.normMap = readFooter(infile, is_synapse, self.master, totalbytes)

    def getNormalizedMatrix(self, chr1, chr2, norm, unit, binsize):

        if not (unit == "BP" or unit == "FRAG"):
            print(
                "Unit specified incorrectly, must be one of <BP/FRAG>\nUsage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>\n")
            return None

        for chrom in [chr1, chr2]:
            if chrom not in self.chromDotSizes.data:
                print(str(chrom) + " wasn't found in the file. Check that the chromosome name matches the genome.\n")
                return None

        chrIndex1 = self.chromDotSizes.getIndex(chr1)
        chrIndex2 = self.chromDotSizes.getIndex(chr2)
        isIntra = chrIndex1 == chrIndex2

        neededToFlipIndices = False
        if chrIndex1 > chrIndex2:
            neededToFlipIndices = True
            chrIndex1, chrIndex2 = chrIndex2, chrIndex1
            chr1, chr2 = chr2, chr1

        executor = concurrent.futures.ThreadPoolExecutor()
        if norm != "NONE":
            try:
                c1NormEntry = self.normMap[norm][chrIndex1][unit][binsize]
            except:
                print(
                    "File did not contain {0} norm vectors for chr {1} at {2} {3}\n".format(norm, chr1, binsize, unit))
                return None

            if not isIntra:
                try:
                    c2NormEntry = self.normMap[norm][chrIndex2][unit][binsize]
                except:
                    print("File did not contain {0} norm vectors for chr {1} at {2} {3}\n".format(norm, chr2, binsize,
                                                                                                  unit))
                    return None
            if self.isHttpFile:
                futureNorm1 = executor.submit(readHttpNorm, self.infile, c1NormEntry, self.is_synapse)
                if not isIntra:
                    futureNorm2 = executor.submit(readHttpNorm, self.infile, c2NormEntry, self.is_synapse)
            else:
                futureNorm1 = executor.submit(readLocalNorm, self.infile, c1NormEntry['position'])
                if not isIntra:
                    futureNorm2 = executor.submit(readLocalNorm, self.infile, c2NormEntry['position'])

        blockMap = dict()
        key = str(chrIndex1) + "_" + str(chrIndex2)
        if key not in self.myFilePositions:
            print("File doesn't have the given {0} map\n".format(key))
            return None
        myFilePos = self.myFilePositions[key][0]
        if self.isHttpFile:
            headers = getHttpHeader('bytes={0}-'.format(myFilePos), self.is_synapse)
            s = requests.Session()
            r = s.get(self.infile, headers=headers, stream=True)
            futureMatrix = executor.submit(readMatrix, r.raw, unit, binsize, blockMap)
        else:
            req = open(self.infile, 'rb')
            req.seek(myFilePos)
            futureMatrix = executor.submit(readMatrix, req, unit, binsize, blockMap)

        if norm != "NONE":
            c1Norm = futureNorm1.result()
            if isIntra:
                c2Norm = c1Norm
            else:
                c2Norm = futureNorm2.result()
        else:
            c1Norm, c2Norm = None, None

        blockBinCount, blockColumnCount = futureMatrix.result()
        return normalizedmatrix(self.infile, self.is_synapse, binsize, isIntra, neededToFlipIndices, blockBinCount,
                                blockColumnCount, blockMap, norm, c1Norm, c2Norm, self.version)


class normalizedmatrix:
    def __init__(self, infile, is_synapse, binsize, isIntra, neededToFlipIndices, blockBinCount, blockColumnCount,
                 blockMap, norm, c1Norm, c2Norm, version):
        self.infile = infile
        self.is_synapse = is_synapse
        self.isHttpFile = infile.startswith("http")
        self.binsize = binsize
        self.isIntra = isIntra
        self.neededToFlipIndices = neededToFlipIndices
        self.blockBinCount = blockBinCount
        self.blockColumnCount = blockColumnCount
        self.norm = norm
        self.c1Norm = c1Norm
        self.c2Norm = c2Norm
        self.blockMap = blockMap
        self.version = version

    def getDataFromBinRegion(self, X1, X2, Y1, Y2):
        binsize = self.binsize
        if self.neededToFlipIndices:
            X1, X2, Y1, Y2 = Y1, Y2, X1, X2
        binPositionsBox = []
        binPositionsBox.append(int(X1))
        binPositionsBox.append(int(X2))
        binPositionsBox.append(int(Y1))
        binPositionsBox.append(int(Y2))

        blockNumbers = getBlockNumbersForRegionFromBinPosition(binPositionsBox, self.blockBinCount,
                                                               self.blockColumnCount, self.isIntra)
        yActual = []
        xActual = []
        counts = []

        executor = concurrent.futures.ProcessPoolExecutor()
        futures = [
            executor.submit(readBlockWorker, self.infile, self.is_synapse, bNum, binsize, self.blockMap, self.norm, \
                            self.c1Norm, self.c2Norm, binPositionsBox, self.isIntra, self.version) for bNum in
            blockNumbers]

        for future in futures:
            xTemp, yTemp, cTemp = future.result()
            xActual.extend(xTemp)
            yActual.extend(yTemp)
            counts.extend(cTemp)
        return [xActual, yActual, counts]

    def getDataFromGenomeRegion(self, X1, X2, Y1, Y2):
        binsize = self.binsize
        return self.getDataFromBinRegion(X1 / binsize, math.ceil(X2 / binsize), Y1 / binsize, math.ceil(Y2 / binsize))

    def getBatchedDataFromGenomeRegion(self, listOfCoordinates):
        executor = concurrent.futures.ThreadPoolExecutor()
        futures = [executor.submit(self.getDataFromGenomeRegion, a, b, c, d) for (a, b, c, d) in listOfCoordinates]
        finalResults = list()
        for future in futures:
            finalResults.append(future.result())
        return finalResults


def printme(norm, infile, chr1loc, chr2loc, unit, binsize, outfile):
    """ Reads a .hic file and extracts and prints the given contact matrix
    to a text file

    Args:
       norm(str): Normalization type, one of VC, KR, VC_SQRT, or NONE
       infile(str): File name or URL of .hic file
       chr1loc(str): Chromosome name and (optionally) range, i.e. "1" or "1:10000:25000"
       chr2loc(str): Chromosome name and (optionally) range, i.e. "1" or "1:10000:25000"
       unit(str): One of BP or FRAG
       binsize(int): Resolution, i.e. 25000 for 25K
       outfile(str): Name of text file to write to
    """
    f = open(outfile, 'w')
    strawObj = straw(infile)

    chr1, X1, X2 = strawObj.chromDotSizes.figureOutEndpoints(chr1loc)
    chr2, Y1, Y2 = strawObj.chromDotSizes.figureOutEndpoints(chr2loc)

    matrxObj = strawObj.getNormalizedMatrix(chr1, chr2, norm, unit, binsize)

    result = matrxObj.getDataFromGenomeRegion(X1, X2, Y1, Y2)

    for i in range(len(result[0])):
        f.write("{0}\t{1}\t{2}\n".format(result[0][i], result[1][i], result[2][i]))
    f.close()