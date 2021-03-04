/*
  The MIT License (MIT)

  Copyright (c) 2011-2016 Broad Institute, Aiden Lab

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
#include <cstring>
#include <iostream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <cmath>
#include <set>
#include <vector>
#include <streambuf>
#include <curl/curl.h>
#include "zlib.h"
#include "straw.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
using namespace std;

/*
  Straw: fast C++ implementation of dump. Not as fully featured as the
  Java version. Reads the .hic file, finds the appropriate matrix and slice
  of data, and outputs as text in sparse upper triangular format.

  Currently only supporting matrices.

  Usage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>
 */
// this is for creating a stream from a byte array for ease of use
struct membuf : std::streambuf
{
  membuf(char* begin, char* end) {
    this->setg(begin, begin, end);
  }
};

// for holding data from URL call
struct MemoryStruct {
  char *memory;
  size_t size;
};

// version number
int version;

// map of block numbers to pointers


long total_bytes;

size_t hdf(char* b, size_t size, size_t nitems, void *userdata) {
  size_t numbytes = size * nitems;
  b[numbytes+1]='\0';
  string s(b);
  int found = s.find("Content-Range");
  if (found !=  string::npos) {
    int found2 = s.find("/");
    //Content-Range: bytes 0-100000/891471462
    if (found2 != string::npos) {
      string total=s.substr(found2+1);
      total_bytes = stol(total);
    }
  }

  return numbytes;
}

// callback for libcurl. data written to this buffer
static size_t
WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp)
{
  size_t realsize = size * nmemb;
  struct MemoryStruct *mem = (struct MemoryStruct *)userp;

  mem->memory = static_cast<char*>(realloc(mem->memory, mem->size + realsize + 1));
  if(mem->memory == NULL) {
    /* out of memory! */
    printf("not enough memory (realloc returned NULL)\n");
    return 0;
  }

  std::memcpy(&(mem->memory[mem->size]), contents, realsize);
  mem->size += realsize;
  mem->memory[mem->size] = 0;

  return realsize;
}

// get a buffer that can be used as an input stream from the URL
char* getData(CURL *curl, long position, int chunksize) {
  std::ostringstream oss;
  struct MemoryStruct chunk;

  chunk.memory = static_cast<char*>(malloc(1));
  chunk.size = 0;    /* no data at this point */
  oss << position << "-" << position + chunksize;
  curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *)&chunk);
  curl_easy_setopt(curl, CURLOPT_RANGE, oss.str().c_str());
  CURLcode res = curl_easy_perform(curl);
  if (res != CURLE_OK) {
    fprintf(stderr, "curl_easy_perform() failed: %s\n",
	    curl_easy_strerror(res));
  }
  //  printf("%lu bytes retrieved\n", (long)chunk.size);

  return chunk.memory;
}

// initialize the CURL stream
CURL* initCURL(const char* url) {
  CURL* curl = curl_easy_init();
  if(curl) {
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteMemoryCallback);
    curl_easy_setopt(curl, CURLOPT_URL, url);
    //curl_easy_setopt (curl, CURLOPT_VERBOSE, 1L);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, hdf);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "straw");
  }
  return curl;
}

// returns whether or not this is valid HiC file
bool readMagicString(istream& fin) {
  string str;
  getline(fin, str, '\0' );
  return str[0]=='H' && str[1]=='I' && str[2]=='C';
}

// reads the header, storing the positions of the normalization vectors and returning the master pointer
map<string, chromosome> readHeader(istream& fin, long& master) {
  map<string, chromosome> chromosomeMap;
  if (!readMagicString(fin)) {
    cerr << "Hi-C magic string is missing, does not appear to be a hic file" << endl;
    master=-1;
    return chromosomeMap;
  }

  fin.read((char*)&version, sizeof(int));
  if (version < 6) {
    cerr << "Version " << version << " no longer supported" << endl;
    master=-1;
    return chromosomeMap;
  }
  fin.read((char*)&master, sizeof(long));
  string genome;
  getline(fin, genome, '\0' );
  int nattributes;
  fin.read((char*)&nattributes, sizeof(int));

  // reading and ignoring attribute-value dictionary
  for (int i=0; i<nattributes; i++) {
    string key, value;
    getline(fin, key, '\0');
    getline(fin, value, '\0');
  }
  int nChrs;
  fin.read((char*)&nChrs, sizeof(int));
  // chromosome map for finding matrix
  for (int i=0; i<nChrs; i++) {
    string name;
    int length;
    getline(fin, name, '\0');
    fin.read((char*)&length, sizeof(int));
    chromosome chr;
    chr.index=i;
    chr.name=name;
    chr.length=length;
    chromosomeMap[name]=chr;
  }
  return chromosomeMap;
}

// reads the footer from the master pointer location. takes in the chromosomes,
// norm, unit (BP or FRAG) and resolution or binsize, and sets the file
// position of the matrix and the normalization vectors for those chromosomes
// at the given normalization and resolution
bool readFooter(istream& fin, long master, int c1, int c2, string norm, string unit, int resolution, long &myFilePos, indexEntry &c1NormEntry, indexEntry &c2NormEntry) {
  int nBytes;
  fin.read((char*)&nBytes, sizeof(int));
  stringstream ss;
  ss << c1 << "_" << c2;
  string key = ss.str();

  int nEntries;
  fin.read((char*)&nEntries, sizeof(int));
  bool found = false;
  for (int i=0; i<nEntries; i++) {
    string str;
    getline(fin, str, '\0');
    long fpos;
    fin.read((char*)&fpos, sizeof(long));
    int sizeinbytes;
    fin.read((char*)&sizeinbytes, sizeof(int));
    if (str == key) {
      myFilePos = fpos;
      found=true;
    }
  }
  if (!found) {
    cerr << "File doesn't have the given chr_chr map " << key << endl;
    return false;
  }

  if (norm=="NONE") return true; // no need to read norm vector index

  // read in and ignore expected value maps; don't store; reading these to
  // get to norm vector index
  int nExpectedValues;
  fin.read((char*)&nExpectedValues, sizeof(int));
  for (int i=0; i<nExpectedValues; i++) {
    string str;
    getline(fin, str, '\0'); //unit
    int binSize;
    fin.read((char*)&binSize, sizeof(int));

    int nValues;
    fin.read((char*)&nValues, sizeof(int));
    for (int j=0; j<nValues; j++) {
      double v;
      fin.read((char*)&v, sizeof(double));
    }

    int nNormalizationFactors;
    fin.read((char*)&nNormalizationFactors, sizeof(int));
    for (int j=0; j<nNormalizationFactors; j++) {
      int chrIdx;
      fin.read((char*)&chrIdx, sizeof(int));
      double v;
      fin.read((char*)&v, sizeof(double));
    }
  }
  fin.read((char*)&nExpectedValues, sizeof(int));
  for (int i=0; i<nExpectedValues; i++) {
    string str;
    getline(fin, str, '\0'); //typeString
    getline(fin, str, '\0'); //unit
    int binSize;
    fin.read((char*)&binSize, sizeof(int));

    int nValues;
    fin.read((char*)&nValues, sizeof(int));
    for (int j=0; j<nValues; j++) {
      double v;
      fin.read((char*)&v, sizeof(double));
    }
    int nNormalizationFactors;
    fin.read((char*)&nNormalizationFactors, sizeof(int));
    for (int j=0; j<nNormalizationFactors; j++) {
      int chrIdx;
      fin.read((char*)&chrIdx, sizeof(int));
      double v;
      fin.read((char*)&v, sizeof(double));
    }
  }
  // Index of normalization vectors
  fin.read((char*)&nEntries, sizeof(int));
  bool found1 = false;
  bool found2 = false;
  for (int i = 0; i < nEntries; i++) {
    string normtype;
    getline(fin, normtype, '\0'); //normalization type
    int chrIdx;
    fin.read((char*)&chrIdx, sizeof(int));
    string unit1;
    getline(fin, unit1, '\0'); //unit
    int resolution1;
    fin.read((char*)&resolution1, sizeof(int));
    long filePosition;
    fin.read((char*)&filePosition, sizeof(long));
    int sizeInBytes;
    fin.read((char*)&sizeInBytes, sizeof(int));
    if (chrIdx == c1 && normtype == norm && unit1 == unit && resolution1 == resolution) {
      c1NormEntry.position=filePosition;
      c1NormEntry.size=sizeInBytes;
      found1 = true;
    }
    if (chrIdx == c2 && normtype == norm && unit1 == unit && resolution1 == resolution) {
      c2NormEntry.position=filePosition;
      c2NormEntry.size=sizeInBytes;
      found2 = true;
    }
  }
  if (!found1 || !found2) {
    cerr << "File did not contain " << norm << " normalization vectors for one or both chromosomes at " << resolution << " " << unit << endl;
  }
  return true;
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count
map <int, indexEntry> readMatrixZoomData(istream& fin, string myunit, int mybinsize, int &myBlockBinCount, int &myBlockColumnCount, bool &found) {

  map <int, indexEntry> blockMap;
  string unit;
  getline(fin, unit, '\0' ); // unit
  int tmp;
  fin.read((char*)&tmp, sizeof(int)); // Old "zoom" index -- not used
  float tmp2;
  fin.read((char*)&tmp2, sizeof(float)); // sumCounts
  fin.read((char*)&tmp2, sizeof(float)); // occupiedCellCount
  fin.read((char*)&tmp2, sizeof(float)); // stdDev
  fin.read((char*)&tmp2, sizeof(float)); // percent95
  int binSize;
  fin.read((char*)&binSize, sizeof(int));
  int blockBinCount;
  fin.read((char*)&blockBinCount, sizeof(int));
  int blockColumnCount;
  fin.read((char*)&blockColumnCount, sizeof(int));

  found = false;
  if (myunit==unit && mybinsize==binSize) {
    myBlockBinCount = blockBinCount;
    myBlockColumnCount = blockColumnCount;
    found = true;
  }

  int nBlocks;
  fin.read((char*)&nBlocks, sizeof(int));

  for (int b = 0; b < nBlocks; b++) {
    int blockNumber;
    fin.read((char*)&blockNumber, sizeof(int));
    long filePosition;
    fin.read((char*)&filePosition, sizeof(long));
    int blockSizeInBytes;
    fin.read((char*)&blockSizeInBytes, sizeof(int));
    indexEntry entry;
    entry.size = blockSizeInBytes;
    entry.position = filePosition;
    if (found) blockMap[blockNumber] = entry;
  }
  return blockMap;
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count
map <int, indexEntry> readMatrixZoomDataHttp(CURL* curl, long &myFilePosition, string myunit, int mybinsize, int &myBlockBinCount, int &myBlockColumnCount, bool &found) {

  map <int, indexEntry> blockMap;
  char* buffer;
  int header_size = 5*sizeof(int)+4*sizeof(float);
  char* first;
  first = getData(curl, myFilePosition, 1);
  if (first[0]=='B') {
    header_size+=3;
  }
  else if (first[0]=='F') {
    header_size+=5;
  }
  else {
    cerr << "Unit not understood" << endl;
    return blockMap;
  }
  buffer = getData(curl, myFilePosition, header_size);
  membuf sbuf(buffer, buffer + header_size);
  istream fin(&sbuf);

  string unit;
  getline(fin, unit, '\0' ); // unit
  int tmp;
  fin.read((char*)&tmp, sizeof(int)); // Old "zoom" index -- not used
  float tmp2;
  fin.read((char*)&tmp2, sizeof(float)); // sumCounts
  fin.read((char*)&tmp2, sizeof(float)); // occupiedCellCount
  fin.read((char*)&tmp2, sizeof(float)); // stdDev
  fin.read((char*)&tmp2, sizeof(float)); // percent95
  int binSize;
  fin.read((char*)&binSize, sizeof(int));
  int blockBinCount;
  fin.read((char*)&blockBinCount, sizeof(int));
  int blockColumnCount;
  fin.read((char*)&blockColumnCount, sizeof(int));

  found = false;
  if (myunit==unit && mybinsize==binSize) {
    myBlockBinCount = blockBinCount;
    myBlockColumnCount = blockColumnCount;
    found = true;
  }

  int nBlocks;
  fin.read((char*)&nBlocks, sizeof(int));

  if (found) {
    buffer = getData(curl, myFilePosition+header_size, nBlocks*(sizeof(int)+sizeof(long)+sizeof(int)));
    membuf sbuf2(buffer, buffer + nBlocks*(sizeof(int)+sizeof(long)+sizeof(int)));
    istream fin2(&sbuf2);
    for (int b = 0; b < nBlocks; b++) {
      int blockNumber;
      fin2.read((char*)&blockNumber, sizeof(int));
      long filePosition;
      fin2.read((char*)&filePosition, sizeof(long));
      int blockSizeInBytes;
      fin2.read((char*)&blockSizeInBytes, sizeof(int));
      indexEntry entry;
      entry.size = blockSizeInBytes;
      entry.position = filePosition;
      blockMap[blockNumber] = entry;
    }
  }
  else {
    myFilePosition = myFilePosition+header_size+(nBlocks*(sizeof(int)+sizeof(long)+sizeof(int)));
  }
  delete buffer;
  return blockMap;
}

// goes to the specified file pointer in http and finds the raw contact matrix at specified resolution, calling readMatrixZoomData.
// sets blockbincount and blockcolumncount
map <int, indexEntry> readMatrixHttp(CURL *curl, long myFilePosition, string unit, int resolution, int &myBlockBinCount, int &myBlockColumnCount) {
  char * buffer;
  int size = sizeof(int)*3;
  buffer = getData(curl, myFilePosition, size);
  membuf sbuf(buffer, buffer + size);
  istream bufin(&sbuf);

  int c1,c2;
  bufin.read((char*)&c1, sizeof(int)); //chr1
  bufin.read((char*)&c2, sizeof(int)); //chr2
  int nRes;
  bufin.read((char*)&nRes, sizeof(int));
  int i=0;
  bool found=false;
  myFilePosition=myFilePosition+size;
  delete buffer;
  map <int, indexEntry> blockMap;

  while (i<nRes && !found) {
    // myFilePosition gets updated within call
    blockMap = readMatrixZoomDataHttp(curl, myFilePosition, unit, resolution, myBlockBinCount, myBlockColumnCount, found);
    i++;
  }
  if (!found) {
    cerr << "Error finding block data" << endl;
  }
  return blockMap;
}

// goes to the specified file pointer and finds the raw contact matrix at specified resolution, calling readMatrixZoomData.
// sets blockbincount and blockcolumncount
map <int, indexEntry> readMatrix(istream& fin, long myFilePosition, string unit, int resolution, int &myBlockBinCount, int &myBlockColumnCount) {
  map <int, indexEntry> blockMap;

  fin.seekg(myFilePosition, ios::beg);
  int c1,c2;
  fin.read((char*)&c1, sizeof(int)); //chr1
  fin.read((char*)&c2, sizeof(int)); //chr2
  int nRes;
  fin.read((char*)&nRes, sizeof(int));
  int i=0;
  bool found=false;
  while (i<nRes && !found) {
    blockMap = readMatrixZoomData(fin, unit, resolution, myBlockBinCount, myBlockColumnCount, found);
    i++;
  }
  if (!found) {
    cerr << "Error finding block data" << endl;
  }
  return blockMap;
}
// gets the blocks that need to be read for this slice of the data.  needs blockbincount, blockcolumncount, and whether
// or not this is intrachromosomal.
set<int> getBlockNumbersForRegionFromBinPosition(int* regionIndices, int blockBinCount, int blockColumnCount, bool intra) {
   int col1 = regionIndices[0] / blockBinCount;
   int col2 = (regionIndices[1] + 1) / blockBinCount;
   int row1 = regionIndices[2] / blockBinCount;
   int row2 = (regionIndices[3] + 1) / blockBinCount;

   set<int> blocksSet;
   // first check the upper triangular matrix
   for (int r = row1; r <= row2; r++) {
     for (int c = col1; c <= col2; c++) {
       int blockNumber = r * blockColumnCount + c;
       blocksSet.insert(blockNumber);
     }
   }
   // check region part that overlaps with lower left triangle
   // but only if intrachromosomal
   if (intra) {
     for (int r = col1; r <= col2; r++) {
       for (int c = row1; c <= row2; c++) {
	 int blockNumber = r * blockColumnCount + c;
	 blocksSet.insert(blockNumber);
       }
     }
   }

   return blocksSet;
}

// this is the meat of reading the data.  takes in the block number and returns the set of contact records corresponding to
// that block.  the block data is compressed and must be decompressed using the zlib library functions
vector<contactRecord> readBlock(istream& fin, CURL* curl, bool isHttp, indexEntry idx) {
  if (idx.size == 0) {
    vector<contactRecord> v;
    return v;
  }
  char* compressedBytes = new char[idx.size];
  char* uncompressedBytes = new char[idx.size*10]; //biggest seen so far is 3

  if (isHttp) {
    compressedBytes = getData(curl, idx.position, idx.size);
  }
  else {
    fin.seekg(idx.position, ios::beg);
    fin.read(compressedBytes, idx.size);
  }
  // Decompress the block
  // zlib struct
  z_stream infstream;
  infstream.zalloc = Z_NULL;
  infstream.zfree = Z_NULL;
  infstream.opaque = Z_NULL;
  infstream.avail_in = (uInt)(idx.size); // size of input
  infstream.next_in = (Bytef *)compressedBytes; // input char array
  infstream.avail_out = (uInt)idx.size*10; // size of output
  infstream.next_out = (Bytef *)uncompressedBytes; // output char array
  // the actual decompression work.
  inflateInit(&infstream);
  inflate(&infstream, Z_NO_FLUSH);
  inflateEnd(&infstream);
  int uncompressedSize=infstream.total_out;

  // create stream from buffer for ease of use
  membuf sbuf(uncompressedBytes, uncompressedBytes + uncompressedSize);
  istream bufferin(&sbuf);
  int nRecords;
  bufferin.read((char*)&nRecords, sizeof(int));
  vector<contactRecord> v(nRecords);
  // different versions have different specific formats
  if (version < 7) {
    for (int i = 0; i < nRecords; i++) {
      int binX, binY;
      bufferin.read((char*)&binX, sizeof(int));
      bufferin.read((char*)&binY, sizeof(int));
      float counts;
      bufferin.read((char*)&counts, sizeof(float));
      contactRecord record;
      record.binX = binX;
      record.binY = binY;
      record.counts = counts;
      v[i] = record;
    }
  }
  else {
    int binXOffset, binYOffset;
    bufferin.read((char*)&binXOffset, sizeof(int));
    bufferin.read((char*)&binYOffset, sizeof(int));
    char useShort;
    bufferin.read((char*)&useShort, sizeof(char));
    char type;
    bufferin.read((char*)&type, sizeof(char));
    int index=0;
    if (type == 1) {
      // List-of-rows representation
      short rowCount;
      bufferin.read((char*)&rowCount, sizeof(short));
      for (int i = 0; i < rowCount; i++) {
	short y;
	bufferin.read((char*)&y, sizeof(short));
	int binY = y + binYOffset;
	short colCount;
	bufferin.read((char*)&colCount, sizeof(short));
	for (int j = 0; j < colCount; j++) {
	  short x;
	  bufferin.read((char*)&x, sizeof(short));
	  int binX = binXOffset + x;
	  float counts;
	  if (useShort == 0) { // yes this is opposite of usual
	    short c;
	    bufferin.read((char*)&c, sizeof(short));
	    counts = c;
	  }
	  else {
	    bufferin.read((char*)&counts, sizeof(float));
	  }
	  contactRecord record;
	  record.binX = binX;
	  record.binY = binY;
	  record.counts = counts;
	  v[index]=record;
	  index++;
	}
      }
    }
    else if (type == 2) {
      int nPts;
      bufferin.read((char*)&nPts, sizeof(int));
      short w;
      bufferin.read((char*)&w, sizeof(short));

      for (int i = 0; i < nPts; i++) {
	//int idx = (p.y - binOffset2) * w + (p.x - binOffset1);
	int row = i / w;
	int col = i - row * w;
	int bin1 = binXOffset + col;
	int bin2 = binYOffset + row;

	float counts;
	if (useShort == 0) { // yes this is opposite of the usual
	  short c;
	  bufferin.read((char*)&c, sizeof(short));
	  if (c != -32768) {
	    contactRecord record;
	    record.binX = bin1;
	    record.binY = bin2;
	    record.counts = c;
	    v[index]=record;
	    index++;
	  }
	}
	else {
	  bufferin.read((char*)&counts, sizeof(float));
	  if (!isnan(counts)) {
	    contactRecord record;
	    record.binX = bin1;
	    record.binY = bin2;
	    record.counts = counts;
	    v[index]=record;
	    index++;
	  }
	}
      }
    }
  }
  delete[] compressedBytes;
  delete[] uncompressedBytes; // don't forget to delete your heap arrays in C++!
  return v;
}

int readSize(istream& fin, CURL* curl, bool isHttp, indexEntry idx) {
  if (idx.size == 0) {
    return 0;
  }
  char* compressedBytes = new char[idx.size];
  char* uncompressedBytes = new char[idx.size*10];

  if (isHttp) {
    compressedBytes = getData(curl, idx.position, idx.size);
  }
  else {
    fin.seekg(idx.position, ios::beg);
    fin.read(compressedBytes, idx.size);
  }
  // Decompress the block
  // zlib struct
  z_stream infstream;
  infstream.zalloc = Z_NULL;
  infstream.zfree = Z_NULL;
  infstream.opaque = Z_NULL;
  infstream.avail_in = (uInt)(idx.size); // size of input
  infstream.next_in = (Bytef *)compressedBytes; // input char array
  infstream.avail_out = (uInt)idx.size*10; // size of output
  infstream.next_out = (Bytef *)uncompressedBytes; // output char array
  // the actual decompression work.
  inflateInit(&infstream);
  inflate(&infstream, Z_NO_FLUSH);
  inflateEnd(&infstream);
  int uncompressedSize=infstream.total_out;

  // create stream from buffer for ease of use
  membuf sbuf(uncompressedBytes, uncompressedBytes + uncompressedSize);
  istream bufferin(&sbuf);
  int nRecords;
  bufferin.read((char*)&nRecords, sizeof(int));
  delete[] compressedBytes;
  delete[] uncompressedBytes;
  return nRecords;
}


// reads the normalization vector from the file at the specified location
vector<double> readNormalizationVector(istream& bufferin) {
  int nValues;
  bufferin.read((char*)&nValues, sizeof(int));
  vector<double> values(nValues);
  //  bool allNaN = true;

  for (int i = 0; i < nValues; i++) {
    double d;
    bufferin.read((char*)&d, sizeof(double));
    values[i] = d;
    /* if (!Double.isNaN(values[i])) {
      allNaN = false;
      }*/
  }
  //  if (allNaN) return null;
  return values;
}

vector<contactRecord> straw(string norm, string fname, string chr1loc, string chr2loc, string unit, int binsize)
{
  if (!(unit=="BP"||unit=="FRAG")) {
    cerr << "Norm specified incorrectly, must be one of <BP/FRAG>" << endl;
    cerr << "Usage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>" << endl;
    vector<contactRecord> v;
    return v;
  }

  // HTTP code
  string prefix="http";
  bool isHttp = false;
  ifstream fin;

  // read header into buffer; 100K should be sufficient
  CURL *curl;

  long master;
  map <string, chromosome> chromosomeMap;

  if (std::strncmp(fname.c_str(), prefix.c_str(), prefix.size()) == 0) {
    isHttp = true;
    char * buffer;
    curl = initCURL(fname.c_str());
    if (curl) {
      buffer = getData(curl, 0, 100000);
    }
    else {
      cerr << "URL " << fname << " cannot be opened for reading" << endl;
      vector<contactRecord> v;
      return v;
    }
    membuf sbuf(buffer, buffer + 100000);
    istream bufin(&sbuf);
    chromosomeMap = readHeader(bufin, master);
    delete buffer;
  }
  else {
    fin.open(fname, fstream::in);
    if (!fin) {
      cerr << "File " << fname << " cannot be opened for reading" << endl;
      vector<contactRecord> v;
      return v;
    }
    chromosomeMap = readHeader(fin, master);
  }

  // parse chromosome positions
  stringstream ss(chr1loc);
  string chr1, chr2, x, y;
  int c1pos1=-100, c1pos2=-100, c2pos1=-100, c2pos2=-100;
  getline(ss, chr1, ':');
  if (chromosomeMap.count(chr1) == 0) {
    cerr << chr1 << " not found in the file." << endl;
    vector<contactRecord> v;
    return v;
  }

  if (getline(ss, x, ':') && getline(ss, y, ':')) {
    c1pos1 = stoi(x);
    c1pos2 = stoi(y);
  }
  else {
    c1pos1 = 0;
    c1pos2 = chromosomeMap[chr1].length;
  }
  stringstream ss1(chr2loc);
  getline(ss1, chr2, ':');
  if (chromosomeMap.count(chr2) == 0) {
    cerr << chr2 << " not found in the file." << endl;
    vector<contactRecord> v;
    return v;
  }

  if (getline(ss1, x, ':') && getline(ss1, y, ':')) {
    c2pos1 = stoi(x);
    c2pos2 = stoi(y);
  }
  else {
    c2pos1 = 0;
    c2pos2 = chromosomeMap[chr2].length;
  }

  // from header have size of chromosomes, set region to read
  int c1=min(chromosomeMap[chr1].index,chromosomeMap[chr2].index);
  int c2=max(chromosomeMap[chr1].index,chromosomeMap[chr2].index);
  int origRegionIndices[4]; // as given by user
  int regionIndices[4]; // used to find the blocks we need to access
  // reverse order if necessary
  if (chromosomeMap[chr1].index > chromosomeMap[chr2].index) {
    origRegionIndices[0] = c2pos1;
    origRegionIndices[1] = c2pos2;
    origRegionIndices[2] = c1pos1;
    origRegionIndices[3] = c1pos2;
    regionIndices[0] = c2pos1 / binsize;
    regionIndices[1] = c2pos2 / binsize;
    regionIndices[2] = c1pos1 / binsize;
    regionIndices[3] = c1pos2 / binsize;
  }
  else {
    origRegionIndices[0] = c1pos1;
    origRegionIndices[1] = c1pos2;
    origRegionIndices[2] = c2pos1;
    origRegionIndices[3] = c2pos2;
    regionIndices[0] = c1pos1 / binsize;
    regionIndices[1] = c1pos2 / binsize;
    regionIndices[2] = c2pos1 / binsize;
    regionIndices[3] = c2pos2 / binsize;
  }

  indexEntry c1NormEntry, c2NormEntry;
  long myFilePos;

  long bytes_to_read = total_bytes - master;
  bool foundFooter = false;
  if (isHttp) {
    char* buffer2;
    buffer2 = getData(curl, master, bytes_to_read);
    membuf sbuf2(buffer2, buffer2 + bytes_to_read);
    istream bufin2(&sbuf2);
    foundFooter = readFooter(bufin2, master, c1, c2, norm, unit, binsize, myFilePos, c1NormEntry, c2NormEntry);
    delete buffer2;
  }
  else {
    fin.seekg(master, ios::beg);
    foundFooter = readFooter(fin, master, c1, c2, norm, unit, binsize, myFilePos, c1NormEntry, c2NormEntry);
  }
  // readFooter will assign the above variables

  if (!foundFooter) {
    vector<contactRecord> v;
    return v;
  }

  vector<double> c1Norm;
  vector<double> c2Norm;

  if (norm != "NONE") {
    char* buffer3;
    if (isHttp) {
      buffer3 = getData(curl, c1NormEntry.position, c1NormEntry.size);
    }
    else {
      buffer3 = new char[c1NormEntry.size];
      fin.seekg(c1NormEntry.position, ios::beg);
      fin.read(buffer3, c1NormEntry.size);
    }
    membuf sbuf3(buffer3, buffer3 + c1NormEntry.size);
    istream bufferin(&sbuf3);
    c1Norm = readNormalizationVector(bufferin);

    char* buffer4;
    if (isHttp) {
      buffer4 = getData(curl, c2NormEntry.position, c2NormEntry.size);
    }
    else {
      buffer4 = new char[c2NormEntry.size];
      fin.seekg(c2NormEntry.position, ios::beg);
      fin.read(buffer4, c2NormEntry.size);
    }
    membuf sbuf4(buffer4, buffer4 + c2NormEntry.size);
    istream bufferin2(&sbuf4);
    c2Norm = readNormalizationVector(bufferin2);
    delete buffer3;
    delete buffer4;
  }

  int blockBinCount, blockColumnCount;
  map<int, indexEntry> blockMap;

  if (isHttp) {
    // readMatrix will assign blockBinCount and blockColumnCount
    blockMap = readMatrixHttp(curl, myFilePos, unit, binsize, blockBinCount, blockColumnCount);
  }
  else {
    // readMatrix will assign blockBinCount and blockColumnCount
    blockMap = readMatrix(fin, myFilePos, unit, binsize, blockBinCount, blockColumnCount);
  }

  set<int> blockNumbers = getBlockNumbersForRegionFromBinPosition(regionIndices, blockBinCount, blockColumnCount, c1==c2);

  // getBlockIndices
  vector<contactRecord> records;
  vector<contactRecord> tmp_records;
  for (set<int>::iterator it=blockNumbers.begin(); it!=blockNumbers.end(); ++it) {
    // get contacts in this block
    tmp_records = readBlock(fin, curl, isHttp, blockMap[*it]);
    for (vector<contactRecord>::iterator it2=tmp_records.begin(); it2!=tmp_records.end(); ++it2) {
      contactRecord rec = *it2;

      int x = rec.binX * binsize;
      int y = rec.binY * binsize;
      float c = rec.counts;
      if (norm != "NONE") {
	c = c / (c1Norm[rec.binX] * c2Norm[rec.binY]);
      }

      if ((x >= origRegionIndices[0] && x <= origRegionIndices[1] &&
	   y >= origRegionIndices[2] && y <= origRegionIndices[3]) ||
	  // or check regions that overlap with lower left
	  ((c1==c2) && y >= origRegionIndices[0] && y <= origRegionIndices[1] && x >= origRegionIndices[2] && x <= origRegionIndices[3])) {
	contactRecord record;
	record.binX = x;
	record.binY = y;
	record.counts = c;
	records.push_back(record);
	//printf("%d\t%d\t%.14g\n", x, y, c);
      }
    }
  }
      //      free(chunk.memory);
      /* always cleanup */
      // curl_easy_cleanup(curl);
      //    curl_global_cleanup();
  return records;
}


int getSize(string norm, string fname, string chr1loc, string chr2loc, string unit, int binsize)
{
  if (!(unit=="BP"||unit=="FRAG")) {
    cerr << "Norm specified incorrectly, must be one of <BP/FRAG>" << endl;
    cerr << "Usage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>" << endl;
    return 0;
  }

  // HTTP code
  string prefix="http";
  bool isHttp = false;
  ifstream fin;

  // read header into buffer; 100K should be sufficient
  CURL *curl;

  long master;
  map<string, chromosome> chromosomeMap;
  if (std::strncmp(fname.c_str(), prefix.c_str(), prefix.size()) == 0) {
    isHttp = true;
    char * buffer;
    curl = initCURL(fname.c_str());
    if (curl) {
      buffer = getData(curl, 0, 100000);
    }
    else {
      cerr << "URL " << fname << " cannot be opened for reading" << endl;
      return 0;
    }
    membuf sbuf(buffer, buffer + 100000);
    istream bufin(&sbuf);
    chromosomeMap = readHeader(bufin, master);
    delete buffer;
  }
  else {
    fin.open(fname, fstream::in);
    if (!fin) {
      cerr << "File " << fname << " cannot be opened for reading" << endl;
      return 0;
    }
    chromosomeMap = readHeader(fin, master);
  }

  // parse chromosome positions
  stringstream ss(chr1loc);
  string chr1, chr2, x, y;
  int c1pos1=-100, c1pos2=-100, c2pos1=-100, c2pos2=-100;
  getline(ss, chr1, ':');
  if (chromosomeMap.count(chr1) == 0) {
    cerr << chr1 << " not found in the file." << endl;
    return 0;
  }

  if (getline(ss, x, ':') && getline(ss, y, ':')) {
    c1pos1 = stoi(x);
    c1pos2 = stoi(y);
  }
  else {
    c1pos1 = 0;
    c1pos2 = chromosomeMap[chr1].length;
  }
  stringstream ss1(chr2loc);
  getline(ss1, chr2, ':');
  if (chromosomeMap.count(chr2) == 0) {
    cerr << chr2 << " not found in the file." << endl;
    return 0;
  }
  if (getline(ss1, x, ':') && getline(ss1, y, ':')) {
    c2pos1 = stoi(x);
    c2pos2 = stoi(y);
  }
  else {
    c2pos1 = 0;
    c2pos2 = chromosomeMap[chr2].length;
  }

  // from header have size of chromosomes, set region to read
  int c1=min(chromosomeMap[chr1].index,chromosomeMap[chr2].index);
  int c2=max(chromosomeMap[chr1].index,chromosomeMap[chr2].index);
  int origRegionIndices[4]; // as given by user
  int regionIndices[4]; // used to find the blocks we need to access
  // reverse order if necessary
  if (chromosomeMap[chr1].index > chromosomeMap[chr2].index) {
    origRegionIndices[0] = c2pos1;
    origRegionIndices[1] = c2pos2;
    origRegionIndices[2] = c1pos1;
    origRegionIndices[3] = c1pos2;
    regionIndices[0] = c2pos1 / binsize;
    regionIndices[1] = c2pos2 / binsize;
    regionIndices[2] = c1pos1 / binsize;
    regionIndices[3] = c1pos2 / binsize;
  }
  else {
    origRegionIndices[0] = c1pos1;
    origRegionIndices[1] = c1pos2;
    origRegionIndices[2] = c2pos1;
    origRegionIndices[3] = c2pos2;
    regionIndices[0] = c1pos1 / binsize;
    regionIndices[1] = c1pos2 / binsize;
    regionIndices[2] = c2pos1 / binsize;
    regionIndices[3] = c2pos2 / binsize;
  }

  indexEntry c1NormEntry, c2NormEntry;
  long myFilePos;

  long bytes_to_read = total_bytes - master;
  bool foundFooter = false;

  if (isHttp) {
    char* buffer2;
    buffer2 = getData(curl, master, bytes_to_read);
    membuf sbuf2(buffer2, buffer2 + bytes_to_read);
    istream bufin2(&sbuf2);
    foundFooter = readFooter(bufin2, master, c1, c2, norm, unit, binsize, myFilePos, c1NormEntry, c2NormEntry);
    delete buffer2;
  }
  else {
    fin.seekg(master, ios::beg);
    foundFooter = readFooter(fin, master, c1, c2, norm, unit, binsize, myFilePos, c1NormEntry, c2NormEntry);
  }
  // readFooter will assign the above variables

  if (!foundFooter) return 0;

  vector<double> c1Norm;
  vector<double> c2Norm;

  if (norm != "NONE") {
    char* buffer3;
    if (isHttp) {
      buffer3 = getData(curl, c1NormEntry.position, c1NormEntry.size);
    }
    else {
      buffer3 = new char[c1NormEntry.size];
      fin.seekg(c1NormEntry.position, ios::beg);
      fin.read(buffer3, c1NormEntry.size);
    }
    membuf sbuf3(buffer3, buffer3 + c1NormEntry.size);
    istream bufferin(&sbuf3);
    c1Norm = readNormalizationVector(bufferin);

    char* buffer4;
    if (isHttp) {
      buffer4 = getData(curl, c2NormEntry.position, c2NormEntry.size);
    }
    else {
      buffer4 = new char[c2NormEntry.size];
      fin.seekg(c2NormEntry.position, ios::beg);
      fin.read(buffer4, c2NormEntry.size);
    }
    membuf sbuf4(buffer4, buffer4 + c2NormEntry.size);
    istream bufferin2(&sbuf4);
    c2Norm = readNormalizationVector(bufferin2);
    delete buffer3;
    delete buffer4;
  }

  int blockBinCount, blockColumnCount;
  map<int, indexEntry> blockMap;
  if (isHttp) {
    // readMatrix will assign blockBinCount and blockColumnCount
    blockMap = readMatrixHttp(curl, myFilePos, unit, binsize, blockBinCount, blockColumnCount);
  }
  else {
    // readMatrix will assign blockBinCount and blockColumnCount
    blockMap = readMatrix(fin, myFilePos, unit, binsize, blockBinCount, blockColumnCount);
  }
  set<int> blockNumbers = getBlockNumbersForRegionFromBinPosition(regionIndices, blockBinCount, blockColumnCount, c1==c2);

  // getBlockIndices
  vector<contactRecord> tmp_records;
  int count=0;
  for (set<int>::iterator it=blockNumbers.begin(); it!=blockNumbers.end(); ++it) {
    // get contacts in this block
    count += readSize(fin, curl, isHttp, blockMap[*it]);
  }
  return count;
}


namespace py = pybind11;

PYBIND11_MODULE(strawC, m) {
  m.doc() = R"pbdoc(
        New straw with pybind
        -----------------------

        .. currentmodule:: straw

        .. autosummary::
           :toctree: _generate

           straw
Straw enables programmatic access to .hic files.
.hic files store the contact matrices from Hi-C experiments and the
normalization and expected vectors, along with meta-data in the header.
The main function, straw, takes in the normalization, the filename or URL,
chromosome1 (and optional range), chromosome2 (and optional range),
whether the bins desired are fragment or base pair delimited, and bin size.
It then reads the header, follows the various pointers to the desired matrix
and normalization vector, and stores as [x, y, count]
Usage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>

Example:
>>>import strawC
>>>result = strawC.strawC('NONE', 'HIC001.hic', 'X', 'X', 'BP', 1000000)
>>>for i in range(len(result)):
...   print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))
See https://github.com/theaidenlab/straw/wiki/Python for more documentation
    )pbdoc";

  m.def("strawC", &straw, R"pbdoc(
        Straw: fast C++ implementation of dump.

        Bound with pybind
Usage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>
    )pbdoc");

  py::class_<contactRecord>(m, "contactRecord")
    .def(py::init<>())
    .def_readwrite("binX", &contactRecord::binX)
    .def_readwrite("binY", &contactRecord::binY)
    .def_readwrite("counts", &contactRecord::counts)
    ;

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}
