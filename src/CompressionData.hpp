
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

#ifndef OPTIMIZEDDATA_HPP
#define OPTIMIZEDDATA_HPP

#include <cmath>
#include <cstdint>
#include "zfp.h"

/*!
  The number of bits per value to compress at.
*/
const int bitRate = 32;

/*!
  The number of bytes per block to compress at.
*/
const int blockRate = ceil(bitRate*4/8.0);

struct CompressionData_STRUCT {
  zfp_stream * zfp;
  bitstream * stream;
  uint8_t * buffer;
  local_int_t bufferSize;

  local_int_t * rowStarts;
  local_int_t * diagonalIndices;
};
typedef struct CompressionData_STRUCT CompressionData;

inline void DeleteCompressionData(CompressionData & data) {
  if (data.zip) {
    zfp_stream_close(data.zfp);
    data.zfp = 0;
  }
  if (data.stream) {
    stream_close(data.stream);
    data.stream = 0
  }
  if (data.buffer) {
    delete [] data.buffer;
    data.buffer = 0;
    bufferSize = 0;
  }
  if (data.diagonalIndices) {
    delete [] data.diagonalIndices;
    data.diagonalIndices = 0;
  }
  if (data.rowStarts) {
    delete [] data.rowStarts;
    data.rowStarts = 0;
  }
}


#endif //OPTIMIZEDDATA_HPP
