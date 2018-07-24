
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

#ifndef DECODEBLOCK_HPP
#define DECODEBLOCK_HPP


#include <cassert>
#include <cstring>
#include "CompressionData.hpp"
#include "SparseMatrix.hpp"

/*!
  Decompresses the requested index

  @param [in] data                 The CompressionData being used
  @param [in] row                  The current local matrix row
  @param [inout] value             Takes the previous value and becomes the next value
  @param [inout] pos               The current position in the stream
*/
inline void DecodeIndex(const CompressionData & data, const local_int_t row, local_int_t & value, uint64_t & pos) {

  uint8_t * arr = data.compressed[row];

  if (pos == 0) {
    std::memcpy(&value, arr, sizeof(local_int_t));
    value &= data.firstIndMask;
    pos += data.firstIndBits;
  } else {
    huffmanResult_t * huffmanTable = data.tables[0];

    while (true) {
      uint64_t val = peek(huffmanTable, arr, pos);
      if (val & RESULT_DONE_MASK) {
        pos += (val & RESULT_SIZE_MASK) >> RESULT_SIZE_OFFSET;
        value += (val & RESULT_VALUE_MASK);
        break;
      } else {
        pos += WINDOW_SIZE;
        huffmanTable = data.tables[val & RESULT_VALUE_MASK];
      }
    }
  }

}


#endif // DECODEBLOCK_HPP
