
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
#include "CompressionData.hpp"
#include "SparseMatrix.hpp"

/*!
  Decompresses the next index of the matrix's foward iteration.

  @param [in] mat                  The matrix being used
  @param [in] id                   The index to fetch
  @param [inout] uncompressedCount The number of previously uncompressed values, incremented as nessacery
  @param [inout] value             The previous value decoded, replaced with the new value

  @return 0 on success, otherwise nonzero
*/
inline void DecodeNextIndex(const SparseMatrix & mat, local_int_t & id, local_int_t & uncompressedCount, local_int_t & value, bool isForward) {
  const uint8_t * compressed;
  if (isForward) {
    compressed = ((CompressionData*)mat.optimizationData)->fIndsCompressed;
  } else {
    compressed = ((CompressionData*)mat.optimizationData)->bIndsCompressed;
  }

  if (compressed[id/8] & (INCREMENT<<(id%8))) {
    //INCREMENT
    value++;
  } else {
    //UNCOMPRESSED
    const local_int_t * uncompressedArray;
    if (isForward) {
      uncompressedArray = ((CompressionData*)mat.optimizationData)->fIndsUncompressed;
    } else {
      uncompressedArray = ((CompressionData*)mat.optimizationData)->bIndsUncompressed;
    }
    value = uncompressedArray[uncompressedCount];
    uncompressedCount++;
  }
}


#endif // DECODEBLOCK_HPP
