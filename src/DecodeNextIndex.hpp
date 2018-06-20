
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

#include "Geometry.hpp"

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
inline void DecodeNextIndex_Forward(const SparseMatrix & mat, local_int_t & id, local_int_t & uncompressedCount, local_int_t & value) {
  const unsigned char * compressed = ((CompressionData*)mat.optimizationData)->forwardCompressed;

  if (compressed[id/8] & (INCREMENT<<(id%8))) {
    //INCREMENT
    value++;
  } else {
    //UNCOMPRESSED
    const local_int_t * uncompressedArray = ((CompressionData*)mat.optimizationData)->forwardUncompressed;
    value = uncompressedArray[uncompressedCount];
    uncompressedCount++;
  }
  id++;
}

/*!
  Decompresses the next index of the matrix's backwards iteration.

  @param [in] mat                  The matrix being used
  @param [in] id                   The index to fetch
  @param [inout] uncompressedCount The number of previously uncompressed values, incremented as nessacery
  @param [inout] value             The previous value decoded, replaced with the new value

  @return 0 on success, otherwise nonzero
*/
inline void DecodeNextIndex_Backward(const SparseMatrix & mat, local_int_t & id, local_int_t & uncompressedCount, local_int_t & value) {
  const unsigned char * compressed = ((CompressionData*)mat.optimizationData)->backwardCompressed;

  if (compressed[id/8] & (INCREMENT<<(id%8))) {
    //INCREMENT
    value--;
  } else {
    //UNCOMPRESSED
    const local_int_t * uncompressedArray =  ((CompressionData*)mat.optimizationData)->backwardUncompressed;
    value = uncompressedArray[uncompressedCount];
    uncompressedCount++;
  }
  id++;
}


#endif // DECODEBLOCK_HPP
