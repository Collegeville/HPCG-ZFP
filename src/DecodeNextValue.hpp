
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

#ifndef DECODENEXTVALUE_HPP
#define DECODENEXTVALUE_HPP

#include "SparseMatrix.hpp"

/*!
  Decompresses the next value of the matrix

  @param [in] mat                  The matrix being used
  @param [in] id                   The index to fetch
  @param [inout] uncompressedCount The number of previously uncompressed values, incremented as nessacery
  @param [inout] value             The previous value decoded, replaced with the new value
  @param [inout] previously        The previous value from the last call to this function
*/
inline void DecodeNextValue(const SparseMatrix & mat, local_int_t & id, local_int_t & uncompressedCount, double & value, double & previous, bool isForward) {
  const uint8_t * compressed;
  if (isForward) {
    compressed = ((CompressionData*)mat.optimizationData)->fValsCompressed;
  } else {
    compressed = ((CompressionData*)mat.optimizationData)->bValsCompressed;
  }

  switch ((compressed[id/VALUES_PER_COMPRESSED_BYTE] >> ((id%VALUES_PER_COMPRESSED_BYTE) * COMPRESSED_BITS)) & COMPRESSED_VALUE_MASK) {
    case NEIGHBOR:
      previous = value;
      break;
    case OVER_NEIGHBOR: {
      double temp = previous;
      previous = value;
      value = temp;
      break;
    }
    default: //UNCOMPRESSED
      const double * uncompressedArray;
      if (isForward) {
        uncompressedArray = ((CompressionData*)mat.optimizationData)->fValsUncompressed;
      } else {
        uncompressedArray = ((CompressionData*)mat.optimizationData)->bValsUncompressed;
      }
      previous = value;
      value = uncompressedArray[uncompressedCount];
      uncompressedCount++;
      break;
  }
  id++;
}

#endif // DECODENEXTVALUE_HPP
