
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

#ifndef DECODENEXTINDEX_HPP
#define DECODENEXTINDEX_HPP

#include <cstdint>
#include "SparseMatrix.hpp"

inline uint32_t getBits(const uint8_t * buffer, int bitPos) {
  // based off https://stackoverflow.com/a/5723250/6353993
  uint64_t bits = *(uint64_t*)(buffer+bitPos/8);
  return (uint32_t)(bits >> bitPos%8);
}

// count number of zeros in the front of the read location
// Because this is only used to decode the Elias gamma encoded length of the primary value, there should be at most 5 zeros (log2(32)-1)
inline int countZeros(uint32_t bits) {
  //copied from https://graphics.stanford.edu/~seander/bithacks.html
  if (bits & 0x1) {
    return 0;
  } else {
    int count = 1;
    if ((bits & 0xF) == 0) {
      bits >>= 4;
      count += 4;
    } else if ((bits & 0x3) == 0) {
      bits >>= 2;
      count += 2;
    } else {
      //if less than 2 zeros, must be 1
      return count;
    }
    return count - (bits&0x1);
  }
}

/*!
  Decompresses the next index of the matrix's

  @param [in] mat                  The matrix being used
  @param [in] id                   The index to fetch
  @param [inout] uncompressedCount The number of previously uncompressed values, incremented as nessacery
  @param [inout] value             The previous value decoded, replaced with the new value

  @return 0 on success, otherwise nonzero
*/
inline void DecodeNextIndex(const uint8_t * buffer, int & bitPosition, local_int_t & previousIndex) {
  bool rowStart = bitPosition == 0;
  uint32_t bits = getBits(buffer, bitPosition);

  //length of N
  int numZeros = countZeros(bits);
  int explicitBitsMask = (1 << numZeros)-1;
  int N = ((bits>>(numZeros+1)) & explicitBitsMask) | (1 << numZeros);

  bitPosition += numZeros*2+1;

  //encoding of N plus encoding of value could go over 32 bits (if the encoded value has 23 significant bits or so)
  bits = getBits(buffer, bitPosition);
  bits &= (1 << (N-1)) - 1; //only N-1 bits
  bits |= 1 << (N-1); //implicit 1
  bitPosition += N-1;

  if (rowStart) {
    //the first value is stored as the value+1 (for when rows start with 0)
    previousIndex = bits - 1;
  } else {
    //values are stored as difference from the previous index
    previousIndex = bits + previousIndex;
  }
}


#endif // DECODENEXTINDEX_HPP
