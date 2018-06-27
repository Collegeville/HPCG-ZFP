
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

//only the bottom 57 bits are garenteed to be accurate
inline uint64_t getBits(const uint8_t * buffer, int bitPos) {
  // based off https://stackoverflow.com/a/5723250/6353993
  uint64_t bits = *(uint64_t*)(buffer + bitPos/8);
  return bits >> bitPos%8;
}

// count number of zeros in the front of the read location
inline int countZeros(uint32_t bits) {
  //copied from https://graphics.stanford.edu/~seander/bithacks.html
  if (bits & 0x1) {
    return 0;
  } else {
    int count = 1;
    if ((bits & 0xFFFF) == 0) {
      bits >>= 16;
      count += 16;
    }
    if ((bits & 0xFF) == 0) {
      bits >>= 8;
      count += 8;
    }
    if ((bits & 0xF) == 0) {
      bits >>= 4;
      count += 4;
    }
    if ((bits & 0x3) == 0) {
      bits >>= 2;
      count += 2;
    }
    return count - (bits&0x1);
  }
}

/*!
  Decompresses the next index of the matrix's row
*/
inline void DecodeNextIndex(const uint8_t * buffer, int & bitPosition, local_int_t & previousIndex) {
  uint64_t bits = getBits(buffer, bitPosition);

  //length of N
  int Nsub1 = countZeros(bits);
  bits >>= Nsub1+1;

  bits &= (1 << Nsub1) - 1; //only N-1 bits
  bits |= 1 << Nsub1; //implicit 1
  bitPosition += 2*Nsub1 + 1;

  //values are stored as difference from the previous index
  previousIndex = bits + previousIndex;
}


#endif // DECODENEXTINDEX_HPP
