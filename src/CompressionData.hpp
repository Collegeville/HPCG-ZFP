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

#ifndef COMPRESSIONDATA_HPP
#define COMPRESSIONDATA_HPP

#include <cmath>
#include <cstdint>
#include "Geometry.hpp"

//error bound modes
constexpr const int PW_REL = 0b0001;
constexpr const int ABS = 0b0010;
//constexpr const int REL = 0b0100;
constexpr const int AND_MODES = 0b1000;
constexpr const int ABS_AND_PW_REL = AND_MODES | ABS | PW_REL;
constexpr const int ABS_OR_PW_REL = ABS | PW_REL;

/*!
  The error bound mode as per SZ
*/
constexpr const int errorMode = ABS;

//error bounds as per SZ
constexpr const double ABS_ERROR_BOUND = 1e-10;
constexpr const double PWREL_BOUND_RATIO = 1.0e-10;

//number of bits each value takes up when compressed
constexpr const int COMPRESSED_BITS = 2;
constexpr const int VALUES_PER_COMPRESSED_BYTE = 8/COMPRESSED_BITS;
constexpr const int COMPRESSED_VALUE_MASK = (1 << COMPRESSED_BITS) - 1;

//Compression modes
constexpr const int UNCOMPRESSED = 0x00;
constexpr const int NEIGHBOR = 0x01;
constexpr const int OVER_NEIGHBOR = 0x02; //two elements before

//default values.
constexpr const int INITIAL_NEIGHBOR = -1;
constexpr const int INITIAL_OVER_NEIGHBOR = 1;


struct CompressionData_STRUCT {
  uint8_t * fValsCompressed;
  double * fValsUncompressed;

  uint8_t * bValsCompressed;
  double * bValsUncompressed;

  double * diagonalValues;

  uint8_t ** mtxIndsL;
  local_int_t localNumberOfRows;
};
typedef CompressionData_STRUCT CompressionData;


inline void DeleteCompressionData(CompressionData & data) {
  if (data.fValsCompressed) {
    delete [] data.fValsCompressed;
    data.fValsCompressed = 0;
  }
  if (data.fValsUncompressed) {
    delete [] data.fValsUncompressed;
    data.fValsUncompressed = 0;
  }
  if (data.bValsCompressed) {
    delete [] data.bValsCompressed;
    data.bValsCompressed = 0;
  }
  if (data.bValsUncompressed) {
    delete [] data.bValsUncompressed;
    data.bValsUncompressed = 0;
  }
  //diagonalValues re-uses allocation that is deallocated with matrix deallocation
  //if(data.diagonalValues) {
  //  delete [] data.diagonalValues;
  //  data.diagonalValues = 0;
  //}
  for (int i = 0; i<data.localNumberOfRows; i++) {
    delete [] data.mtxIndsL[i];
  }
  delete [] data.mtxIndsL;
}


#endif // COMPRESSIONDATA_HPP
