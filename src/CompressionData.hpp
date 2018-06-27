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

/*!
 @file CompressionData.hpp

 HPCG data structure for optimization data
 */

#ifndef COMPRESSIONDATA_HPP
#define COMPRESSIONDATA_HPP

struct CompressionData_Struct {
  float ** matrixValues;
  float ** matrixDiagonal;
  local_int_t numRows;

  uint8_t ** mtxIndL;
};
typedef struct CompressionData_Struct CompressionData;

inline void DeleteCompressionData(CompressionData & data) {
  if (data.matrixDiagonal) {
    delete [] data.matrixDiagonal;
    data.matrixDiagonal = NULL;
  }
    if (data.matrixValues) {
    for (local_int_t i = 0; i< data.numRows; ++i) {
      delete [] data.matrixValues[i];
    }
    delete [] data.matrixValues;
    data.matrixValues = NULL;
  }
  if (data.mtxIndL) {
    for (int i = 0; i<data.numRows; i++) {
      delete [] data.mtxIndL[i];
    }
    delete [] data.mtxIndL;
    data.mtxIndL = NULL;
  }
}

#endif //COMPRESSIONDATA_HPP
