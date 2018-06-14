
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

#include "zfparray1.h"

struct OptimizedData_STRUCT {
  zfp::array1d * matrixValues;
  zfp::array1d::pointer * matrixDiagonal;
};
typedef struct OptimizedData_STRUCT OptimizedData;

inline void DeleteOptimizedData(OptimizedData & data) {
  delete [] data.matrixDiagonal;
  data.matrixDiagonal = 0;
  delete [] data.matrixValues;
  data.matrixValues = 0;
}


#endif //OPTIMIZEDDATA_HPP
