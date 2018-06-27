
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

#ifndef ENCODEVALUES_HPP
#define ENCODEVALUES_HPP

//Circular dependacy between EncodeValues.hpp and SparseMatrix.hpp
#include "Geometry.hpp"
struct SparseMatrix_STRUCT;
typedef struct SparseMatrix_STRUCT SparseMatrix;

void EncodeValues(SparseMatrix & mat, const bool forward);

#endif // ENCODEVALUES_HPP
