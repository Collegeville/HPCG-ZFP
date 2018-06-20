
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

#ifndef ENCODEBLOCK_HPP
#define ENCODEBLOCK_HPP

#include "Geometry.hpp"

//Circular dependancy between Vector.hpp and DecodeBlock.hpp
struct SparseMatrix_STRUCT;
typedef struct SparseMatrix_STRUCT SparseMatrix;

void EncodeIndices(SparseMatrix & mat, const bool forward);

#endif // ENCODEBLOCK_HPP
