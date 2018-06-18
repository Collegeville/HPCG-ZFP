
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

//Circular dependancy between Vector.hpp and DecodeBlock.hpp
struct SparseMatrix_STRUCT;
typedef struct SparseMatrix_STRUCT SparseMatrix;

void DecodeNextIndex_Forward(const SparseMatrix & mat, local_int_t & id, local_int_t & uncompressedCount, local_int_t & value);
void DecodeNextIndex_Backward(const SparseMatrix & mat, local_int_t & id, local_int_t & uncompressedCount, local_int_t & value);

#endif // DECODEBLOCK_HPP
