
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
#include "zfp.h"

//Circular dependancy between SparseMatrix.hpp and DecodeBlock.hpp
struct SparseMatrix_STRUCT;
typedef struct SparseMatrix_STRUCT SparseMatrix;

int DecodeBlock(const SparseMatrix & mat, local_int_t blockid, double * block);
int DecodeBlock(zfp_stream * zfp, local_int_t arrayLength, local_int_t blockid, double * block);
int DecodeFullBlock(const SparseMatrix & mat, const local_int_t blockid, double * block);
int DecodeFullBlock(zfp_stream * zfp, const local_int_t blockid, double * block);

#endif // DECODEBLOCK_HPP
