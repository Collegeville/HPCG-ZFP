
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
#include "zfp.h"

//Circular dependancy between SparseMatrix.hpp and EncodeBlock.hpp
struct SparseMatrix_STRUCT;
typedef struct SparseMatrix_STRUCT SparseMatrix;

int EncodeBlock(SparseMatrix & mat, const local_int_t blockid, double * block);
int EncodeBlock(zfp_stream * zfp, local_int_t arrayLength, local_int_t blockid, double * block);
int EncodeFullBlock(SparseMatrix & mat, const local_int_t blockid, double * block);
int EncodeFullBlock(zfp_stream * zfp, local_int_t blockid, double * block);

#endif // ENCODEBLOCK_HPP
