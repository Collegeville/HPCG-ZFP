
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
struct Vector_STRUCT;
typedef struct Vector_STRUCT Vector;

int DecodeBlock(const Vector & stream, local_int_t blockid, double * block);
int PartialDecodeBlock(const Vector & vect, const local_int_t blockid, const local_int_t blockLength, double * block);

#endif // DECODEBLOCK_HPP
