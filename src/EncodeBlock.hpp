
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
struct Vector_STRUCT;
typedef struct Vector_STRUCT Vector;

int EncodeBlock(Vector & vect, const local_int_t blockid, const double * block);

#endif // ENCODEBLOCK_HPP
