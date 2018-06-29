
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

#ifndef GENERATEPROBLEM_HPP
#define GENERATEPROBLEM_HPP
#include "SparseMatrix.hpp"
#include "Vector.hpp"

void GenerateProblem(SparseMatrix & A, Vector<b_type> * b, Vector<x_type> * x, Vector<x_type> * xexact);
#endif // GENERATEPROBLEM_HPP
