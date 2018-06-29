
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

#ifndef CHECKPROBLEM_HPP
#define CHECKPROBLEM_HPP
#include "SparseMatrix.hpp"
#include "Vector.hpp"

void CheckProblem(const SparseMatrix & A, Vector<b_type> * b, Vector<x_type> * x, Vector<x_type> * xexact);
#endif // CHECKPROBLEM_HPP
