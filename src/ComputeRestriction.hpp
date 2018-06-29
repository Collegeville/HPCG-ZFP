
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

#ifndef COMPUTERESTRICTION_HPP
#define COMPUTERESTRICTION_HPP

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include "SparseMatrix.hpp"
#include "Vector.hpp"

/*!
  Routine to compute the coarse residual vector.

  @param[inout]  A - Sparse matrix object containing pointers to mgData->Axf, the fine grid matrix-vector product and mgData->rc the coarse residual vector.
  @param[in]    rf - Fine grid RHS.


  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
*/
template<class T>
int ComputeRestriction(const SparseMatrix & A, const Vector<T> & rf) {

  double * Axfv = (double*)A.mgData->Axf->optimizationData;
  T * rfv = (T*)rf.optimizationData;
  double * rcv = (double*)A.mgData->rc->optimizationData;
  local_int_t * f2c = A.mgData->f2cOperator;
  local_int_t nc = A.mgData->rc->localLength;

#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
  for (local_int_t i=0; i<nc; ++i) rcv[i] = (double)rfv[f2c[i]] - (double)Axfv[f2c[i]];

  return 0;
}
#endif // COMPUTERESTRICTION_HPP
