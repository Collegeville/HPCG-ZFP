
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

/*!
 @file ComputeRestriction.cpp

 HPCG routine
 */

#include "ComputeRestriction.hpp"

/*!
  Routine to compute the coarse residual vector.

  @param[inout]  A - Sparse matrix object containing pointers to mgData->Axf, the fine grid matrix-vector product and mgData->rc the coarse residual vector.
  @param[in]    rf - Fine grid RHS.


  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
*/
int ComputeRestriction(const SparseMatrix & A, const Vector & rf) {

  zfp::array1d & Axfv = *(zfp::array1d*)A.mgData->Axf->optimizationData;
  zfp::array1d & rfv = *(zfp::array1d*)rf.optimizationData;
  zfp::array1d & rcv = *(zfp::array1d*)A.mgData->rc->optimizationData;
  local_int_t * f2c = A.mgData->f2cOperator;
  local_int_t nc = A.mgData->rc->localLength;

  for (local_int_t i=0; i<nc; ++i) rcv[i] = rfv[f2c[i]] - Axfv[f2c[i]];

  return 0;
}
