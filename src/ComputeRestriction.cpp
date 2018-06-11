
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

  Vector & Axf = *A.mgData->Axf;
  Vector & rc = *A.mgData->rc;
  local_int_t * f2c = A.mgData->f2cOperator;
  local_int_t nc = A.mgData->rc->localLength;


  #ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
  #endif
  for (local_int_t block = 0; block<nc/BLOCK_SIZE; block++){
    double rcBlock[BLOCK_SIZE];
    double rfBlock[BLOCK_SIZE];
    double AxfBlock[BLOCK_SIZE];

    local_int_t i = block*BLOCK_SIZE;
    for (local_int_t j = 0; j < BLOCK_SIZE; j++) {
      local_int_t dest = f2c[i+j];
      DecodeBlock(rf, dest/BLOCK_SIZE, rfBlock);
      DecodeBlock(Axf, dest/BLOCK_SIZE, AxfBlock);
      rcBlock[j] = rfBlock[dest%BLOCK_SIZE] - AxfBlock[dest%BLOCK_SIZE];
    }
    EncodeBlock(rc, block, rcBlock);
  }
  if (nc%BLOCK_SIZE != 0) {
    double rcBlock[BLOCK_SIZE];
    double rfBlock[BLOCK_SIZE];
    double AxfBlock[BLOCK_SIZE];

    local_int_t block = nc/BLOCK_SIZE;
    local_int_t i = block*BLOCK_SIZE;

    for (local_int_t j = 0; j < nc%BLOCK_SIZE; j++) {
      local_int_t dest = f2c[i+j];
      DecodeBlock(rf, dest/BLOCK_SIZE, rfBlock);
      DecodeBlock(Axf, dest/BLOCK_SIZE, AxfBlock);
      rcBlock[j] = rfBlock[dest%BLOCK_SIZE] - AxfBlock[dest%BLOCK_SIZE];
    }
    EncodeBlock(rc, nc/BLOCK_SIZE, rcBlock);
  }



  return 0;
}
