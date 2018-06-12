
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
 @file ComputeSPMV.cpp

 HPCG routine
 */

#include "ComputeSPMV.hpp"

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>

/*!
  Routine to compute sparse matrix vector product y = Ax where:
  Precondition: First call exchange_externals to get off-processor values of x

  This routine calls the reference SpMV implementation by default, but
  can be replaced by a custom, optimized routine suited for
  the target system.

  @param[in]  A the known system matrix
  @param[in]  x the known vector
  @param[out] y the On exit contains the result: Ax.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSPMV_ref
*/
int ComputeSPMV( const SparseMatrix & A, Vector & x, Vector & y) {

  assert(x.localLength>=A.localNumberOfColumns); // Test vector lengths
  assert(y.localLength>=A.localNumberOfRows);

#ifndef HPCG_NO_MPI
    ExchangeHalo(A,x);
#endif
  const local_int_t nrow = A.localNumberOfRows;

  #ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
  #endif
  for (local_int_t block=0; block < nrow/BLOCK_SIZE; block++)  {
    double yBlock[BLOCK_SIZE];
    double xBlock[BLOCK_SIZE];
    local_int_t lastDecodeBlock = -1;
    local_int_t lastDecodePosition;

    local_int_t i = block*BLOCK_SIZE;
    for (local_int_t k = 0; k < BLOCK_SIZE; k++) {
      double sum = 0.0;
      const double * const cur_vals = A.matrixValues[i+k];
      const local_int_t * const cur_inds = A.mtxIndL[i+k];
      const int cur_nnz = A.nonzerosInRow[i+k];

      for (int j=0; j < cur_nnz; j++) {
        local_int_t src = cur_inds[j];
        local_int_t srcBlock = src/BLOCK_SIZE;
        local_int_t srcPos = src%BLOCK_SIZE;

        if (srcBlock != lastDecodeBlock) {
          PartialDecodeBlock(x, srcBlock, srcPos+1, xBlock);
          lastDecodeBlock = srcBlock;
          lastDecodePosition = srcPos+1;
        } else if(srcPos >= lastDecodePosition) {
          ResumePartialDecodeBlock(x, srcBlock, srcPos+1, lastDecodePosition, xBlock);
          lastDecodePosition = srcPos+1;
        }
        sum += cur_vals[j]*xBlock[srcPos];
      }
      yBlock[k] = sum;
    }
    EncodeBlock(y, block, yBlock);
  }

  if (nrow%BLOCK_SIZE) {
    double yBlock[BLOCK_SIZE];
    double xBlock[BLOCK_SIZE];
    local_int_t lastDecodeBlock = -1;
    local_int_t lastDecodePosition = 0;

    local_int_t block = nrow/BLOCK_SIZE;
    local_int_t i = block*BLOCK_SIZE;
    for (local_int_t k = 0; k < nrow%BLOCK_SIZE; k++) {
      double sum = 0.0;
      const double * const cur_vals = A.matrixValues[i+k];
      const local_int_t * const cur_inds = A.mtxIndL[i+k];
      const int cur_nnz = A.nonzerosInRow[i+k];

      for (int j=0; j< cur_nnz; j++) {
        local_int_t src = cur_inds[j];
        local_int_t srcBlock = src/BLOCK_SIZE;
        local_int_t srcPos = src%BLOCK_SIZE;
        if (srcBlock != lastDecodeBlock) {
          PartialDecodeBlock(x, srcBlock, srcPos+1, xBlock);
          lastDecodeBlock = srcBlock;
          lastDecodePosition = srcPos+1;
        } else if(srcPos >= lastDecodePosition) {
          ResumePartialDecodeBlock(x, srcBlock, srcPos+1, lastDecodePosition, xBlock);
          lastDecodePosition = srcPos+1;
        }
        sum += cur_vals[j]*xBlock[srcPos];
      }
      yBlock[k] = sum;
    }
    EncodeBlock(y, block, yBlock);
  }
  return 0;
}
