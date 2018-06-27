
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
 @file ComputeSYMGS.cpp

 HPCG routine
 */

#include "ComputeSYMGS.hpp"
#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#endif
#include <cassert>
#include "DecodeNextIndex.hpp"
#include "DecodeNextValue.hpp"

/*!
  Routine to one step of symmetrix Gauss-Seidel:

  Assumption about the structure of matrix A:
  - Each row 'i' of the matrix has nonzero diagonal value whose address is matrixDiagonal[i]
  - Entries in row 'i' are ordered such that:
       - lower triangular terms are stored before the diagonal element.
       - upper triangular terms are stored after the diagonal element.
       - No other assumptions are made about entry ordering.

  Symmetric Gauss-Seidel notes:
  - We use the input vector x as the RHS and start with an initial guess for y of all zeros.
  - We perform one forward sweep.  Since y is initially zero we can ignore the upper triangular terms of A.
  - We then perform one back sweep.
       - For simplicity we include the diagonal contribution in the for-j loop, then correct the sum after

  @param[in]  A the known system matrix
  @param[in]  x the input vector
  @param[out] y On exit contains the result of one symmetric GS sweep with x as the RHS.

  @return returns 0 upon success and non-zero otherwise

  @warning Early versions of this kernel (Version 1.1 and earlier) had the r and x arguments in reverse order, and out of sync with other kernels.

  @see ComputeSYMGS_ref
*/
int ComputeSYMGS( const SparseMatrix & A, const Vector & r, Vector & x) {

  assert(x.localLength==A.localNumberOfColumns); // Make sure x contain space for halo values

#ifndef HPCG_NO_MPI
  ExchangeHalo(A,x);
#endif

  const local_int_t nrow = A.localNumberOfRows;
  double * matrixDiagonal = ((CompressionData*)A.optimizationData)->diagonalValues;
  const double * const rv = r.values;
  double * const xv = x.values;

  uint8_t ** indices = ((CompressionData*)A.optimizationData)->mtxIndL;

  local_int_t indexId = 0;
  local_int_t uCount = 0;
  double curVal = INITIAL_NEIGHBOR;
  double prevVal = INITIAL_OVER_NEIGHBOR;

  for (local_int_t i = 0; i < nrow; i++) {
    const uint8_t * const cur_inds = indices[i];
    int bitPosition = 0;
    local_int_t curCol = -1;
    const int currentNumberOfNonzeros = A.nonzerosInRow[i];
    const double currentDiagonal = matrixDiagonal[i]; // Current diagonal value
    double sum = rv[i]; // RHS value

    for (int j = 0; j < currentNumberOfNonzeros; j++) {
      DecodeNextIndex(cur_inds, bitPosition, curCol);
      DecodeNextValue(A, indexId, uCount, curVal, prevVal, true);
      sum -= curVal*xv[curCol];
    }
    sum += xv[i]*currentDiagonal; // Remove diagonal contribution from previous loop

    xv[i] = sum/currentDiagonal;

  }

  // Now the back sweep.

  indexId = 0;
  uCount = 0;
  curVal = INITIAL_NEIGHBOR;
  prevVal = INITIAL_OVER_NEIGHBOR;

  for (local_int_t i=nrow-1; i>=0; i--) {
    const uint8_t * const cur_inds = indices[i];
    int bitPosition = 0;
    local_int_t curCol = -1;
    const int currentNumberOfNonzeros = A.nonzerosInRow[i];
    const double currentDiagonal = matrixDiagonal[i]; // Current diagonal value
    double sum = rv[i]; // RHS value

    for (int j = 0; j < currentNumberOfNonzeros; j++) {
      DecodeNextIndex(cur_inds, bitPosition, curCol);
      DecodeNextValue(A, indexId, uCount, curVal, prevVal, false);
      sum -= curVal*xv[curCol];
    }
    sum += xv[i]*currentDiagonal; // Remove diagonal contribution from previous loop

    xv[i] = sum/currentDiagonal;
  }

  return 0;
}
