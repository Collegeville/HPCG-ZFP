
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

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#endif
#include "ComputeSYMGS.hpp"
#include <cassert>

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
  local_int_t * rowStarts = ((CompressionData*)A.optimizationData)->rowStarts;
  local_int_t * diagonalIndices = ((CompressionData*)A.optimizationData)->diagonalIndices;  // An array of pointers to the diagonal entries A.matrixValues
  double currentBlock[4];
  double diagCurrentBlock[4];
  const double * const rv = r.values;
  double * const xv = x.values;

  for (local_int_t i=0; i< nrow; i++) {
    local_int_t matrixRowStart = rowStarts[i];
    const local_int_t * const currentColIndices = A.mtxIndL[i];
    const int currentNumberOfNonzeros = A.nonzerosInRow[i];
    DecodeBlock(A, diagonalIndices[i]/4, diagCurrentBlock);
    const double currentDiagonal = diagCurrentBlock[(diagonalIndices[i])%4]; // Current diagonal value
    double sum = rv[i]; // RHS value

    for (int j=0; j< currentNumberOfNonzeros; j++) {
      if ((matrixRowStart+j)%4 == 0) {
        DecodeBlock(A, (matrixRowStart+j)/4, currentBlock);
      }
      local_int_t curCol = currentColIndices[j];
      sum -= currentBlock[(matrixRowStart+j)%4] * xv[curCol];
    }
    sum += xv[i]*currentDiagonal; // Remove diagonal contribution from previous loop

    xv[i] = sum/currentDiagonal;

  }

  // Now the back sweep.

  for (local_int_t i=nrow-1; i>=0; i--) {
    local_int_t matrixRowStart = rowStarts[i];
    const local_int_t * const currentColIndices = A.mtxIndL[i];
    const int currentNumberOfNonzeros = A.nonzerosInRow[i];
    DecodeBlock(A, diagonalIndices[i]/4, diagCurrentBlock);
    const double currentDiagonal = diagCurrentBlock[(diagonalIndices[i])%4]; // Current diagonal value
    double sum = rv[i]; // RHS value

    for (int j = currentNumberOfNonzeros-1; j >= 0; j--) {
      if ((matrixRowStart+j)%4 == 3) {
        DecodeBlock(A, (matrixRowStart+j)/4, currentBlock);
      }
      local_int_t curCol = currentColIndices[j];
      sum -= currentBlock[(matrixRowStart+j)%4]*xv[curCol];
    }
    sum += xv[i]*currentDiagonal; // Remove diagonal contribution from previous loop

    xv[i] = sum/currentDiagonal;
  }

  return 0;

}
