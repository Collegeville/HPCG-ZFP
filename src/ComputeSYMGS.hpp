
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

#ifndef COMPUTESYMGS_HPP
#define COMPUTESYMGS_HPP

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#endif
#include "ComputeSYMGS.hpp"
#include <cassert>
#include "DecodeNextIndex.hpp"
#include "DecodeNextValue.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"

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
template<class T, class U>
int ComputeSYMGS( const SparseMatrix & A, const Vector<T> & r, Vector<U> & x) {

  assert(x.localLength==A.localNumberOfColumns); // Make sure x contain space for halo values

#ifndef HPCG_NO_MPI
  ExchangeHalo(A,x);
#endif

  const local_int_t nrow = A.localNumberOfRows;
  float * matrixDiagonal = ((CompressionData*)A.optimizationData)->diagonalValues;
  const T * const rv = (T*)r.optimizationData;
  U * const xv = (U*)x.optimizationData;

  local_int_t index = 0;
  local_int_t valsUCount = 0;
  local_int_t indsUCount = 0;
  double curVal = INITIAL_NEIGHBOR;
  local_int_t curCol = INITIAL_NEIGHBOR;

  for (local_int_t i=0; i< nrow; i++) {
    const int currentNumberOfNonzeros = A.nonzerosInRow[i];
    const double  currentDiagonal = matrixDiagonal[i]; // Current diagonal value
    double sum = rv[i]; // RHS value

    for (int j=0; j< currentNumberOfNonzeros; j++) {
      DecodeNextIndex(A, index, indsUCount, curCol, true);
      DecodeNextValue(A, index, valsUCount, curVal, true);
      index++;
      sum -= curVal * (double)xv[curCol];
    }
    sum += xv[i]*currentDiagonal; // Remove diagonal contribution from previous loop

    xv[i] = sum/currentDiagonal;

  }

  // Now the back sweep.
  index = 0;
  valsUCount = 0;
  indsUCount = 0;
  curVal = INITIAL_NEIGHBOR;
  curCol = INITIAL_NEIGHBOR;
  for (local_int_t i = nrow-1; i >= 0; i--) {
    const int currentNumberOfNonzeros = A.nonzerosInRow[i];
    const double  currentDiagonal = matrixDiagonal[i]; // Current diagonal value
    double sum = rv[i]; // RHS value

    for (int j = 0; j< currentNumberOfNonzeros; j++) {
      DecodeNextIndex(A, index, indsUCount, curCol, false);
      DecodeNextValue(A, index, valsUCount, curVal, false);
      index++;
      sum -= curVal * (double)xv[curCol];
    }
    sum += xv[i]*currentDiagonal; // Remove diagonal contribution from previous loop

    xv[i] = sum/currentDiagonal;
  }

  return 0;

}


#endif // COMPUTESYMGS_HPP
