
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
#include "ComputeSYMGS_ref.hpp"
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
  assert(r.optimizationData);
  assert(x.optimizationData);

#ifndef HPCG_NO_MPI
  ExchangeHalo(A, x);
#endif

  const local_int_t nrow = A.localNumberOfRows;
  double ** matrixDiagonal = A.matrixDiagonal;  // An array of pointers to the diagonal entries A.matrixValues


  double rBlock [BLOCK_SIZE];
  double xBlock[BLOCK_SIZE];

  double tempBlock [BLOCK_SIZE];//cache to most recently decoded block
  int lastDecodeBlock = -1;
  int lastDecodePosition;

  for (local_int_t block=0; block < nrow/BLOCK_SIZE; block++)  {
    PartialDecodeBlock(r, block, BLOCK_SIZE, rBlock);
    PartialDecodeBlock(x, block, BLOCK_SIZE,xBlock);

    local_int_t i = block*BLOCK_SIZE;
    for (local_int_t k = 0; k < BLOCK_SIZE; k++) {
      const double * const currentValues = A.matrixValues[i+k];
      const local_int_t * const currentColIndices = A.mtxIndL[i+k];
      const int currentNumberOfNonzeros = A.nonzerosInRow[i+k];
      const double  currentDiagonal = matrixDiagonal[i+k][0]; // Current diagonal value
      double sum = rBlock[k]; // RHS value

      for (int j=0; j< currentNumberOfNonzeros; j++) {
        local_int_t curCol = currentColIndices[j];
        local_int_t curBlock = curCol/BLOCK_SIZE;
        local_int_t curPos = curCol%BLOCK_SIZE;
        if (curBlock != block) {
          if (curBlock != lastDecodeBlock) {
            PartialDecodeBlock(x, curBlock, curPos+1, tempBlock);
            lastDecodeBlock = curBlock;
            lastDecodePosition = curPos+1;
          } else if(curPos >= lastDecodePosition) {
            ResumePartialDecodeBlock(x, curBlock, curPos+1, lastDecodePosition, tempBlock);
            lastDecodePosition = curPos+1;
          }
          sum -= currentValues[j] * tempBlock[curPos];
        } else if (curPos != k) {// Exclude diagonal contribution
          sum -= currentValues[j] * xBlock[curPos];
        }
      }

      xBlock[k] = sum/currentDiagonal;
    }
    EncodeBlock(x, block, xBlock);
  }
  if (nrow%BLOCK_SIZE) {
    local_int_t block = nrow/BLOCK_SIZE;
    local_int_t i = block*BLOCK_SIZE;
    DecodeBlock(r, block, rBlock);
    DecodeBlock(x, block, xBlock);
    for (local_int_t k = 0; k < nrow%BLOCK_SIZE; k++) {
      const double * const currentValues = A.matrixValues[i+k];
      const local_int_t * const currentColIndices = A.mtxIndL[i+k];
      const int currentNumberOfNonzeros = A.nonzerosInRow[i+k];
      const double  currentDiagonal = matrixDiagonal[i+k][0]; // Current diagonal value
      double sum = rBlock[k]; // RHS value

      for (int j=0; j< currentNumberOfNonzeros; j++) {
        local_int_t curCol = currentColIndices[j];
        local_int_t curBlock = curCol/BLOCK_SIZE;
        local_int_t curPos = curCol%BLOCK_SIZE;

        if (curBlock != block) {
          if (curBlock != lastDecodeBlock) {
            PartialDecodeBlock(x, curBlock, curPos+1, tempBlock);
            lastDecodeBlock = curBlock;
            lastDecodePosition = curPos+1;
          } else if(curPos >= lastDecodePosition) {
            ResumePartialDecodeBlock(x, curBlock, curPos+1, lastDecodePosition, tempBlock);
            lastDecodePosition = curPos+1;
          }
          sum -= currentValues[j] * tempBlock[curPos];
        } else if (curPos != k) {// Exclude diagonal contribution
          sum -= currentValues[j] * xBlock[curPos];
        }
      }

      xBlock[k] = sum/currentDiagonal;
    }
//    EncodeBlock(x, block, xBlock);
//  }

  // Now the back sweep.

//  if (nrow%BLOCK_SIZE) {
//    local_int_t block = nrow/BLOCK_SIZE;
//    local_int_t i = block*BLOCK_SIZE;
//    DecodeBlock(r, block, rBlock);
//    DecodeBlock(x, block, xBlock);
    for (local_int_t k = nrow%BLOCK_SIZE-1; k >=0 ; k--) {
      const double * const currentValues = A.matrixValues[i+k];
      const local_int_t * const currentColIndices = A.mtxIndL[i+k];
      const int currentNumberOfNonzeros = A.nonzerosInRow[i+k];
      const double  currentDiagonal = matrixDiagonal[i+k][0]; // Current diagonal value
      double sum = rBlock[k]; // RHS value

      for (int j=0; j< currentNumberOfNonzeros; j++) {
        local_int_t curCol = currentColIndices[j];
        local_int_t curBlock = curCol/BLOCK_SIZE;
        local_int_t curPos = curCol%BLOCK_SIZE;

        if (curBlock != block) {
          if (curBlock != lastDecodeBlock) {
            PartialDecodeBlock(x, curBlock, curPos+1, tempBlock);
            lastDecodeBlock = curBlock;
            lastDecodePosition = curPos+1;
          } else if(curPos >= lastDecodePosition) {
            ResumePartialDecodeBlock(x, curBlock, curPos+1, lastDecodePosition, tempBlock);
            lastDecodePosition = curPos+1;
          }
          sum -= currentValues[j] * tempBlock[curPos];
        } else if (curPos != k) {// Exclude diagonal contribution
          sum -= currentValues[j] * xBlock[curPos];
        }
      }

      xBlock[k] = sum/currentDiagonal;
    }
    EncodeBlock(x, block, xBlock);
  }
  for (local_int_t block=(nrow/BLOCK_SIZE) - 1; block >= 0; block--)  {
    PartialDecodeBlock(r, block, BLOCK_SIZE, rBlock);
    PartialDecodeBlock(x, block, BLOCK_SIZE, xBlock);

    local_int_t i = block*BLOCK_SIZE;
    for (local_int_t k = BLOCK_SIZE-1; k >=0; k--) {
      const double * const currentValues = A.matrixValues[i+k];
      const local_int_t * const currentColIndices = A.mtxIndL[i+k];
      const int currentNumberOfNonzeros = A.nonzerosInRow[i+k];
      const double  currentDiagonal = matrixDiagonal[i+k][0]; // Current diagonal value
      double sum = rBlock[k]; // RHS value

      for (int j=0; j< currentNumberOfNonzeros; j++) {
        local_int_t curCol = currentColIndices[j];
        local_int_t curBlock = curCol/BLOCK_SIZE;
        local_int_t curPos = curCol%BLOCK_SIZE;

        if (curBlock != block) {
          if (curBlock != lastDecodeBlock) {
            PartialDecodeBlock(x, curBlock, curPos+1, tempBlock);
            lastDecodeBlock = curBlock;
            lastDecodePosition = curPos+1;
          } else if(curPos >= lastDecodePosition) {
            ResumePartialDecodeBlock(x, curBlock, curPos+1, lastDecodePosition, tempBlock);
            lastDecodePosition = curPos+1;
          }
          sum -= currentValues[j] * tempBlock[curPos];
        } else if (curPos != k) {// Exclude diagonal contribution
          sum -= currentValues[j] * xBlock[curPos];
        }
      }

      xBlock[k] = sum/currentDiagonal;
    }
    EncodeBlock(x, block, xBlock);
  }

  return 0;

}
