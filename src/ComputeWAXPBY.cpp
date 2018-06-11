
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
 @file ComputeWAXPBY.cpp

 HPCG routine
 */

#include "ComputeWAXPBY.hpp"
#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>

/*!
  Routine to compute the update of a vector with the sum of two
  scaled vectors where: w = alpha*x + beta*y

  This routine calls the reference WAXPBY implementation by default, but
  can be replaced by a custom, optimized routine suited for
  the target system.

  @param[in] n the number of vector elements (on this processor)
  @param[in] alpha, beta the scalars applied to x and y respectively.
  @param[in] x, y the input vectors
  @param[out] w the output vector
  @param[out] isOptimized should be set to false if this routine uses the reference implementation (is not optimized); otherwise leave it unchanged

  @return returns 0 upon success and non-zero otherwise

  @see ComputeWAXPBY_ref
*/
int ComputeWAXPBY(const local_int_t n, const double alpha, const Vector & x,
    const double beta, const Vector & y, Vector & w, bool & isOptimized) {

  assert(x.localLength>=n); // Test vector lengths
  assert(y.localLength>=n);
  assert(w.localLength>=n);


  if (alpha==1.0) {
    #ifndef HPCG_NO_OPENMP
        #pragma omp parallel for
    #endif
    for (local_int_t block = 0; block<n/BLOCK_SIZE; block++) {
      double xBlock[BLOCK_SIZE];
      double yBlock[BLOCK_SIZE];
      double wBlock[BLOCK_SIZE];

      DecodeBlock(x, block, xBlock);
      DecodeBlock(y, block, yBlock);
      for (local_int_t j = 0; j < BLOCK_SIZE; j++) {
        wBlock[j] = xBlock[j] + beta * yBlock[j];
      }
      EncodeBlock(w, block, wBlock);
    }
    if (n%BLOCK_SIZE != 0) {
      double xBlock[BLOCK_SIZE];
      double yBlock[BLOCK_SIZE];
      double wBlock[BLOCK_SIZE];

      DecodeBlock(x, n/BLOCK_SIZE, xBlock);
      DecodeBlock(y, n/BLOCK_SIZE, yBlock);
      for (local_int_t j = 0; j < n%BLOCK_SIZE; j++) {
        wBlock[j] = xBlock[j] + beta * yBlock[j];
      }
      EncodeBlock(w, n/BLOCK_SIZE, wBlock);
    }
  } else if (beta==1.0) {
    #ifndef HPCG_NO_OPENMP
      #pragma omp parallel for
    #endif
    for (local_int_t block = 0; block<n/BLOCK_SIZE; block++) {
      double xBlock[BLOCK_SIZE];
      double yBlock[BLOCK_SIZE];
      double wBlock[BLOCK_SIZE];

      DecodeBlock(x, block, xBlock);
      DecodeBlock(y, block, yBlock);
      for (local_int_t j = 0; j < BLOCK_SIZE; j++) {
        wBlock[j] = alpha * xBlock[j] + yBlock[j];
      }
      EncodeBlock(w, block, wBlock);
    }
    if (n%BLOCK_SIZE != 0) {
      double xBlock[BLOCK_SIZE];
      double yBlock[BLOCK_SIZE];
      double wBlock[BLOCK_SIZE];

      DecodeBlock(x, n/BLOCK_SIZE, xBlock);
      DecodeBlock(y, n/BLOCK_SIZE, yBlock);
      for (local_int_t j = 0; j < n%BLOCK_SIZE; j++) {
        wBlock[j] = alpha * xBlock[j] + yBlock[j];
      }
      EncodeBlock(w, n/BLOCK_SIZE, wBlock);
    }
  } else {
    #ifndef HPCG_NO_OPENMP
      #pragma omp parallel for
    #endif
    for (local_int_t block = 0; block<n/BLOCK_SIZE; block++) {
      double xBlock[BLOCK_SIZE];
      double yBlock[BLOCK_SIZE];
      double wBlock[BLOCK_SIZE];

      DecodeBlock(x, block, xBlock);
      DecodeBlock(y, block, yBlock);
      for (local_int_t j = 0; j < BLOCK_SIZE; j++) {
        wBlock[j] = alpha * xBlock[j] + beta * yBlock[j];
      }
      EncodeBlock(w, block, wBlock);
    }
    if (n%BLOCK_SIZE != 0) {
      double xBlock[BLOCK_SIZE];
      double yBlock[BLOCK_SIZE];
      double wBlock[BLOCK_SIZE];

      DecodeBlock(x, n/BLOCK_SIZE, xBlock);
      DecodeBlock(y, n/BLOCK_SIZE, yBlock);
      for (local_int_t j = 0; j < n%BLOCK_SIZE; j++) {
        wBlock[j] = alpha * xBlock[j] + beta * yBlock[j];
      }
      EncodeBlock(w, n/BLOCK_SIZE, wBlock);
    }
  }

    return 0;
}
