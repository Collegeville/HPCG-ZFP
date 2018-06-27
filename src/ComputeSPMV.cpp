
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
#include "DecodeNextIndex.hpp"
#include "DecodeNextValue.hpp"


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
  const double * const xv = x.values;
  double * const yv = y.values;
  const local_int_t nrow = A.localNumberOfRows;
  uint8_t ** indices = ((CompressionData*)A.optimizationData)->mtxIndL;

  local_int_t indexId = 0;
  local_int_t uCount = 0;
  double curVal = INITIAL_NEIGHBOR;
  double prevVal = INITIAL_OVER_NEIGHBOR;

//#ifndef HPCG_NO_OPENMP
//  #pragma omp parallel for
//#endif
  for (local_int_t i=0; i< nrow; i++)  {
    double sum = 0.0;
    const uint8_t * const cur_inds = indices[i];
    int bitPosition = 0;
    local_int_t cur_col = -1;
    const int cur_nnz = A.nonzerosInRow[i];

    for (int j = 0; j < cur_nnz; j++){
      DecodeNextIndex(cur_inds, bitPosition, cur_col);
      DecodeNextValue(A, indexId, uCount, curVal, prevVal, true);
      sum += curVal*xv[cur_col];
    }
    yv[i] = sum;
  }
  return 0;
}
