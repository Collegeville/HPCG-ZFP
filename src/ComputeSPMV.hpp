
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

#ifndef COMPUTESPMV_HPP
#define COMPUTESPMV_HPP

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>
#include "SparseMatrix.hpp"
#include "Vector.hpp"
#include "DecodeNextIndex.hpp"


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
template<class T, class U>
int ComputeSPMV( const SparseMatrix & A, Vector<T> & x, Vector<U> & y) {


  assert(x.localLength>=A.localNumberOfColumns); // Test vector lengths
  assert(y.localLength>=A.localNumberOfRows);

#ifndef HPCG_NO_MPI
    ExchangeHalo(A,x);
#endif
  const T * const xv = (T*)x.optimizationData;
  U * const yv = (U*)y.optimizationData;
  const local_int_t nrow = A.localNumberOfRows;
  uint8_t ** indices = ((CompressionData*)A.optimizationData)->mtxIndL;
  float ** matrixValues = ((CompressionData*)A.optimizationData)->matrixValues;
#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for (local_int_t i=0; i< nrow; i++)  {
    double sum = 0.0;
    const float * cur_vals = matrixValues[i];
    const uint8_t * const cur_inds = indices[i];
    int bitPosition = 0;
    local_int_t cur_col = -1;
    const int cur_nnz = A.nonzerosInRow[i];

    for (int j=0; j< cur_nnz; j++) {
      DecodeNextIndex(cur_inds, bitPosition, cur_col);
      sum += cur_vals[j]*(double)xv[cur_col];
    }
    yv[i] = sum;
  }
  return 0;
}

#endif  // COMPUTESPMV_HPP
