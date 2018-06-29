
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

#ifndef COMPUTEWAXPBY_HPP
#define COMPUTEWAXPBY_HPP

#include "ComputeWAXPBY.hpp"
#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>
#include "Vector.hpp"

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
template<class T, class U, class V>
int ComputeWAXPBY(const local_int_t n, const double alpha, const Vector<T> & x,
    const double beta, const Vector<U> & y, Vector<V> & w, bool & isOptimized) {

  assert(x.localLength>=n); // Test vector lengths
  assert(y.localLength>=n);

	const T * const xv = (T*)x.optimizationData;
	const U * const yv = (U*)y.optimizationData;
	V * const wv = (V*)w.optimizationData;
  if (alpha==1.0) {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
#endif
    for (local_int_t i=0; i<n; i++) wv[i] = (double)xv[i] + beta * (double)yv[i];
  } else if (beta==1.0) {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
#endif
    for (local_int_t i=0; i<n; i++) wv[i] = alpha * (double)xv[i] + (double)yv[i];
  } else  {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
#endif
    for (local_int_t i=0; i<n; i++) wv[i] = alpha * (double)xv[i] + beta * (double)yv[i];
  }

  return 0;
}


#endif // COMPUTEWAXPBY_HPP
