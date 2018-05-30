
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
  assert(x.optimizationData);
  assert(y.optimizationData);

  const zfp::array1d & xv = *(zfp::array1d*)x.optimizationData;
  const zfp::array1d & yv = *(zfp::array1d*)y.optimizationData;
  zfp::array1d & wv = *(zfp::array1d*)w.optimizationData;
if (alpha==1.0) {
  for (local_int_t i=0; i<n; i++) wv[i] = xv[i] + beta * yv[i];
} else if (beta==1.0) {
  for (local_int_t i=0; i<n; i++) wv[i] = alpha * xv[i] + yv[i];
} else  {
  for (local_int_t i=0; i<n; i++) wv[i] = alpha * xv[i] + beta * yv[i];
}

return 0;
}
