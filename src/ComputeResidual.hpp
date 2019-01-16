
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

#ifndef COMPUTERESIDUAL_HPP
#define COMPUTERESIDUAL_HPP

#ifndef HPCG_NO_MPI
#include <mpi.h>
#endif
#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#include "Vector.hpp"

#ifdef HPCG_DETAILED_DEBUG
#include <fstream>
#include "hpcg.hpp"
#endif

#include <cmath>  // needed for fabs
#include "ComputeResidual.hpp"
#ifdef HPCG_DETAILED_DEBUG
#include <iostream>
#endif

/*!
  Routine to compute the inf-norm difference between two vectors where:

  @param[in]  n        number of vector elements (local to this processor)
  @param[in]  v1, v2   input vectors
  @param[out] residual pointer to scalar value; on exit, will contain result: inf-norm difference

  @return Returns zero on success and a non-zero value otherwise.
*/
template<class T, class U>
int ComputeResidual(const local_int_t n, const Vector<T> & v1, const Vector<U> & v2, double & residual) {

  T * v1v = (T*)v1.optimizationData;
  U * v2v = (U*)v2.optimizationData;
  double local_residual = 0.0;

#ifndef HPCG_NO_OPENMP
  #pragma omp parallel default(none) shared(local_residual, v1v, v2v)
  {
    double threadlocal_residual = 0.0;
    #pragma omp for
    for (local_int_t i=0; i<n; i++) {
      double diff = std::fabs(v1v[i] - v2v[i]);
      threadlocal_residual += diff*diff;
    }
    #pragma omp critical
    {
      local_residual += threadlocal_residual;
    }
  }
#else // No threading
  for (local_int_t i=0; i<n; i++) {
    double diff = std::fabs(v1v[i] - v2v[i]);
    local_residual += diff*diff;
#ifdef HPCG_DETAILED_DEBUG
    HPCG_fout << " Computed, exact, diff = " << v1v[i] << " " << v2v[i] << " " << diff << std::endl;
#endif
  }
#endif

#ifndef HPCG_NO_MPI
  // Use MPI's reduce function to collect all partial sums
  double global_residual = 0;
  MPI_Allreduce(&local_residual, &global_residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  residual = global_residual;
#else
  residual = local_residual;
#endif

  residual = std::sqrt(residual);

  return 0;
}

#endif // COMPUTERESIDUAL_HPP
