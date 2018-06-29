
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
 @file Vector.hpp

 HPCG data structures for dense vectors
 */

#ifndef VECTOR_HPP
#define VECTOR_HPP
#include <cassert>
#include <cstdlib>
#include "Geometry.hpp"
#include <type_traits>

template <class T>
struct Vector_STRUCT {
  local_int_t localLength;  //!< length of local portion of the vector
  double * values;          //!< array of values
  /*!
   This is for storing optimized data structures created in OptimizeProblem and
   used inside optimized ComputeSPMV().
   */
  void * optimizationData;

};

template <class T>
using Vector = Vector_STRUCT<T>;

/*!
  Initializes input vector.

  @param[in] v
  @param[in] localLength Length of local portion of input vector
 */
template <class T>
inline void InitializeVector(Vector<T> & v, local_int_t localLength) {
  v.localLength = localLength;
  v.values = new double[localLength];
  v.optimizationData = 0;
  return;
}

/*!
  Fill the input vector with zero values.

  @param[inout] v - On entrance v is initialized, on exit all its values are zero.
 */
template <class T>
inline void ZeroVector(Vector<T> & v) {
  local_int_t localLength = v.localLength;
  if (v.optimizationData && std::is_same<T, float>::value) {
    float * vv = (float*)v.optimizationData;
    for (int i=0; i<localLength; ++i) vv[i] = 0.0;
  } else {
    double * vv = v.values;
    for (int i=0; i<localLength; ++i) vv[i] = 0.0;
  }
  return;
}
/*!
  Multiply (scale) a specific vector entry by a given value.

  @param[inout] v Vector to be modified
  @param[in] index Local index of entry to scale
  @param[in] value Value to scale by
 */
template <class T>
inline void ScaleVectorValue(Vector<T> & v, local_int_t index, double value) {
  assert(index>=0 && index < v.localLength);
  if (v.optimizationData && std::is_same<T, float>::value) {
    float * vv = (float*)v.optimizationData;
    vv[index] *= value;
  } else {
    double * vv = v.values;
    vv[index] *= value;
  }
  return;
}
/*!
  Fill the input vector with pseudo-random values.

  @param[in] v
 */
template <class T>
inline void FillRandomVector(Vector<T> & v) {
  local_int_t localLength = v.localLength;
  if (v.optimizationData && std::is_same<T, float>::value) {
    float * vv = (float*)v.optimizationData;
    for (int i=0; i<localLength; ++i) vv[i] = rand() / (double)(RAND_MAX) + 1.0;
  } else {
    double * vv = v.values;
    for (int i=0; i<localLength; ++i) vv[i] = rand() / (double)(RAND_MAX) + 1.0;
  }
  return;
}
/*!
  Copy input vector to output vector.

  @param[in] v Input vector
  @param[in] w Output vector
 */
template <class T, class U>
inline void CopyVector(const Vector<T> & v, Vector<U> & w) {
  local_int_t localLength = v.localLength;
  assert(w.localLength >= localLength);
  if (v.optimizationData) {
    T * vv = (T*)v.optimizationData;
    U * wv = (U*)w.optimizationData;
    for (int i=0; i<localLength; ++i) wv[i] = vv[i];
  } else {
    double * vv = v.values;
    double * wv = w.values;
    for (int i=0; i<localLength; ++i) wv[i] = vv[i];
  }
  return;
}


/*!
  Deallocates the members of the data structure of the known system matrix provided they are not 0.

  @param[in] A the known system matrix
 */
template <class T>
inline void DeleteVector(Vector<T> & v) {

  delete [] v.values;
  v.localLength = 0;
  return;
}

#endif // VECTOR_HPP
