
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
#include <cmath>
#include <cstdlib>
#include "Geometry.hpp"
#include "CompressionData.hpp"
#include "EncodeBlock.hpp"
#include "DecodeBlock.hpp"

struct Vector_STRUCT {
  local_int_t localLength;  //!< length of local portion of the vector
  double * values;          //!< array of values
  /*!
   This is for storing optimized data structures created in OptimizeProblem and
   used inside optimized ComputeSPMV().
   */
  void * optimizationData;

};
typedef struct Vector_STRUCT Vector;

/*!
  Initializes input vector.

  @param[in] v
  @param[in] localLength Length of local portion of input vector
 */
inline void InitializeVector(Vector & v, local_int_t localLength) {
  v.localLength = localLength;
  v.values = new double[localLength];
  v.optimizationData = 0;
  return;
}

/*!
  Fill the input vector with zero values.

  @param[inout] v - On entrance v is initialized, on exit all its values are zero.
 */
inline void ZeroVector(Vector & v) {
  local_int_t localLength = v.localLength;
  if (v.optimizationData) {
    double vBlock [BLOCK_SIZE] = {};
    local_int_t numBlocks = ceil(localLength/(double)BLOCK_SIZE);


    #ifndef HPCG_NO_OPENMP
      #pragma omp parallel for
    #endif
    for (local_int_t block = 0; block < numBlocks; block++){
      char * data = (char*)v.optimizationData+block*BLOCK_BYTES;
      
      data[0] = UNCOMPRESSED << 6 | NEIGHBOR << 4 | NEIGHBOR << 2 | NEIGHBOR;
      for(int i = 1; i < COMPRESSED_BYTES; i++) {
        data[i] = NEIGHBOR << 6 | NEIGHBOR << 4 | NEIGHBOR << 2 | NEIGHBOR;
      }
      ((double*)((void*)((char*)data+COMPRESSED_BYTES)))[0] = 0;
    }
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
inline void ScaleVectorValue(Vector & v, local_int_t index, double value) {
  assert(index>=0 && index < v.localLength);
  if (v.optimizationData) {
    local_int_t block = index/BLOCK_SIZE;
    double vBlock[BLOCK_SIZE];
    DecodeBlock(v, block, vBlock);
    vBlock[index%BLOCK_SIZE] *= value;
    EncodeBlock(v, block, vBlock);
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
inline void FillRandomVector(Vector & v) {
  local_int_t localLength = v.localLength;
  if (v.optimizationData) {

    local_int_t numBlocks = ceil(localLength/(double)BLOCK_SIZE);

    #ifndef HPCG_NO_OPENMP
      #pragma omp parallel for
    #endif
    for (local_int_t block = 0; block<numBlocks; block++){
      double vBlock [BLOCK_SIZE];
      for (local_int_t j = 0; j < BLOCK_SIZE; j++) {
        vBlock[j] = rand() / (double)(RAND_MAX) + 1.0;
      }
      EncodeBlock(v, block, vBlock);
    }
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
inline void CopyVector(const Vector & v, Vector & w) {
  local_int_t localLength = v.localLength;
  assert(w.localLength >= localLength);
  assert(!v.optimizationData == !w.optimizationData);
  if (v.optimizationData) {
    const char* vBuf = (char*)v.optimizationData;
    char* wBuf = (char*)w.optimizationData;
    local_int_t bytes = ceil(v.localLength/(double)BLOCK_SIZE)*BLOCK_BYTES;
    #ifndef HPCG_NO_OPENMP
      #pragma omp parallel for
    #endif
    for (int i=0; i<bytes; ++i) wBuf[i] = vBuf[i];

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
inline void DeleteVector(Vector & v) {

  delete [] v.values;
  v.localLength = 0;

  if (v.optimizationData) {
    //allocated with malloc/aligned_alloc
    free(v.optimizationData);
    v.optimizationData = NULL;
  }

  return;
}

#endif // VECTOR_HPP
