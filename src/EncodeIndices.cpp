
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

#include "EncodeIndices.hpp"

#include "CompressionData.hpp"
#include "SparseMatrix.hpp"

/*!
  Compresses a block into the given Vector.
  Note that the last block may be partial.

  @param[inout] vect   The vector to read from, must have optimization data.
  @param[in] blockid   The index of the block to read
  @param[in] block    The memory with the new content of the block.
*/
void EncodeIndices(SparseMatrix & mat, const bool forward) {
  unsigned char * compressed;
  local_int_t * uncompressedArray;
  int direction;
  if (forward) {
    compressed = ((CompressionData*)mat.optimizationData)->forwardCompressed;
    uncompressedArray = ((CompressionData*)mat.optimizationData)->forwardUncompressed;
    direction = 1;
  } else {
    compressed = ((CompressionData*)mat.optimizationData)->backwardCompressed;
    uncompressedArray = ((CompressionData*)mat.optimizationData)->backwardUncompressed;
    direction = -1;
  }

  //start compressed section as all 0's
  for (int i = 0; i < ceil(mat.localNumberOfNonzeros/8.0); i++) {
    compressed[i] = 0;
  }

  local_int_t ** indices = mat.mtxIndL;
  const char * nnzInRow = mat.nonzerosInRow;

  local_int_t index = 0;
  //indices should be non-negative
  local_int_t previousElement = -5;
  local_int_t uncompressedCount = 0;

  if (forward) {
    for (local_int_t row = 0; row < mat.localNumberOfRows; row++) {
      for (int j = 0; j < nnzInRow[row]; j++){
        local_int_t value = indices[row][j];
        if (value == previousElement+1) {
          compressed[index/8] |= INCREMENT << index%8;
        } else {
          compressed[index/8] |= UNCOMPRESSED << index%8;
          uncompressedArray[uncompressedCount] = value;
          uncompressedCount++;
        }
        previousElement = value;
        index++;
      }
    }
  } else {
    for (local_int_t row = mat.localNumberOfRows; row-- > 0; row) {
      for (int j = nnzInRow[row]; j-- > 0; ){
        local_int_t value = indices[row][j];
        if (value == previousElement-1) {
          compressed[index/8] |= INCREMENT << index%8;
        } else {
          compressed[index/8] |= UNCOMPRESSED << index%8;
          uncompressedArray[uncompressedCount] = value;
          uncompressedCount++;
        }
        previousElement = value;
        index++;
      }
    }
  }
}
