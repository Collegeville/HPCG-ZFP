
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

#include "EncodeValues.hpp"

#include "CompressionData.hpp"
#include "SparseMatrix.hpp"

inline bool withinTolerance(const double value, const double error) {
  if (errorMode & AND_MODES) {
    int withinTol = true;
    if (errorMode & PW_REL) {
      withinTol &= error <= fabs(value)*PWREL_BOUND_RATIO;
    }
    if(errorMode & ABS) {
      withinTol &= error <= ABS_ERROR_BOUND;
    }
    return withinTol;
  } else {
    int withinTol = false;
    if (errorMode & PW_REL) {
      withinTol |= error <= fabs(value)*PWREL_BOUND_RATIO;
    }
    if(errorMode & ABS) {
      withinTol |= error <= ABS_ERROR_BOUND;
    }
    return withinTol;
  }
}

inline void compressValue(local_int_t & index, const float value, float & neighbor,
                        uint8_t * compressed, float * uncompressedArray, local_int_t & uncompressedCount) {
  double neighborErr = fabs(value-neighbor);

  if (withinTolerance(value, neighborErr)) {
    compressed[index/VALUES_PER_COMPRESSED_BYTE] |= NEIGHBOR << (index%VALUES_PER_COMPRESSED_BYTE * VAL_COMPRESSED_BITS);
  } else {
    compressed[index/VALUES_PER_COMPRESSED_BYTE] |= UNCOMPRESSED << (index%VALUES_PER_COMPRESSED_BYTE * VAL_COMPRESSED_BITS);
    uncompressedArray[uncompressedCount] = value;
    uncompressedCount++;
    neighbor = value;
  }
  index++;
}

/*!
  Compresses the matrix's values using SZ compression.
  Three modes are used: uncompressed, neighbor and neighbor's neighbor.

  @param[inout] vect   The vector to read from, must have optimization data.
  @param[in] blockid   The index of the block to read
  @param[in] block    The memory with the new content of the block.
*/
void EncodeValues(SparseMatrix & mat, const bool forward) {
  uint8_t * compressed;
  float * uncompressedArray;
  if (forward) {
    compressed = ((CompressionData*)mat.optimizationData)->fValsCompressed;
    uncompressedArray = ((CompressionData*)mat.optimizationData)->fValsUncompressed;
  } else {
    compressed = ((CompressionData*)mat.optimizationData)->bValsCompressed;
    uncompressedArray = ((CompressionData*)mat.optimizationData)->bValsUncompressed;
  }

  //start compressed section as all 0's
  for (int i = 0; i < ceil(mat.localNumberOfNonzeros/(double)VALUES_PER_COMPRESSED_BYTE); i++) {
    compressed[i] = 0;
  }

  double ** values = mat.matrixValues;
  const char * nnzInRow = mat.nonzerosInRow;

  local_int_t index = 0;
  float neighbor = INITIAL_NEIGHBOR;
  local_int_t uncompressedCount = 0;

  if (forward) {
    for (local_int_t row = 0; row < mat.localNumberOfRows; row++) {
      for (int j = 0; j < nnzInRow[row]; j++){
        compressValue(index, values[row][j], neighbor,
                      compressed, uncompressedArray, uncompressedCount);
      }
    }
  } else {
    for (local_int_t row = mat.localNumberOfRows-1; row >= 0; row--) {
      for (int j = 0; j < nnzInRow[row]; j++){
        compressValue(index, values[row][j], neighbor,
                      compressed, uncompressedArray, uncompressedCount);
      }
    }
  }
}
