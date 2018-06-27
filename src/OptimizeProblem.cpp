
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
 @file OptimizeProblem.cpp

 HPCG routine
 */

#include "OptimizeProblem.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>

double optimizationAllocation = 0.0;

// computes floor(log_2(v))+1 where v is positive
inline int numBits(local_int_t v) {
  //copied from https://graphics.stanford.edu/~seander/bithacks.html
  const unsigned int b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000};
  const unsigned int S[] = {1, 2, 4, 8, 16};

  register unsigned int r = 0; // result of log2(v) will go here
  for (int i = 4; i >= 0; i--) {
    if (v & b[i]) {
      v >>= S[i];
      r |= S[i];
    }
  }
  //post loop r == floor(log2(v))
  return r+1;
}

//write bits to the given buffer
inline void putBits(uint8_t * buffer, uint32_t val, int bitPos) {
  // based off https://stackoverflow.com/a/5723250/6353993
  *(uint64_t*)(buffer + bitPos/8) |= val << (bitPos%8);
}

/*!
  helper function to create compressionData

  @param[inout] mat   The matrix to optimize

  @return the number of bytes used
*/
int CreateCompressedArray(SparseMatrix & mat){
  assert(sizeof(local_int_t) == 4);

  //temp[j] = inds[j]-inds[j-1] for j > 0; temp[0] = temp[0]-1
  //then compress with Elias delta coding

  uint8_t temp[213]; //maximum needed space, will allocate the exact amount for each row

  uint8_t ** mtxIndL = new uint8_t*[mat.localNumberOfRows];

  double totalMemory = sizeof(*mtxIndL)*mat.localNumberOfRows;
  for (local_int_t i = 0; i < mat.localNumberOfRows; i++) {
    int nonzerosInRow = mat.nonzerosInRow[i];
    local_int_t * inds = mat.mtxIndL[i];
    double * vals = mat.matrixValues[i];
    double *& diagPtr = mat.matrixDiagonal[i];

    // need row indices in ascending order
    if (!std::is_sorted(inds, inds+nonzerosInRow)) {
      for (int j = 0; j<nonzerosInRow-1; j++) {
        local_int_t * nextElement = std::min_element(inds+j, inds+nonzerosInRow);
        int index = nextElement - inds;
        std::swap(vals[j], vals[index]);
        std::swap(inds[j], inds[index]);
        if (vals+index == diagPtr) {
          diagPtr = vals+j;
        } else if (vals+j == diagPtr) {
          diagPtr = vals+index;
        }
      }
    }


    std::fill_n(temp, 142, 0 );
    int position = 0;
    //number of bits for inds[0]
    int N = numBits(inds[0]+1);
    assert(N <= 29); // otherwise we can't always decode with a single read from the buffer
    position += N-1;
    //put last N-1 bits of delta plus seperator bit
    unsigned int explicitBitsMask = (1 << (N-1))-1;
    putBits(temp, (((inds[0]+1)&explicitBitsMask) << 1) | 1, position);
    position += N;

    for (int j = 1; j < nonzerosInRow; j++) {
      local_int_t delta = inds[j]-inds[j-1];

      N = numBits(delta);
      assert(N <= 29); // otherwise we can't always decode with a single read from the buffer
      //gamma code N
      position += N-1;
      //put last N-1 bits of delta
      explicitBitsMask = (1 << (N-1))-1;
      putBits(temp, ((delta&explicitBitsMask) << 1) | 1, position);
      position += N;
    }

    int bytesUsed = ceil(position/8.0);
    mtxIndL[i] = new uint8_t[bytesUsed];
    std::copy(temp, temp+bytesUsed, mtxIndL[i]);
    totalMemory += bytesUsed;
  }

  float ** values = new float*[mat.localNumberOfRows];
  float ** diag = new float*[mat.localNumberOfRows];

  local_int_t elts = 0;
  for (int i = 0; i < mat.localNumberOfRows; i++) {
    values[i] = new float[mat.nonzerosInRow[i]];
    for (int j = 0; j < mat.nonzerosInRow[i]; j++){
      values[i][j] = (float)mat.matrixValues[i][j];
    }
    diag[i] = values[i] + (mat.matrixDiagonal[i]-mat.matrixValues[i]);
    totalMemory += sizeof(float) * mat.nonzerosInRow[i];
  }
  CompressionData * data = new CompressionData();
  data->matrixValues = values;
  data->matrixDiagonal = diag;
  data->numRows = mat.localNumberOfRows;
  data->mtxIndL = mtxIndL;

  mat.optimizationData = data;
  return sizeof(data)
        + mat.localNumberOfRows*sizeof(*values)*2 //arrays of pointers
        + totalMemory;
}

/*!
  Optimizes the data structures used for CG iteration to increase the
  performance of the benchmark version of the preconditioned CG algorithm.

  @param[inout] A      The known system matrix, also contains the MG hierarchy in attributes Ac and mgData.
  @param[inout] data   The data structure with all necessary CG vectors preallocated
  @param[inout] b      The known right hand side vector
  @param[inout] x      The solution vector to be computed in future CG iteration
  @param[inout] xexact The exact solution vector

  @return returns 0 upon success and non-zero otherwise

  @see GenerateGeometry
  @see GenerateProblem
*/
int OptimizeProblem(SparseMatrix & A, CGData & data, Vector & b, Vector & x, Vector & xexact) {

  // This function can be used to completely transform any part of the data structures.
  // Right now it does nothing, so compiling with a check for unused variables results in complaints

#if defined(HPCG_USE_MULTICOLORING)
  const local_int_t nrow = A.localNumberOfRows;
  std::vector<local_int_t> colors(nrow, nrow); // value `nrow' means `uninitialized'; initialized colors go from 0 to nrow-1
  int totalColors = 1;
  colors[0] = 0; // first point gets color 0

  // Finds colors in a greedy (a likely non-optimal) fashion.

  for (local_int_t i=1; i < nrow; ++i) {
    if (colors[i] == nrow) { // if color not assigned
      std::vector<int> assigned(totalColors, 0);
      int currentlyAssigned = 0;
      const local_int_t * const currentColIndices = A.mtxIndL[i];
      const int currentNumberOfNonzeros = A.nonzerosInRow[i];

      for (int j=0; j< currentNumberOfNonzeros; j++) { // scan neighbors
        local_int_t curCol = currentColIndices[j];
        if (curCol < i) { // if this point has an assigned color (points beyond `i' are unassigned)
          if (assigned[colors[curCol]] == 0)
            currentlyAssigned += 1;
          assigned[colors[curCol]] = 1; // this color has been used before by `curCol' point
        } // else // could take advantage of indices being sorted
      }

      if (currentlyAssigned < totalColors) { // if there is at least one color left to use
        for (int j=0; j < totalColors; ++j)  // try all current colors
          if (assigned[j] == 0) { // if no neighbor with this color
            colors[i] = j;
            break;
          }
      } else {
        if (colors[i] == nrow) {
          colors[i] = totalColors;
          totalColors += 1;
        }
      }
    }
  }

  std::vector<local_int_t> counters(totalColors);
  for (local_int_t i=0; i<nrow; ++i)
    counters[colors[i]]++;

  local_int_t old, old0;
  for (int i=1; i < totalColors; ++i) {
    old0 = counters[i];
    counters[i] = counters[i-1] + old;
    old = old0;
  }optimization_allocation
  counters[0] = 0;

  // translate `colors' into a permutation
  for (local_int_t i=0; i<nrow; ++i) // for each color `c'
    colors[i] = counters[colors[i]]++;
#endif

  optimizationAllocation = 0;

  SparseMatrix * Anext = &A;
  while (Anext) {
    optimizationAllocation += CreateCompressedArray(*Anext);
    Anext = Anext->Ac;
  }

  return 0;
}

// Helper function (see OptimizeProblem.hpp for details)
double OptimizeProblemMemoryUse(const SparseMatrix & A) {

  return optimizationAllocation;

}
